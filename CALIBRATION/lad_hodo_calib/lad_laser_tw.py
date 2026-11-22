#!/usr/bin/env python3
"""Time walk from the LAD laser attenuation scan (lad_laser_histos.C outputs, one file per run).

For every PMT and run, the (TDC - photodiode) peak is measured in slices of FADC amplitude. All PMTs and
runs are then fit together:
    t[pmt, run](A) = c1[pmt, session] + c3[pmt] * (A/thr)^(-c2[pmt]) + delta[run]
delta[run] (one per run, delta[first run of each session] = 0) absorbs run-to-run shifts of the photodiode
time; c1 is separate per laser session (the PMT offsets change between sessions).

usage: lad_laser_tw.py out_prefix laser_RUN.root [laser_RUN.root ...] [--thr 120] [--amin 10] [--amax 900]
writes out_prefix.csv (per-PMT fit), out_prefix_points.csv, out_prefix.pdf
"""
import argparse
import csv
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_tw_fit import PLANES, NPAD, SIDES, peak_fit  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


SESSIONS = [(22350, 22360), (22820, 22840), (23240, 23260), (23790, 23830)]


def session(run):
    for i, (lo, hi) in enumerate(SESSIONS):
        if lo <= run <= hi:
            return i
    return len(SESSIONS)


def run_points(fn, amin, amax, nmin=150):
    """{(plane, side, paddle): [(A, t, et), ...]} for one laser run."""
    f = ROOT.TFile.Open(fn)
    out = {}
    for ipl, pl in enumerate(PLANES):
        for s in SIDES:
            for ip in range(1, NPAD[ipl] + 1):
                h = f.Get(f"{pl}/h_las_{pl}_{s}_{ip}")
                if not h or h.GetEntries() < 500:
                    continue
                px = h.ProjectionX("_px", 1, h.GetNbinsY())
                if px.GetMean() < amin or px.GetMean() > amax:
                    continue
                # amplitude slices: quantiles of this run's amplitude distribution
                qs = np.zeros(6)
                px.GetQuantiles(6, qs, np.array([0.02, 0.2, 0.4, 0.6, 0.8, 0.98]))
                pts = []
                for lo, hi in zip(qs[:-1], qs[1:]):
                    ax = h.GetXaxis()
                    b1, b2 = ax.FindBin(lo), ax.FindBin(hi)
                    if b2 <= b1:
                        continue
                    py = h.ProjectionY("_py", b1, b2)
                    if py.GetEntries() < nmin:
                        continue
                    r = peak_fit(py)
                    if r is None:
                        continue
                    pxs = h.ProjectionX("_pxs", 1, h.GetNbinsY())
                    pxs.GetXaxis().SetRange(b1, b2)
                    A = pxs.GetMean()
                    if amin <= A <= amax:
                        pts.append((A, r[0], max(r[1], 0.003)))
                if pts:
                    out[(pl, s, ip)] = pts
    f.Close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("files", nargs="+")
    ap.add_argument("--thr", type=float, default=120.0)
    ap.add_argument("--amin", type=float, default=10.0)
    ap.add_argument("--amax", type=float, default=900.0)
    ap.add_argument("--sys", type=float, default=0.01, help="ns added in quadrature to every point")
    ap.add_argument("--no-delta", action="store_true", help="fix delta[run] = 0 (photodiode time independent of attenuation)")
    a = ap.parse_args()

    runs, data = [], {}
    for fn in a.files:
        m = re.search(r"laser_(\d+)\.root", fn)
        run = int(m.group(1)) if m else len(runs)
        pts = run_points(fn, a.amin, a.amax)
        if not pts:
            print(f"{fn}: no points")
            continue
        runs.append(run)
        for k, v in pts.items():
            for A, t, e in v:
                data.setdefault(k, []).append((run, A, t, e))
        print(f"run {run}: {len(pts)} PMTs")
    runs = sorted(runs)
    ridx = {r: i for i, r in enumerate(runs)}
    pmts = sorted(k for k, v in data.items() if len({p[0] for p in v}) >= 3)
    print(f"{len(pmts)} PMTs with >= 3 runs, {len(runs)} runs")

    # arrays
    P, R, A, T, E = [], [], [], [], []
    for i, k in enumerate(pmts):
        for run, Ai, ti, ei in data[k]:
            P.append(i)
            R.append(ridx[run])
            A.append(Ai)
            T.append(ti)
            E.append(np.hypot(ei, a.sys))
    P, R = np.array(P, dtype=int), np.array(R, dtype=int)
    sess = sorted({session(r) for r in runs})
    sidx = {s_: i for i, s_ in enumerate(sess)}
    RS = np.array([sidx[session(runs[r])] for r in R], dtype=int)  # session index of each point
    ns = len(sess)
    first_of_session = {sidx[session(r)]: ridx[r] for r in reversed(runs)}
    free_runs = [] if a.no_delta else [i for i in range(len(runs)) if i not in first_of_session.values()]
    A, T, E = np.array(A, dtype=float), np.array(T, dtype=float), np.array(E, dtype=float)
    if len(P) == 0:
        sys.exit("no points")
    npm, nr = len(pmts), len(runs)
    thr = a.thr

    nf = len(free_runs)

    def unpack(x):
        c1 = x[:npm * ns].reshape(npm, ns)
        c2 = x[npm * ns:npm * ns + npm]
        c3 = x[npm * ns + npm:npm * ns + 2 * npm]
        d = np.zeros(nr)
        d[free_runs] = x[npm * ns + 2 * npm:]
        return c1, c2, c3, d

    def resid(x):
        c1, c2, c3, d = unpack(x)
        return (c1[P, RS] + c3[P] * np.power(A / thr, -c2[P]) + d[R] - T) / E

    # start: c1 = high-amplitude time of each PMT
    c1_0 = np.array([np.min([p[2] for p in data[k]]) for k in pmts])
    x0 = np.concatenate([np.repeat(c1_0, ns), np.full(npm, 0.3), np.full(npm, 2.0), np.zeros(nf)])
    lo = np.concatenate([np.full(npm * ns, -1e4), np.full(npm, 0.01), np.full(npm, 0.01), np.full(nf, -5)])
    hi = np.concatenate([np.full(npm * ns, 1e4), np.full(npm, 4.0), np.full(npm, 50.0), np.full(nf, 5)])
    res = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
    c1, c2, c3, d = unpack(res.x)
    r = res.fun
    print(f"chi2/ndf = {np.sum(r**2):.0f}/{len(r) - len(res.x)}")
    for run, dd in zip(runs, d):
        print(f"  run {run}: delta = {dd:+.3f} ns")

    with open(a.prefix + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["plane", "side", "paddle"] + [f"c1_s{SESSIONS[s_][0] if s_ < len(SESSIONS) else 'x'}" for s_ in sess]
                   + ["c2", "c3", "npts", "chi2", "amin", "amax"])
        for i, k in enumerate(pmts):
            sel = P == i
            w.writerow([*k, *[f"{c1[i, j]:.4f}" for j in range(ns)], f"{c2[i]:.4f}", f"{c3[i]:.4f}", int(sel.sum()),
                        f"{np.sum(r[sel]**2):.1f}", f"{A[sel].min():.1f}", f"{A[sel].max():.1f}"])
    with open(a.prefix + "_points.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["plane", "side", "paddle", "run", "amp", "t", "et", "t_minus_delta", "pull"])
        for j in range(len(P)):
            k = pmts[P[j]]
            w.writerow([*k, runs[R[j]], f"{A[j]:.2f}", f"{T[j]:.4f}", f"{E[j]:.4f}", f"{T[j]-d[R[j]]:.4f}", f"{r[j]:.2f}"])

    # plots: one page per plane/side, t - delta - c1 vs A with the fit
    c = ROOT.TCanvas("c", "", 1600, 1000)
    pdf = a.prefix + ".pdf"
    c.Print(pdf + "[")
    keep = []
    for ipl, pl in enumerate(PLANES):
        for s in SIDES:
            c.Clear()
            c.Divide(4, 3)
            for ip in range(1, NPAD[ipl] + 1):
                c.cd(ip)
                ROOT.gPad.SetLogx()
                if (pl, s, ip) not in pmts:
                    continue
                i = pmts.index((pl, s, ip))
                sel = P == i
                x = A[sel]
                y = T[sel] - d[R[sel]] - c1[i, RS[sel]]
                g = ROOT.TGraphErrors(int(sel.sum()), x.astype(float), y.astype(float), np.zeros(int(sel.sum())),
                                      E[sel].astype(float))
                g.SetTitle(f"{pl} {s} {ip};FADC amplitude (mV);t - c1 - #delta_{{run}} (ns)")
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.4)
                g.Draw("AP")
                f1 = ROOT.TF1(f"f{i}", f"[1]*pow(x/{thr},-[0])", max(x.min() * 0.9, 1), x.max() * 1.1)
                f1.SetParameters(c2[i], c3[i])
                f1.SetLineColor(ROOT.kRed)
                f1.Draw("same")
                t = ROOT.TLatex()
                t.SetNDC()
                t.SetTextSize(0.05)
                t.DrawLatex(0.4, 0.85, f"c2={c2[i]:.2f} c3={c3[i]:.2f}")
                keep += [g, f1, t]
            c.cd(12)
            t = ROOT.TLatex()
            t.SetTextSize(0.08)
            t.DrawLatexNDC(0.05, 0.5, f"laser TW: plane {pl} {s}")
            keep.append(t)
            c.Print(pdf)
    c.Print(pdf + "]")


if __name__ == "__main__":
    main()
