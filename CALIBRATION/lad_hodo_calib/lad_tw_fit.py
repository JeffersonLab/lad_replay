#!/usr/bin/env python3
"""Fit the LAD time walk from lad_tw_histos.C output.

For each PMT the (TDC - FADC time) peak is found in slices of FADC amplitude, and the peak
positions are fit with
    t(A) = c1 + c3 * (A / thr)^(-c2)          thr = TDC_threshold (mV)
LADlib applies tw(A) = c3 * [(A/thr)^(-c2) - (200/thr)^(-c2)] (c3 = 1 unless lladhodo_c3_* is set).
Two fits are made: c3 fixed to 1 (the original form) and c3 free.

usage: lad_tw_fit.py histos.root out_prefix [--thr 120] [--amin 15] [--amax 700] [--run N]
writes out_prefix.param (LADlib param file), out_prefix.csv (all fits) and out_prefix.pdf.
"""
import argparse
import csv
import math
import sys

import numpy as np
from scipy.optimize import curve_fit

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

PLANES = ["000", "001", "100", "101", "200", "REFBAR"]
NPAD = [11, 11, 11, 11, 11, 1]
SIDES = ["Top", "Btm"]

# amplitude slice edges (mV): fine where the walk changes fastest
EDGES = np.concatenate([np.arange(10, 40, 4), np.arange(40, 100, 8), np.arange(100, 300, 20), np.arange(300, 1001, 50)])


def peak_fit(h1):
    """Peak position of a 1D histogram: iterative Gaussian fit on [mu - 1.2 s, mu + 1.2 s]."""
    if h1.GetEntries() < 60:
        return None
    mu = h1.GetBinCenter(h1.GetMaximumBin())
    s = 0.6
    f = ROOT.TF1("g", "gaus", mu - 2, mu + 2)
    ok = None
    for _ in range(4):
        f.SetRange(mu - 1.2 * s, mu + 1.2 * s)
        f.SetParameters(h1.GetMaximum(), mu, s)
        r = h1.Fit(f, "QNRS")
        if int(r) != 0:
            break
        mu_new, s_new, emu = f.GetParameter(1), abs(f.GetParameter(2)), f.GetParError(1)
        if not (0.03 < s_new < 4.0) or abs(mu_new - mu) > 3:
            break
        mu, s, ok = mu_new, s_new, (mu_new, emu, s_new)
    return ok


def model(A, c1, c2, c3, thr):
    return c1 + c3 * np.power(A / thr, -c2)


def fit_pmt(h2, thr, amin, amax):
    pts = []
    ax = h2.GetXaxis()
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        if lo < amin or hi > amax:
            continue
        b1, b2 = ax.FindBin(lo + 1e-6), ax.FindBin(hi - 1e-6)
        py = h2.ProjectionY("_py", b1, b2)
        r = peak_fit(py)
        if r is None:
            continue
        # amplitude: mean of the slice
        px = h2.ProjectionX("_px", 1, h2.GetNbinsY())
        px.GetXaxis().SetRange(b1, b2)
        pts.append((px.GetMean(), r[0], max(r[1], 0.01), r[2], py.GetEntries()))
    if len(pts) < 5:
        return None, pts
    P = np.array(pts)
    A, t, et = P[:, 0], P[:, 1], P[:, 2]
    out = {}
    try:
        p, cov = curve_fit(lambda a, c1, c2: model(a, c1, c2, 1.0, thr), A, t, p0=[t[-1], 0.6], sigma=et,
                           absolute_sigma=True, bounds=([-1e4, 0.01], [1e4, 5]))
        res = (t - model(A, p[0], p[1], 1.0, thr)) / et
        out["fix"] = dict(c1=p[0], c2=p[1], c3=1.0, chi2=float(np.sum(res**2)), ndf=len(A) - 2)
    except Exception:
        out["fix"] = None
    try:
        p, cov = curve_fit(lambda a, c1, c2, c3: model(a, c1, c2, c3, thr), A, t, p0=[t[-1], 0.6, 1.0], sigma=et,
                           absolute_sigma=True, bounds=([-1e4, 0.01, 0.01], [1e4, 5, 50]))
        res = (t - model(A, *p, thr)) / et
        out["free"] = dict(c1=p[0], c2=p[1], c3=p[2], chi2=float(np.sum(res**2)), ndf=len(A) - 3)
    except Exception:
        out["free"] = None
    return out, pts


def walk(A, c2, c3, thr):
    """The correction LADlib applies (relative to 200 mV)."""
    return c3 * (np.power(A / thr, -c2) - np.power(200.0 / thr, -c2))


def read_fallback(fn):
    """{(plane, side, paddle): {"c1","c2","c3"}} from an existing LAD TW param file (c3 = 1 if absent)."""
    import re
    txt = open(fn).read()
    out = {}
    for par in ["c1", "c2", "c3"]:
        for s in SIDES:
            m = re.search(rf"^\s*lladhodo_{par}_{s}\s*=\s*(.*?)(?=^\s*;|^\s*\w+\s*=|\Z)", txt, re.S | re.M)
            if not m:
                continue
            vals = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(1))]
            for ip in range(11):
                for ipl, pl in enumerate(PLANES):
                    k = ip * 6 + ipl
                    if k < len(vals) and ip < NPAD[ipl]:
                        out.setdefault((pl, s, ip + 1), {"c1": 0.0, "c2": 0.0, "c3": 1.0})[par] = vals[k]
    return out


def write_param(fn, res, which, thr, run_label):
    with open(fn, "w") as f:
        f.write(f";LAD Hodoscopes time-walk parameters: {run_label}\n")
        f.write(f";written by CALIBRATION/lad_hodo_calib/lad_tw_fit.py ({which} fit)\n")
        f.write(";tw(A) = c3 * [(A/thr)^-c2 - (200/thr)^-c2], A = FADC pulse amplitude (mV); c1 is the TDC-FADC offset (not used)\n\n")
        f.write(f"lTDC_threshold={thr:.1f} ;units of mV\n\n")
        for par in ["c1", "c2", "c3"]:
            if par == "c3" and which == "fix":
                continue
            for s in SIDES:
                f.write(f";Param {par}-{s}\n;" + "".join(f"{p:>15s}" for p in PLANES) + "\n")
                f.write(f"lladhodo_{par}_{s} = ")
                for ip in range(11):
                    vals = []
                    for ipl, pl in enumerate(PLANES):
                        r = res.get((pl, s, ip + 1))
                        default = {"c1": 0.0, "c2": 0.0, "c3": 1.0}[par]
                        v = r[which][par] if (r and r.get(which)) else default
                        if ip >= NPAD[ipl]:
                            v = default
                        vals.append(f"{v:12.6f}")
                    f.write(("" if ip == 0 else " " * (len(f"lladhodo_{par}_{s} = "))) + ", ".join(vals) + "\n")
                f.write("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histos")
    ap.add_argument("prefix")
    ap.add_argument("--thr", type=float, default=120.0)
    ap.add_argument("--amin", type=float, default=15.0)
    ap.add_argument("--amax", type=float, default=700.0)
    ap.add_argument("--label", default="")
    ap.add_argument("--fallback", default="", help="TW param file for PMTs without a good fit (and REFBAR)")
    ap.add_argument("--min-entries", type=int, default=20000)
    a = ap.parse_args()

    fin = ROOT.TFile.Open(a.histos)
    res, allpts = {}, {}
    rows = []
    for ipl, pl in enumerate(PLANES):
        for s in SIDES:
            for ip in range(1, NPAD[ipl] + 1):
                h = fin.Get(f"{pl}/h_tw_{pl}_{s}_{ip}")
                if not h or h.GetEntries() < 500:
                    print(f"{pl} {s} {ip}: too few entries ({h.GetEntries() if h else 0})")
                    continue
                r, pts = fit_pmt(h, a.thr, a.amin, a.amax)
                allpts[(pl, s, ip)] = pts
                if r is None:
                    print(f"{pl} {s} {ip}: fit failed ({len(pts)} points)")
                    continue
                res[(pl, s, ip)] = r
                row = dict(plane=pl, side=s, paddle=ip, entries=int(h.GetEntries()), npts=len(pts))
                for w in ["fix", "free"]:
                    for k in ["c1", "c2", "c3", "chi2", "ndf"]:
                        row[f"{w}_{k}"] = r[w][k] if r.get(w) else float("nan")
                rows.append(row)

    if a.fallback:
        fb = read_fallback(a.fallback)
        for k, v in fb.items():
            r = res.get(k)
            row = next((x for x in rows if (x["plane"], x["side"], x["paddle"]) == k), None)
            bad = (r is None or k[0] == "REFBAR" or (row and row["entries"] < a.min_entries)
                   or not r.get("free") or r["free"]["c3"] >= 49.9)
            if bad:
                res[k] = {"fix": dict(v, chi2=float("nan"), ndf=0), "free": dict(v, chi2=float("nan"), ndf=0)}
                print(f"{k}: using fallback c2={v['c2']:.3f} c3={v['c3']:.3f}")
    with open(a.prefix + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    write_param(a.prefix + "_fix.param", res, "fix", a.thr, a.label)
    write_param(a.prefix + "_free.param", res, "free", a.thr, a.label)

    # summary
    for w in ["fix", "free"]:
        c2 = [r[w]["c2"] for r in res.values() if r.get(w)]
        c3 = [r[w]["c3"] for r in res.values() if r.get(w)]
        chi = [r[w]["chi2"] / max(r[w]["ndf"], 1) for r in res.values() if r.get(w) and r[w]["ndf"] > 0]
        print(f"{w}: n={len(c2)} c2 median {np.median(c2):.3f} [{np.min(c2):.3f},{np.max(c2):.3f}]  "
              f"c3 median {np.median(c3):.3f} [{np.min(c3):.3f},{np.max(c3):.3f}]  chi2/ndf median {np.median(chi):.1f}")

    # plots: one page per plane/side
    c = ROOT.TCanvas("c", "", 1600, 1000)
    c.Print(a.prefix + ".pdf[")
    keep = []
    for ipl, pl in enumerate(PLANES):
        for s in SIDES:
            c.Clear()
            c.Divide(4, 3)
            for ip in range(1, NPAD[ipl] + 1):
                c.cd(ip)
                ROOT.gPad.SetLogz()
                h = fin.Get(f"{pl}/h_tw_{pl}_{s}_{ip}")
                if not h:
                    continue
                h.GetXaxis().SetRangeUser(0, 800)
                h.Draw("colz")
                pts = allpts.get((pl, s, ip), [])
                if pts:
                    g = ROOT.TGraphErrors(len(pts))
                    for i, p in enumerate(pts):
                        g.SetPoint(i, p[0], p[1])
                        g.SetPointError(i, 0, p[2])
                    g.SetMarkerStyle(20)
                    g.SetMarkerSize(0.5)
                    g.Draw("P same")
                    keep.append(g)
                r = res.get((pl, s, ip))
                for w, col in [("fix", ROOT.kRed), ("free", ROOT.kBlack)]:
                    if r and r.get(w):
                        q = r[w]
                        f1 = ROOT.TF1(f"f_{pl}{s}{ip}{w}", f"[0]+[2]*pow(x/{a.thr},-[1])", a.amin, a.amax)
                        f1.SetParameters(q["c1"], q["c2"], q["c3"])
                        f1.SetLineColor(col)
                        f1.SetLineWidth(1)
                        f1.Draw("same")
                        keep.append(f1)
                if r and r.get("free"):
                    t = ROOT.TLatex()
                    t.SetNDC()
                    t.SetTextSize(0.05)
                    q = r["free"]
                    t.DrawLatex(0.35, 0.85, f"c2={q['c2']:.2f} c3={q['c3']:.2f} #chi^{{2}}/n={q['chi2']/max(q['ndf'],1):.1f}")
                    keep.append(t)
            c.cd(12)
            t = ROOT.TLatex()
            t.SetTextSize(0.08)
            t.DrawLatexNDC(0.05, 0.5, f"Plane {pl} {s}")
            t.DrawLatexNDC(0.05, 0.38, "red: c3=1, black: c3 free")
            keep.append(t)
            c.Print(a.prefix + ".pdf")
    c.Print(a.prefix + ".pdf]")


if __name__ == "__main__":
    main()
