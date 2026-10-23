#!/usr/bin/env python3
"""LAD top-bottom cable offsets and bar-to-bar offsets (LCoeff) from laser runs.

Input: lad_laser_histos.C outputs (TDC - photodiode vs FADC amplitude per PMT) and a time-walk param file
(lad_tw_fit.py output). The laser light enters each bar at its centre, so after the time-walk correction
    cableFit = (t_btm - t_top) / 2
    LCoeff   = T_ref - T_bar,   T = (t_top + t_btm)/2 - cableFit,   ref = plane 000 paddle 1
which is the LADlib convention (timec_top = t_top - tw + LCoeff, timec_btm = t_btm - tw - 2 cableFit + LCoeff).
Several runs are averaged (median per PMT); the run-to-run spread is reported.

usage: lad_laser_offsets.py out_prefix tw.param laser_RUN.root [...] [--vel 16.2635]
writes out_prefix.csv and out_prefix.param (velFit, cableFit, LCoeff; FADC versions copied from --fadc-from)
"""
import argparse
import csv
import re
import sys
import os

import numpy as np

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_tw_fit import PLANES, NPAD, SIDES, peak_fit  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
REF = ("000", 1)


def read_param_table(fn, name):
    """Read a LAD param array (11 rows x 6 planes) into {(plane, paddle): value}."""
    txt = open(fn).read()
    m = re.search(rf"^\s*l{name}\s*=\s*(.*?)(?=^\s*;|^\s*\w+\s*=|\Z)", txt, re.S | re.M)
    if not m:
        return None
    vals = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(1))]
    out = {}
    for ip in range(11):
        for ipl, pl in enumerate(PLANES):
            k = ip * 6 + ipl
            if k < len(vals):
                out[(pl, ip + 1)] = vals[k]
    return out


def read_thr(fn):
    m = re.search(r"lTDC_threshold\s*=\s*([-+\d.eE]+)", open(fn).read())
    return float(m.group(1)) if m else 120.0


def tw(A, c2, c3, thr):
    return c3 * (np.power(A / thr, -c2) - np.power(200.0 / thr, -c2))


def corrected_peak(h2, c2, c3, thr):
    """Peak of TDC - PD - tw(A) from the (A, t) histogram."""
    ax, ay = h2.GetXaxis(), h2.GetYaxis()
    h1 = ROOT.TH1D("h1c", "", 1200, ay.GetXmin() - 10, ay.GetXmax() + 10)
    for i in range(1, h2.GetNbinsX() + 1):
        A = ax.GetBinCenter(i)
        if A < 15:
            continue
        sh = tw(A, c2, c3, thr)
        for j in range(1, h2.GetNbinsY() + 1):
            w = h2.GetBinContent(i, j)
            if w > 0:
                h1.Fill(ay.GetBinCenter(j) - sh, w)
    r = peak_fit(h1)
    amp = h2.ProjectionX("_pxa", 1, h2.GetNbinsY()).GetMean()
    return (r[0], r[1], amp) if r else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("twparam")
    ap.add_argument("files", nargs="+")
    ap.add_argument("--vel", type=float, default=16.2635)
    ap.add_argument("--fadc-from", default="", help="param file to copy the *_FADC arrays from")
    ap.add_argument("--fallback", default="", help="param file for bars without laser data (cableFit, LCoeff)")
    ap.add_argument("--amp-range", default="60,700", help="use a PMT in a run only if its mean amplitude (mV) is in range")
    a = ap.parse_args()

    thr = read_thr(a.twparam)
    c2 = {s: read_param_table(a.twparam, f"ladhodo_c2_{s}") for s in SIDES}
    c3 = {s: (read_param_table(a.twparam, f"ladhodo_c3_{s}") or {}) for s in SIDES}

    alo, ahi = (float(x) for x in a.amp_range.split(","))
    per_run = {}  # (pl, side, pad) -> list of (run, t)
    runs = []
    for fn in a.files:
        m = re.search(r"laser_(\d+)\.root", fn)
        run = int(m.group(1)) if m else len(runs)
        f = ROOT.TFile.Open(fn)
        n = 0
        for ipl, pl in enumerate(PLANES):
            for s in SIDES:
                for ip in range(1, NPAD[ipl] + 1):
                    h = f.Get(f"{pl}/h_las_{pl}_{s}_{ip}")
                    if not h or h.GetEntries() < 1000:
                        continue
                    r = corrected_peak(h, c2[s].get((pl, ip), 0.0), c3[s].get((pl, ip), 1.0), thr)
                    if r and alo <= r[2] <= ahi:
                        per_run.setdefault((pl, s, ip), []).append((run, r[0], r[2]))
                        n += 1
        f.Close()
        runs.append(run)
        print(f"run {run}: {n} PMTs")

    # per-run bar quantities, then median over runs
    bars = {}
    for ipl, pl in enumerate(PLANES):
        for ip in range(1, NPAD[ipl] + 1):
            T = dict((r, t) for r, t, _ in per_run.get((pl, "Top", ip), []))
            B = dict((r, t) for r, t, _ in per_run.get((pl, "Btm", ip), []))
            common = sorted(set(T) & set(B))
            if not common:
                continue
            bars[(pl, ip)] = {r: (T[r], B[r]) for r in common}
    if REF not in bars:
        sys.exit("reference bar missing")
    rows = []
    out = {}
    for (pl, ip), d in bars.items():
        cab, lco, used = [], [], []
        for r, (t, b) in d.items():
            if r not in bars[REF]:
                continue
            c = 0.5 * (b - t)
            tr, br = bars[REF][r]
            cr = 0.5 * (br - tr)
            Tbar = 0.5 * (t + b) - c
            Tref = 0.5 * (tr + br) - cr
            cab.append(c)
            lco.append(Tref - Tbar)
            used.append(r)
        if not cab:
            continue
        out[(pl, ip)] = (float(np.median(cab)), float(np.median(lco)))
        rows.append(dict(plane=pl, paddle=ip, nruns=len(cab), cableFit=np.median(cab), cable_spread=np.std(cab),
                         LCoeff=np.median(lco), LCoeff_spread=np.std(lco),
                         per_run=" ".join(f"{r}:{x:.2f}" for r, x in zip(used, lco))))
    with open(a.prefix + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print(f"{len(rows)} bars; median run spread: cable {np.median([r['cable_spread'] for r in rows]):.3f} ns, "
          f"LCoeff {np.median([r['LCoeff_spread'] for r in rows]):.3f} ns")

    # param file
    fadc = {}
    if a.fadc_from:
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            fadc[nm] = read_param_table(a.fadc_from, nm)

    def table(name, getter, default=0.0):
        s = f"l{name} = "
        lines = []
        for ip in range(11):
            vals = []
            for ipl, pl in enumerate(PLANES):
                v = getter(pl, ip + 1) if ip < NPAD[ipl] else None
                vals.append(f"{(default if v is None else v):12.6f}")
            lines.append(", ".join(vals))
        return s + ("\n" + " " * len(s)).join(lines) + "\n"

    fb_cab = read_param_table(a.fallback, "ladhodo_cableFit") if a.fallback else {}
    fb_lco = read_param_table(a.fallback, "ladhodo_LCoeff") if a.fallback else {}
    missing = [k for k in fb_cab if k not in out and k[0] != "REFBAR"]
    if missing:
        print("no laser data, using " + os.path.basename(a.fallback) + " for: " + ", ".join(f"{p}-{b}" for p, b in missing))
    for k in fb_cab:
        if k not in out:
            out[k] = (fb_cab[k], fb_lco.get(k, 0.0))
    with open(a.prefix + ".param", "w") as f:
        f.write("; LAD hodoscope propagation velocity, top-bottom cable offsets and bar offsets\n")
        f.write(f"; from laser runs {', '.join(str(r) for r in runs)} (CALIBRATION/lad_hodo_calib/lad_laser_offsets.py)\n")
        f.write(f"; time walk: {os.path.basename(a.twparam)}; reference bar: plane {REF[0]} paddle {REF[1]}\n")
        f.write(";" + "".join(f"{p:>14s}" for p in PLANES) + "\n")
        f.write(table("ladhodo_velFit", lambda pl, ip: a.vel, a.vel) + "\n")
        f.write(table("ladhodo_cableFit", lambda pl, ip: out.get((pl, ip), (None, None))[0]) + "\n")
        f.write(table("ladhodo_LCoeff", lambda pl, ip: out.get((pl, ip), (None, None))[1]) + "\n")
        for nm, tab in fadc.items():
            if tab:
                f.write(table(nm, lambda pl, ip, tab=tab: tab.get((pl, ip))) + "\n")


if __name__ == "__main__":
    main()
