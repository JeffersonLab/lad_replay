#!/usr/bin/env python3
"""Back-plane LCoeff offsets for a run group from same-paddle front-back fast pairs, for periods without a photon
flash or a laser run (e.g. 22375-22535, 3-pass 22013-22120).

The fast-pair peak of t_back - t_front (h_dtb_amp_<back>_<paddle> from lad_edep_histos.C, replay constants) is
measured in the run group and in the reference group (runs 22698+). With constants K the corrected difference
changes by (L - c)_K,back - (L - c)_K,front, so the residual of the group with the input constants is
    r = (x_group - x_ref) + [(L-c)_in - (L-c)_ref]_back - [(L-c)_in - (L-c)_ref]_front
and the back bar's LCoeff is shifted by -r. Front planes, plane 200 and all cable offsets are unchanged, so this
fixes the front-back time difference per paddle, not the absolute time of either plane. Bars without a usable peak
(fit error > --max-err) take the plane median of the measured residuals.

usage: lad_pairs_period_offsets.py --group <histos.root>[,...] --dir <period dir> --ref <histos.root> --refdir p5_22698
           --in in.param --refparam ref.param --out out.param [--label "22375-22535"]
"""
import argparse
import os
import sys

import numpy as np
import ROOT
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_laser_offsets import read_param_table  # noqa: E402

ROOT.gROOT.SetBatch(True)
ALLPL = ["000", "001", "100", "101", "200", "REFBAR"]
PAIRS = [("000", "001"), ("100", "101")]


def peak(h, minn=150):
    """Fast-pair peak of the dt projection (flat level from 8.5-11.5 ns removed): (mean, error) or (None, None)."""
    if not h:
        return None, None
    px = h.ProjectionX("_px", 0, -1)
    x = np.array([px.GetBinCenter(i) for i in range(1, px.GetNbinsX() + 1)])
    y = np.array([px.GetBinContent(i) for i in range(1, px.GetNbinsX() + 1)])
    bkg = np.mean(y[(x > 8.5) & (x < 11.5)])
    y = y - bkg
    s = (x > -2) & (x < 5)
    if y[s].sum() < minn:
        return None, None
    mu = x[s][np.argmax(np.convolve(y, np.ones(3) / 3, "same")[s])]
    w = (x > mu - 0.8) & (x < mu + 0.8)
    try:
        p, c = curve_fit(lambda t, A, m, sg: A * np.exp(-0.5 * ((t - m) / sg) ** 2), x[w], y[w],
                         p0=[y[w].max(), mu, 0.4], sigma=np.sqrt(np.maximum(y[w] + bkg, 1)))
    except Exception:
        return None, None
    if not (-2 < p[1] < 5) or not np.isfinite(c[1, 1]):
        return None, None
    return float(p[1]), float(np.sqrt(c[1, 1]))


def summed(files, d, name):
    h = None
    for f in files:
        o = f.Get(f"{d}/{name}")
        if not o:
            continue
        if h is None:
            h = o.Clone(name + "_sum")
            h.SetDirectory(0)
        else:
            h.Add(o)
    return h


def table(name, tab):
    s = f"l{name} = "
    lines = [", ".join(f"{tab.get((p, ip), 0.0):12.6f}" for p in ALLPL) for ip in range(1, 12)]
    return s + ("\n" + " " * len(s)).join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True, help="comma list of lad_edep_histos.C outputs for the run group")
    ap.add_argument("--dir", required=True, help="period directory inside them (p0_lt22375, p1_22375, ...)")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--refdir", default="p5_22698")
    ap.add_argument("--in", dest="inp", required=True, help="constants used for the group so far")
    ap.add_argument("--refparam", required=True, help="constants of the reference group")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--max-err", type=float, default=0.25)
    a = ap.parse_args()

    gfiles = [ROOT.TFile.Open(f) for f in a.group.split(",")]
    rfile = ROOT.TFile.Open(a.ref)
    Li, Ci = read_param_table(a.inp, "ladhodo_LCoeff"), read_param_table(a.inp, "ladhodo_cableFit")
    Lr, Cr = read_param_table(a.refparam, "ladhodo_LCoeff"), read_param_table(a.refparam, "ladhodo_cableFit")
    d = lambda pl, b: (Li[(pl, b)] - Lr[(pl, b)]) - (Ci[(pl, b)] - Cr[(pl, b)])
    res = {}
    for front, back in PAIRS:
        for b in range(1, 12):
            name = f"h_dtb_amp_{back}_{b}"
            x, ex = peak(summed(gfiles, a.dir, name))
            x5, ex5 = peak(rfile.Get(f"{a.refdir}/{name}"))
            if x is None or x5 is None:
                continue
            r = (x - x5) + d(back, b) - d(front, b)
            res[(back, b)] = (r, float(np.hypot(ex, ex5)))
    L = dict(Li)
    notes = []
    for _, back in PAIRS:
        good = [v[0] for (pl, b), v in res.items() if pl == back and v[1] <= a.max_err]
        med = float(np.median(good)) if good else 0.0
        for b in range(1, 12):
            r, e = res.get((back, b), (None, None))
            if r is None or e > a.max_err:
                L[(back, b)] = Li[(back, b)] - med
                notes.append(f"{back}-{b}: no usable pair peak, plane median {med:+.2f} ns")
                print(f"{back}-{b:2d}: plane median {med:+.3f}")
            else:
                L[(back, b)] = Li[(back, b)] - r
                print(f"{back}-{b:2d}: residual {r:+.3f} +- {e:.3f} ns -> LCoeff {Li[(back, b)]:+.3f} -> {L[(back, b)]:+.3f}")
    rr = np.array([v[0] for v in res.values() if v[1] <= a.max_err])
    print(f"bars with a peak: {len(rr)}, mean residual {rr.mean():+.3f}, rms {rr.std():.3f}")

    txt = open(a.inp).read()
    with open(a.out, "w") as f:
        f.write(f"; LAD hodoscope velFit, cableFit, LCoeff for runs {a.label}\n")
        f.write(f"; = {os.path.basename(a.inp)} with the back-plane (001, 101) LCoeff shifted by the same-paddle\n")
        f.write(";   front-back fast-pair residual of this run group against runs 22698+ (lad_pairs_period_offsets.py).\n")
        f.write(";   Fixes t_back - t_front per paddle; front planes, plane 200 and cables as in the input file.\n")
        f.write(f";   {len(rr)} bars measured, mean residual {rr.mean():+.2f} ns, rms {rr.std():.2f} ns\n")
        for n in notes:
            f.write(f";   {n}\n")
        f.write(";\n")
        f.write(table("ladhodo_velFit", read_param_table(a.inp, "ladhodo_velFit")) + "\n")
        f.write(table("ladhodo_cableFit", Ci) + "\n")
        f.write(table("ladhodo_LCoeff", L) + "\n")
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            t = read_param_table(a.inp, nm)
            if t:
                f.write(table(nm, t) + "\n")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
