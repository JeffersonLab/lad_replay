#!/usr/bin/env python3
"""Per-period LAD bar offsets: the reference (laser, P5) calibration corrected with the photon flash and the
hit-position distribution of another PMT-HV period.

Both inputs are lad_bar_timing_check.C outputs made with the same "new" parameter set (the reference one):
  - LCoeff: the photon peak of each bar, relative to the median of all bars, is moved to where it is in the
    reference period: corr = (pk_ref - med_ref) - (pk_per - med_per)
    --mode bar: that correction for every bar; --mode plane: the median correction of the bar's plane for every bar
    of the plane; --mode hybrid (default): the plane median, plus the bar's own deviation from it where that is
    larger than --nsig times its fit error (low-statistics periods: keeps the plane-wide shifts, drops the noise)
  - cableFit: the centre of the y distribution of each bar (midpoint of its 50% edges) is moved to where it is
    in the reference period: c_per = c_ref + (y_per - y_ref) / v; LCoeff then also changes by +dc so that the
    hit time is unchanged by the cable correction.
Bars without a usable photon peak or y edges keep the reference values (listed).

usage: lad_period_offsets.py ref_btc.root per_btc.root ref.param out.param --label "22590-22694" [--mode hybrid] [--no-cable]
"""
import argparse
import os
import sys

import numpy as np

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_bar_timing_peaks import photon_peak  # noqa: E402
from lad_laser_offsets import read_param_table  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
PLANES = ["000", "001", "100", "101", "200"]
ALLPL = PLANES + ["REFBAR"]


def peaks(fn, window):
    """Photon peak of each bar relative to the median of all bars, and its fit error."""
    f = ROOT.TFile.Open(fn)
    out, err = {}, {}
    for pl in PLANES:
        cat = "bveto" if pl in ("001", "101") else "all"
        for b in range(1, 12):
            h = f.Get(f"h_new_{cat}_{pl}_{b}")
            r = photon_peak(h, *window)
            if r and r[1] < 0.3:
                out[(pl, b)], err[(pl, b)] = r[0], r[1]
    med = np.median(list(out.values()))
    return {k: v - med for k, v in out.items()}, err, f


def y_centre(f, pl, b, nmin=2000):
    h = f.Get(f"h_y_new_{pl}")
    py = h.ProjectionY("_y", b, b)
    if py.Integral() < nmin:
        return None
    y = np.array([py.GetBinContent(i) for i in range(1, py.GetNbinsX() + 1)])
    x = np.array([py.GetBinCenter(i) for i in range(1, py.GetNbinsX() + 1)])
    ys = np.convolve(y, np.ones(3) / 3, mode="same")
    plateau = np.median(ys[np.abs(x) < 120])
    if plateau <= 0:
        return None
    above = np.where(ys > 0.5 * plateau)[0]
    if len(above) < 10:
        return None

    def cross(i0, i1):  # linear interpolation of the 50% crossing between bins i0 (below) and i1 (above)
        y0, y1 = ys[i0], ys[i1]
        return x[i0] + (0.5 * plateau - y0) / (y1 - y0) * (x[i1] - x[i0]) if y1 != y0 else x[i1]

    lo, hi = above[0], above[-1]
    if lo == 0 or hi == len(x) - 1:
        return None
    return 0.5 * (cross(lo - 1, lo) + cross(hi + 1, hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_btc")
    ap.add_argument("per_btc")
    ap.add_argument("ref_param")
    ap.add_argument("out_param")
    ap.add_argument("--label", default="")
    ap.add_argument("--window", default="1725.5,1730.5")
    ap.add_argument("--no-cable", action="store_true", help="keep the reference cable offsets")
    ap.add_argument("--mode", default="hybrid", choices=["bar", "plane", "hybrid"])
    ap.add_argument("--nsig", type=float, default=3.0, help="hybrid: significance for a per-bar deviation")
    a = ap.parse_args()
    win = tuple(float(x) for x in a.window.split(","))
    pk_ref, er_ref, fref = peaks(a.ref_btc, win)
    pk_per, er_per, fper = peaks(a.per_btc, win)
    # per-bar photon corrections and their errors; plane medians
    dlb = {k: pk_ref[k] - pk_per[k] for k in pk_ref if k in pk_per}
    edl = {k: float(np.hypot(er_ref[k], er_per[k])) for k in dlb}
    dlp = {pl: float(np.median([v for k, v in dlb.items() if k[0] == pl])) for pl in PLANES
           if any(k[0] == pl for k in dlb)}
    print("plane median corrections (ns): " + "  ".join(f"{pl} {v:+.2f}" for pl, v in dlp.items()))
    vel = read_param_table(a.ref_param, "ladhodo_velFit")
    cab = read_param_table(a.ref_param, "ladhodo_cableFit")
    lco = read_param_table(a.ref_param, "ladhodo_LCoeff")
    notes = []
    newc, newl = dict(cab), dict(lco)
    print(" bar    dL(photon)   dy(cm)   dc(ns)")
    for pl in PLANES:
        for b in range(1, 12):
            k = (pl, b)
            if a.mode == "bar":
                dl = dlb.get(k)
            elif a.mode == "plane":
                dl = dlp.get(pl)
            else:
                dl = dlp.get(pl)
                if k in dlb and dl is not None and abs(dlb[k] - dl) > a.nsig * edl[k]:
                    dl = dlb[k]
                    notes.append(f"{pl}-{b}: own photon correction {dl:+.2f} +- {edl[k]:.2f} (plane {dlp[pl]:+.2f})")
            dc = None
            if not a.no_cable:
                yr, yp = y_centre(fref, pl, b), y_centre(fper, pl, b)
                if yr is not None and yp is not None:
                    dc = (yp - yr) / vel[k]
            if dc is not None:
                newc[k] = cab[k] + dc
            elif not a.no_cable:
                notes.append(f"{pl}-{b}: cable from reference")
            if dl is not None:
                newl[k] = lco[k] + dl + (dc or 0.0)
            else:
                newl[k] = lco[k] + (dc or 0.0)
                notes.append(f"{pl}-{b}: no photon correction, LCoeff from reference")
            print(f"{pl}-{b:2d}  {dl if dl is not None else float('nan'):+7.2f}  "
                  f"{(dc * vel[k]) if dc is not None else float('nan'):+7.1f}  {dc if dc is not None else float('nan'):+6.2f}")

    def table(name, tab):
        s = f"l{name} = "
        lines = []
        for ip in range(1, 12):
            lines.append(", ".join(f"{tab.get((pl, ip), 0.0):12.6f}" for pl in ALLPL))
        return s + ("\n" + " " * len(s)).join(lines) + "\n"

    with open(a.out_param, "w") as f:
        f.write(f"; LAD hodoscope velFit, cableFit, LCoeff for runs {a.label}\n")
        f.write(f"; reference {os.path.basename(a.ref_param)} corrected with the photon flash ({a.mode} mode)"
                + (" and hit-position edges\n" if not a.no_cable else "; cables from the reference\n"))
        f.write(f"; (CALIBRATION/lad_hodo_calib/lad_period_offsets.py; inputs {os.path.basename(a.ref_btc)}, "
                f"{os.path.basename(a.per_btc)})\n")
        for n in notes:
            f.write(f";   {n}\n")
        f.write(";" + "".join(f"{p:>14s}" for p in ALLPL) + "\n")
        f.write(table("ladhodo_velFit", vel) + "\n")
        f.write(table("ladhodo_cableFit", newc) + "\n")
        f.write(table("ladhodo_LCoeff", newl) + "\n")
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            t = read_param_table(a.ref_param, nm)
            if t:
                f.write(table(nm, t) + "\n")
    print(f"wrote {a.out_param} ({len(notes)} notes)")


if __name__ == "__main__":
    main()
