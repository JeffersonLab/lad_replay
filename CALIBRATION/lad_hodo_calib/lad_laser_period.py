#!/usr/bin/env python3
"""LAD bar timing for a PMT-HV period without a usable photon flash (runs before 22536, where the PMT gains were
too low for the photon hits to pass threshold), from a TDC-only laser run of that period, relative to the reference
(P5) laser calibration.

The pre-22590 laser runs have no FADC pulses, so the walk at the laser amplitude cannot be measured there. The
amplitude is estimated as A_per = A_ref * r, with A_ref the FADC amplitude of the same PMT in reference laser runs at
the same attenuation and r the period/reference MIP gain ratio of the bar (lad_edep_transfer.py; the median of the
period where a bar has none). Per PMT,
    D = [t_per - tw(A_per)] - [t_ref - tw(A_ref)]          (t: TDC - photodiode peak without ADC requirement)
and, with the lad_laser_offsets.py conventions (cableFit = (t_btm - t_top)/2, LCoeff = T_ref - T_bar, T = t_top),
    cableFit_per = cableFit_ref + (D_btm - D_top) / 2
    LCoeff_per   = LCoeff_ref - (D_top - median over bars of D_top)
The median reference keeps the average LCoeff (the global LAD timing) of the reference calibration.

usage: lad_laser_period.py out.param tw.param ref_vp.param transfer_amp.csv --period p0_lt22375
           --per laser_22358.root --ref laser_22824.root laser_23249.root laser_23798.root [--label "runs < 22375"]
"""
import argparse
import csv
import os
import sys

import numpy as np

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_laser_offsets import read_param_table, read_thr, tw  # noqa: E402
from lad_tw_fit import peak_fit  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
PLANES = ["000", "001", "100", "101", "200"]
ALLPL = PLANES + ["REFBAR"]


def tdc_peaks(fn):
    """TDC - PD peak (no ADC requirement) and mean FADC amplitude per PMT."""
    f = ROOT.TFile.Open(fn)
    out = {}
    for pl in PLANES:
        for b in range(1, 12):
            for s in ["Top", "Btm"]:
                h = f.Get(f"{pl}/h_tdc_{pl}_{s}_{b}")
                r = peak_fit(h) if h and h.GetEntries() > 500 else None
                h2 = f.Get(f"{pl}/h_las_{pl}_{s}_{b}")
                amp = None
                if h2 and h2.GetEntries() > 200:
                    px = h2.ProjectionX("_px", 1, h2.GetNbinsY())
                    px.GetXaxis().SetRangeUser(15, px.GetXaxis().GetXmax())
                    amp = px.GetMean() if px.Integral() > 200 else None
                if r:
                    out[(pl, b, s)] = (r[0], amp)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("twparam")
    ap.add_argument("ref_vp")
    ap.add_argument("transfer_csv")
    ap.add_argument("--period", required=True)
    ap.add_argument("--per", required=True, help="laser histograms (lad_laser_histos.C) of the period")
    ap.add_argument("--ref", nargs="+", required=True, help="reference-period laser histograms at the same attenuation")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    thr = read_thr(a.twparam)
    c2 = {s: read_param_table(a.twparam, f"ladhodo_c2_{s}") for s in ["Top", "Btm"]}
    c3 = {s: read_param_table(a.twparam, f"ladhodo_c3_{s}") or {} for s in ["Top", "Btm"]}
    ratio = {}
    for r in csv.DictReader(open(a.transfer_csv)):
        if r["period"] == a.period and r["gain_ratio"] not in ("", "nan"):
            ratio[(r["plane"], int(r["paddle"]))] = float(r["gain_ratio"])
    rmed = float(np.median(list(ratio.values())))
    print(f"{a.period}: gain ratio for {len(ratio)} bars, median {rmed:.3f} (used for the others)")

    per = tdc_peaks(a.per)
    refs = [tdc_peaks(fn) for fn in a.ref]
    D, notes = {}, []
    for pl in PLANES:
        for b in range(1, 12):
            r = ratio.get((pl, b), rmed)
            for s in ["Top", "Btm"]:
                k = (pl, b, s)
                if k not in per:
                    continue
                ds = []
                for ref in refs:
                    if k not in ref or ref[k][1] is None:
                        continue
                    t_ref, A_ref = ref[k]
                    A_per = A_ref * r
                    p2, p3 = c2[s][(pl, b)], c3[s].get((pl, b), 1.0)
                    ds.append((per[k][0] - tw(A_per, p2, p3, thr)) - (t_ref - tw(A_ref, p2, p3, thr)))
                if ds:
                    D[k] = float(np.median(ds))
    tops = [v for k, v in D.items() if k[2] == "Top"]
    medtop = float(np.median(tops))
    vel = read_param_table(a.ref_vp, "ladhodo_velFit")
    cab = read_param_table(a.ref_vp, "ladhodo_cableFit")
    lco = read_param_table(a.ref_vp, "ladhodo_LCoeff")
    newc, newl = dict(cab), dict(lco)
    print(" bar     D_top   D_btm   dcable   dLCoeff  (ns)")
    for pl in PLANES:
        for b in range(1, 12):
            kt, kb = (pl, b, "Top"), (pl, b, "Btm")
            if kt in D and kb in D:
                newc[(pl, b)] = cab[(pl, b)] + 0.5 * (D[kb] - D[kt])
                newl[(pl, b)] = lco[(pl, b)] - (D[kt] - medtop)
                print(f"{pl}-{b:2d}  {D[kt]:+6.2f}  {D[kb]:+6.2f}  {0.5 * (D[kb] - D[kt]):+6.2f}  {-(D[kt] - medtop):+6.2f}")
            else:
                notes.append(f"{pl}-{b}: no laser TDC peak on both ends, reference values")
            if (pl, b) not in ratio:
                notes.append(f"{pl}-{b}: no MIP gain ratio, period median {rmed:.3f} used for the laser amplitude")
    dl = np.array([newl[(pl, b)] - lco[(pl, b)] for pl in PLANES for b in range(1, 12)])
    dc = np.array([newc[(pl, b)] - cab[(pl, b)] for pl in PLANES for b in range(1, 12)])
    print(f"rms change: LCoeff {np.std(dl):.2f} ns, cableFit {np.std(dc):.2f} ns")

    def table(name, tab):
        s = f"l{name} = "
        lines = [", ".join(f"{tab.get((pl, ip), 0.0):12.6f}" for pl in ALLPL) for ip in range(1, 12)]
        return s + ("\n" + " " * len(s)).join(lines) + "\n"

    with open(a.out, "w") as f:
        f.write(f"; LAD hodoscope velFit, cableFit, LCoeff for runs {a.label}\n")
        f.write(f"; reference {os.path.basename(a.ref_vp)} shifted per PMT by the TDC-only laser run "
                f"{os.path.basename(a.per)} vs {', '.join(os.path.basename(x) for x in a.ref)}\n")
        f.write(f"; laser amplitude = reference amplitude x {a.period} MIP gain ratio; time walk "
                f"{os.path.basename(a.twparam)} (CALIBRATION/lad_hodo_calib/lad_laser_period.py)\n")
        f.write("; no photon flash in this period to check these (PMT gains too low); global LAD timing kept\n")
        for n in notes:
            f.write(f";   {n}\n")
        f.write(";" + "".join(f"{p:>14s}" for p in ALLPL) + "\n")
        f.write(table("ladhodo_velFit", vel) + "\n")
        f.write(table("ladhodo_cableFit", newc) + "\n")
        f.write(table("ladhodo_LCoeff", newl) + "\n")
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            t = read_param_table(a.ref_vp, nm)
            if t:
                f.write(table(nm, t) + "\n")
    print(f"wrote {a.out} ({len(notes)} notes)")


if __name__ == "__main__":
    main()
