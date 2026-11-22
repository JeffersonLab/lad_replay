#!/usr/bin/env python3
"""LAD energy calibration for every PMT-HV period: the reference-period (proton-scale) constants transferred
bar by bar with the MIP peak, which has plenty of statistics in every period.

    c_X(bar) = c_ref(bar) * MIP_ref(bar) / MIP_X(bar)          (c in MeV per ADC unit)

Bars without a MIP peak in a period (plane 200 before 22590, where the LAD ToF of the existing replays has no
fast-particle peak) are scaled by the median gain ratio of that period and flagged.
MIP peaks: planes 000/001/100/101 from fast same-paddle front-back pairs (t_back - t_front ~ 1.3 ns, no ToF
needed); plane 200 from the 1/beta ~ 1 slice of the "all" category (accidental-subtracted). Reference constants: lad_edep_fit.py
(front, 200) and lad_edep_back.py (back; MIP match by default, --back-method dt for the Delta-t fit).

usage: lad_edep_transfer.py histos.root fit_<var>.csv back_<var>.csv out.csv --var int [--ref p5_22698]
"""
import argparse
import csv
import os
import sys

import numpy as np

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_edep_fit import acc_spectrum, mip_peak  # noqa: E402
from lad_edep_back import flat_sub, mip  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
PLANES = ["000", "001", "100", "101", "200"]
BACK = {"001", "101"}


FRONT_OF = {"000": "001", "100": "101"}  # front plane -> back plane whose pair histograms hold it


def mip_of(d, var, pl, b):
    """MIP peak: fast same-paddle front-back pairs (no ToF needed) for planes 000/001/100/101; the 1/beta ~ 1
    slice of all hits for plane 200 (needs a working LAD ToF)."""
    if pl in BACK or pl in FRONT_OF:
        name = f"h_dtb_{var}_{pl}_{b}" if pl in BACK else f"h_dtf_{var}_{FRONT_OF[pl]}_{b}"
        h = d.Get(name)
        if not h or h.GetEntries() < 500:
            return None
        return mip(flat_sub(h))
    h = d.Get(f"h_{var}_all_{pl}_{b}")
    if not h or h.GetEntries() < 500:
        return None
    return mip_peak(h, acc=acc_spectrum(h))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histos")
    ap.add_argument("fitcsv")
    ap.add_argument("backcsv")
    ap.add_argument("out")
    ap.add_argument("--var", default="int", choices=["int", "amp"])
    ap.add_argument("--ref", default="p5_22698")
    ap.add_argument("--back-method", default="mip", choices=["mip", "dt"])
    ap.add_argument("--skip-mip", default="", help="plane:paddle list whose MIP peaks in these histograms are not "
                    "trusted (e.g. pairs formed with wrong constants); they take the period-median gain ratio")
    ap.add_argument("--ref-histos", default="", help="histogram file for the reference-period MIPs (default: histos)")
    a = ap.parse_args()
    skip = {(x.split(":")[0], int(x.split(":")[1])) for x in a.skip_mip.split(",") if x}

    ref = {}
    for r in csv.DictReader(open(a.fitcsv)):
        if r["period"] == a.ref and r["plane"] not in BACK and r.get("g") not in (None, "", "nan"):
            ref[(r["plane"], int(r["paddle"]))] = (1.0 / float(r["g"]), "tof")
    key, alt = ("g_mip", "g_dt") if a.back_method == "mip" else ("g_dt", "g_mip")
    for r in csv.DictReader(open(a.backcsv)):
        if r["period"] != a.ref:
            continue
        if r.get(key) not in (None, "", "nan"):
            ref[(r["plane"], int(r["paddle"]))] = (1.0 / float(r[key]), a.back_method)
        elif r.get(alt) not in (None, "", "nan"):  # e.g. MIP match impossible when the front bar has no fit
            ref[(r["plane"], int(r["paddle"]))] = (1.0 / float(r[alt]), alt[2:])

    fin = ROOT.TFile.Open(a.histos)
    periods = [k.GetName() for k in fin.GetListOfKeys() if k.ReadObj().InheritsFrom("TDirectory")]
    fref = ROOT.TFile.Open(a.ref_histos) if a.ref_histos else fin
    dref = fref.Get(a.ref)
    mref = {(pl, b): (None if (pl, b) in skip else mip_of(dref, a.var, pl, b)) for pl in PLANES for b in range(1, 12)}
    rows = []
    for per in periods:
        d = fin.Get(per)
        for pl in PLANES:
            for b in range(1, 12):
                k = (pl, b)
                m = (None if k in skip else mip_of(d, a.var, pl, b)) if per != a.ref else mref[k]
                row = dict(period=per, plane=pl, paddle=b, mip=m if m else np.nan, mip_ref=mref[k] or np.nan)
                if k in ref:
                    row["ref_c"], row["ref_method"] = ref[k]
                    if per == a.ref:
                        row["c"] = ref[k][0]
                    elif m and mref[k]:
                        row["c"] = ref[k][0] * mref[k] / m
                    row["gain_ratio"] = (m / mref[k]) if (m and mref[k]) else np.nan
                rows.append(row)
    # bars with a reference constant but no MIP peak in a period: scale by the median gain ratio of that period
    for per in periods:
        gr = [r["gain_ratio"] for r in rows if r["period"] == per and np.isfinite(r.get("gain_ratio", np.nan))]
        if not gr:
            continue
        med = float(np.median(gr))
        for r in rows:
            if r["period"] == per and "ref_c" in r and not r.get("c"):
                r["c"] = r["ref_c"] / med
                r["ref_method"] = r["ref_method"] + "+period_median_ratio"
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["period", "plane", "paddle", "mip", "mip_ref", "gain_ratio", "ref_c",
                                          "ref_method", "c"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    for per in periods:
        gr = [r["gain_ratio"] for r in rows if r["period"] == per and np.isfinite(r.get("gain_ratio", np.nan))]
        nc = sum(1 for r in rows if r["period"] == per and r.get("c"))
        if gr:
            print(f"{per}: {nc} bars calibrated; MIP gain ratio to {a.ref}: median {np.median(gr):.3f}, "
                  f"range {np.min(gr):.2f}-{np.max(gr):.2f}")


if __name__ == "__main__":
    main()
