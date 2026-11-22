#!/usr/bin/env python3
"""Write the LAD edep calibration param files (one per PMT-HV period) from lad_edep_fit.py output.

usage: lad_edep_param.py transfer_int.csv transfer_amp.csv outdir   (lad_edep_transfer.py or lad_edep_fit.py csv files) [--fallback p2_22533:p1_22375,p4_22695:p5_22698]

Each file holds lladhodo_adcInt2MeV (MeV per pC of sqrt(int_top int_btm)) and lladhodo_adcAmp2MeV (MeV per mV of
sqrt(amp_top amp_btm)) in the 11 x 6 layout of the other LAD hodo params (REFBAR = 1). Bars without a fit take the
value from the --fallback period, else the median of their plane in that period, and are listed in the header.
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np

PLANES = ["000", "001", "100", "101", "200", "REFBAR"]
PERIOD_RUNS = {"p0_lt22375": "runs < 22375", "p1_22375": "22375-22532", "p2_22533": "22533-22535",
               "p3_22536": "22536-22694", "p4_22695": "22695-22697", "p5_22698": "22698 and later"}


def load(fn):
    out = defaultdict(dict)
    for r in csv.DictReader(open(fn)):
        v = r.get("MeV_per_unit_final") or r.get("c") or ""
        if v not in ("", "nan"):
            out[r["period"]][(r["plane"], int(r["paddle"]))] = (float(v), r.get("method", r.get("ref_method", "")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fit_int")
    ap.add_argument("fit_amp")
    ap.add_argument("outdir")
    ap.add_argument("--fallback", default="p2_22533:p1_22375,p4_22695:p5_22698")
    ap.add_argument("--max-dev", type=float, default=0.6, help="reject fits more than this fraction from the plane median")
    a = ap.parse_args()
    fb = dict(x.split(":") for x in a.fallback.split(",") if x)
    data = {"int": load(a.fit_int), "amp": load(a.fit_amp)}
    os.makedirs(a.outdir, exist_ok=True)
    for per in PERIOD_RUNS:
        tables, notes = {}, []
        for var in ["int", "amp"]:
            d = data[var].get(per, {})
            tab = {}
            for pl in PLANES[:5]:
                vals = [d[(pl, b)][0] for b in range(1, 12) if (pl, b) in d]
                med = float(np.median(vals)) if vals else np.nan
                for b in range(1, 12):
                    v = d.get((pl, b), (None,))[0]
                    if v is not None and np.isfinite(med) and abs(v / med - 1) > a.max_dev:
                        notes.append(f"{var} {pl}-{b}: fit {v:.3f} rejected (plane median {med:.3f})")
                        v = None
                    if v is None and per in fb and (pl, b) in data[var].get(fb[per], {}):
                        v = data[var][fb[per]][(pl, b)][0]
                        notes.append(f"{var} {pl}-{b}: from {fb[per]}")
                    if v is None and np.isfinite(med):
                        v = med
                        notes.append(f"{var} {pl}-{b}: plane median")
                    tab[(pl, b)] = v if v is not None else 1.0
            tables[var] = tab
        fn = os.path.join(a.outdir, f"ladhodo_edep_{per}.param")
        with open(fn, "w") as f:
            f.write(f"; LAD hodoscope energy calibration, PMT-HV period {per} ({PERIOD_RUNS[per]})\n")
            f.write("; proton scale: a proton that just punches through a bar reads its deposited energy\n")
            f.write("; written by CALIBRATION/lad_edep/lad_edep_param.py from lad_edep_fit.py\n")
            for n in notes:
                f.write(f";   {n}\n")
            f.write(";" + "".join(f"{p:>14s}" for p in PLANES) + "\n")
            for var, name in [("int", "lladhodo_adcInt2MeV"), ("amp", "lladhodo_adcAmp2MeV")]:
                pre = f"{name} = "
                lines = []
                for b in range(1, 12):
                    vals = [tables[var].get((pl, b), 1.0) if pl != "REFBAR" else 1.0 for pl in PLANES]
                    lines.append(", ".join(f"{v:12.6f}" for v in vals))
                f.write(pre + ("\n" + " " * len(pre)).join(lines) + "\n\n")
        print(f"wrote {fn} ({len(notes)} notes)")


if __name__ == "__main__":
    main()
