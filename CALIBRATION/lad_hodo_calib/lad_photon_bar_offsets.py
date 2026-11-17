#!/usr/bin/env python3
"""Per-bar LCoeff corrections from the RF-corrected photon flash of a full replay.

Input is the bar table of the 2026 replay analysis (analysis2026/results/bar_timing.txt: period arm plane paddle peak
err sigma n), with the photon peak of tof_rfcorr - path/c per bar and reference arm. For each bar the SHMS and HMS
peaks are combined (weighted mean; they share the LAD calibration, so they should agree), and a bar whose combined
residual r is measured to better than --max-err and differs from 0 by more than --min-sig sigma gets its LCoeff
shifted by -r (the same convention as lad_pairs_period_offsets.py: a bar late by r gets LCoeff - r). velFit, cableFit
and the FADC tables are copied unchanged.

usage: lad_photon_bar_offsets.py --bars bar_timing.txt --periods P5a P5b --inp ladhodo_Vpcalib_2026_p5.param \
         --out new.param --label "22695+"
"""
import argparse
import os
import re

import numpy as np

PL = ["000", "001", "100", "101", "200", "REFBAR"]


def read_param_table(path, name):
    """{(plane, paddle 1-based): value} of a 6-column (plane) x 11-row (paddle) table."""
    txt = open(path).read()
    m = re.search(r"^l" + name + r"\s*=(.*?)(?=^\s*l\w+\s*=|\Z)", txt, re.M | re.S)
    vals = [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+", re.sub(r";.*", "", m.group(1)))]
    return {(PL[i % 6], i // 6 + 1): v for i, v in enumerate(vals)}


def table(name, d):
    rows = []
    for k in range(1, 12):
        rows.append(", ".join(f"{d.get((p, k), 0.0):12.6f}" for p in PL))
    pad = " " * (len(name) + 4)
    return f"l{name} = " + ("\n" + pad).join(rows) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bars", required=True)
    ap.add_argument("--periods", nargs="+", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--max-err", type=float, default=0.08, help="ns; combined error needed to correct a bar")
    ap.add_argument("--min-sig", type=float, default=3.0, help="correct only bars with |r| > min-sig * err")
    ap.add_argument("--skip", nargs="*", default=[], help="bars to leave unchanged, e.g. 100-10 (peak shape unclear)")
    a = ap.parse_args()

    meas = {}
    for line in open(a.bars):
        f = line.split()
        if len(f) < 8 or f[0] not in a.periods: continue
        try: pk, er = float(f[4]), float(f[5])
        except ValueError: continue
        if not (np.isfinite(pk) and np.isfinite(er) and 0 < er < 0.5): continue
        pl = f[2]
        meas.setdefault((pl, int(f[3])), []).append((pk, er, f[0], f[1]))

    L = read_param_table(a.inp, "ladhodo_LCoeff")
    Lnew = dict(L)
    log = []
    for (pl, k), ms in sorted(meas.items()):
        w = np.array([1 / e ** 2 for _, e, _, _ in ms]); x = np.array([p for p, _, _, _ in ms])
        r = float(np.sum(w * x) / w.sum()); e = float(1 / np.sqrt(w.sum()))
        chi2 = float(np.sum(w * (x - r) ** 2)); nd = len(ms) - 1
        apply = e < a.max_err and abs(r) > a.min_sig * e and f"{pl}-{k}" not in a.skip
        if apply: Lnew[(pl, k)] = L[(pl, k)] - r
        log.append(f"{pl}-{k:2d}: photon residual {r:+.3f} +- {e:.3f} ns (chi2 {chi2:.1f}/{nd}, from "
                   f"{' '.join(f'{pp}/{s}:{p:+.2f}' for p, _, pp, s in ms)})"
                   + (f" -> LCoeff {L[(pl, k)]:+.3f} -> {Lnew[(pl, k)]:+.3f}" if apply else " -> unchanged"))
    print("\n".join(log))

    txt = open(a.inp).read()
    with open(a.out, "w") as f:
        f.write(f"; LAD hodoscope velFit, cableFit, LCoeff for runs {a.label}\n")
        f.write(f"; = {os.path.basename(a.inp)} with LCoeff corrected per bar by the RF-corrected photon flash of the 2026\n")
        f.write(f"; full replay (CALIBRATION/lad_hodo_calib/lad_photon_bar_offsets.py; periods {' '.join(a.periods)};\n")
        f.write(f"; bars with error < {a.max_err} ns and |r| > {a.min_sig} sigma; LCoeff -> LCoeff - r"
                + (f"; left unchanged by hand: {' '.join(a.skip)}" if a.skip else "") + ")\n")
        for l in log: f.write(";   " + l + "\n")
        for nm in ["ladhodo_velFit", "ladhodo_cableFit"]:
            f.write(table(nm, read_param_table(a.inp, nm)) + "\n")
        f.write(table("ladhodo_LCoeff", Lnew) + "\n")
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            if re.search(r"^l" + nm + r"\s*=", txt, re.M): f.write(table(nm, read_param_table(a.inp, nm)) + "\n")


if __name__ == "__main__":
    main()
