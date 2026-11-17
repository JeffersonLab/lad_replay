#!/usr/bin/env python3
"""Write the parameter table for lad_bar_timing_check.C from LADlib param files.

usage: lad_bar_timing_params.py out.txt --old TW.param VP.param --new TW.param VP.param
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_laser_offsets import read_param_table, read_thr  # noqa: E402

PLANES = ["000", "001", "100", "101", "200"]


def rows(setn, twp, vpp):
    out = [f"{setn} thr {read_thr(twp)}"]
    c2 = {s: read_param_table(twp, f"ladhodo_c2_{s}") for s in ["Top", "Btm"]}
    c3 = {s: (read_param_table(twp, f"ladhodo_c3_{s}") or {}) for s in ["Top", "Btm"]}
    cab = read_param_table(vpp, "ladhodo_cableFit")
    lco = read_param_table(vpp, "ladhodo_LCoeff")
    vel = read_param_table(vpp, "ladhodo_velFit")
    for pl in PLANES:
        for ip in range(1, 12):
            for s in ["Top", "Btm"]:
                out.append(f"{setn} {pl} {s} {ip} {c2[s][(pl, ip)]:.6f} {c3[s].get((pl, ip), 1.0):.6f} "
                           f"{cab[(pl, ip)]:.6f} {lco[(pl, ip)]:.6f} {vel[(pl, ip)]:.4f}")
    return out


ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--old", nargs=2, required=True)
ap.add_argument("--new", nargs=2, required=True)
a = ap.parse_args()
with open(a.out, "w") as f:
    f.write("\n".join(rows(0, *a.old) + rows(1, *a.new)) + "\n")
