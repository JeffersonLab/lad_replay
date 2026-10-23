#!/usr/bin/env python3
"""Write the parameter table for lad_photon_sigma.C: one (TW, VP) pair of LADlib param files per set, in order.

usage: lad_photon_sigma_params.py out.txt TW0.param VP0.param [TW1.param VP1.param ...]
"""
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


out, files = sys.argv[1], sys.argv[2:]
assert len(files) % 2 == 0 and files, __doc__
with open(out, "w") as f:
    for i in range(len(files) // 2):
        f.write("\n".join(rows(i, files[2 * i], files[2 * i + 1])) + "\n")
