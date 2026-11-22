#!/usr/bin/env python3
"""Cable offset and bar offset of a back-plane bar (001/101) from same-paddle front-back pairs, for bars without
laser data (101-10). The front bar keeps its (laser) constants; the back bar is moved so that
  - the y_back - y_front peak is at 0:       cableFit += dy / velFit
  - the t_back - t_front peak equals the median peak of the other bars of the plane:
                                             LCoeff   += target - dt + dy / velFit
(the hit time is 0.5 (t_top + t_btm) - cableFit + LCoeff, so the cable change moves it by -dy/velFit).

Input: the summed lad_pair_histos.py output made with the parameter file being corrected as set 1 ("new").

usage: lad_bar_from_pairs.py pairs_sum.npz in.param out.param --bars 101:10
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lad_laser_offsets import read_param_table  # noqa: E402

ALLPL = ["000", "001", "100", "101", "200", "REFBAR"]


def gpeak(h, c, half):
    """Gaussian peak position near the maximum of h (bin centres c), fit within +-half."""
    i = int(np.argmax(np.convolve(h, np.ones(5) / 5, "same")))
    m = np.abs(c - c[i]) < half
    p, cov = curve_fit(lambda x, A, mu, s: A * np.exp(-0.5 * ((x - mu) / s) ** 2), c[m], h[m],
                       p0=[h[m].max(), c[i], half / 2], sigma=np.sqrt(np.maximum(h[m], 1)))
    return float(p[1]), float(np.sqrt(cov[1, 1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs")
    ap.add_argument("inparam")
    ap.add_argument("outparam")
    ap.add_argument("--bars", default="101:10", help="comma list of plane:paddle (back planes only)")
    a = ap.parse_args()
    S = np.load(a.pairs)
    yc = 0.5 * (S["ybins"][1:] + S["ybins"][:-1])
    tc = 0.5 * (S["tbins"][1:] + S["tbins"][:-1])
    vel = read_param_table(a.inparam, "ladhodo_velFit")
    cab = read_param_table(a.inparam, "ladhodo_cableFit")
    lco = read_param_table(a.inparam, "ladhodo_LCoeff")
    bars = [(x.split(":")[0], int(x.split(":")[1])) for x in a.bars.split(",")]
    notes = []
    for pl, b in bars:
        dts = {k: gpeak(S[f"new_dt_{pl}_{k}"], tc, 0.8)[0] for k in range(1, 12) if (pl, k) not in bars}
        target = float(np.median(list(dts.values())))
        dy, edy = gpeak(S[f"new_dy_{pl}_{b}"], yc, 15)
        dt, edt = gpeak(S[f"new_dt_{pl}_{b}"], tc, 0.8)
        dc = dy / vel[(pl, b)]
        dl = target - dt + dc
        print(f"{pl}-{b}: dy {dy:+.1f} +- {edy:.1f} cm, dt {dt:+.2f} +- {edt:.2f} ns (target {target:+.2f}): "
              f"cableFit {cab[(pl, b)]:+.3f} -> {cab[(pl, b)] + dc:+.3f}, LCoeff {lco[(pl, b)]:+.3f} -> {lco[(pl, b)] + dl:+.3f}")
        notes.append(f"{pl}-{b}: cableFit and LCoeff from front-back pairs (lad_bar_from_pairs.py): dy {dy:+.1f} cm, "
                     f"dt {dt:+.2f} ns -> {target:+.2f}")
        cab[(pl, b)] += dc
        lco[(pl, b)] += dl

    def table(name, tab):
        s = f"l{name} = "
        lines = [", ".join(f"{tab.get((p, ip), 0.0):12.6f}" for p in ALLPL) for ip in range(1, 12)]
        return s + ("\n" + " " * len(s)).join(lines) + "\n"

    txt = open(a.inparam).read()
    head = [l for l in txt.splitlines() if l.startswith(";")]
    with open(a.outparam, "w") as f:
        for l in head[:-1]:
            f.write(l + "\n")
        for n in notes:
            f.write(f";   {n}\n")
        f.write(head[-1] + "\n")
        f.write(table("ladhodo_velFit", vel) + "\n")
        f.write(table("ladhodo_cableFit", cab) + "\n")
        f.write(table("ladhodo_LCoeff", lco) + "\n")
        for nm in ["ladhodo_velFit_FADC", "ladhodo_cableFit_FADC", "ladhodo_LCoeff_FADC"]:
            t = read_param_table(a.inparam, nm)
            if t:
                f.write(table(nm, t) + "\n")
    print("wrote", a.outparam)


if __name__ == "__main__":
    main()
