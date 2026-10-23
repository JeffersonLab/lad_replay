#!/usr/bin/env python3
"""Per-bar photon-flash peaks from lad_bar_timing_check.C output, for the old and new calibrations.

The photon peak is fit (iterative Gaussian core) in t_hit - t_vertex,RF - L/c for every bar, using back-plane
hits with a front veto for planes 001/101 and all hits for planes 000/100/200. Prints the closure test and the
per-bar peak positions relative to the median of all bars, and writes out_prefix.csv.

usage: lad_bar_timing_peaks.py btc.root out_prefix [--window 1700,1760]
"""
import argparse
import csv

import numpy as np

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
PLANES = ["000", "001", "100", "101", "200"]


def photon_peak(h, lo, hi, rebin=3):
    """Photon peak: maximum of the rebinned histogram in [lo, hi], then a least-squares Gaussian + linear
    background fit on [mu - 5, mu + 2] ns (the slow-particle tail is on the late side), iterated twice."""
    from scipy.optimize import curve_fit

    h = h.Clone()
    h.Rebin(rebin)
    ax = h.GetXaxis()
    x = np.array([ax.GetBinCenter(i) for i in range(1, h.GetNbinsX() + 1)])
    y = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])
    win = (x > lo) & (x < hi)
    if y[win].sum() < 200:
        return None
    ys = np.convolve(y, np.ones(3) / 3, mode="same")
    mu = x[win][np.argmax(ys[win])]
    res = None
    for _ in range(2):
        sel = (x > mu - 5) & (x < mu + 2)
        xs, yy = x[sel], y[sel]
        b0 = np.median(y[(x > mu - 8) & (x < mu - 4)]) if np.any((x > mu - 8) & (x < mu - 4)) else yy.min()
        f = lambda t, A, m, sg, c0, c1: A * np.exp(-0.5 * ((t - m) / sg) ** 2) + c0 + c1 * (t - mu)
        try:
            p, cov = curve_fit(f, xs, yy, p0=[max(yy.max() - b0, 1), mu, 0.8, b0, 0.0],
                               sigma=np.sqrt(np.maximum(yy, 1)), bounds=([0, mu - 3, 0.15, -np.inf, -np.inf],
                                                                         [np.inf, mu + 3, 4.0, np.inf, np.inf]))
        except Exception:
            break
        mu = p[1]
        res = (p[1], float(np.sqrt(cov[1, 1])), p[2], float(y[win].sum()), p[0] / max(b0, 1))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("btc")
    ap.add_argument("prefix")
    ap.add_argument("--window", default="1725.5,1730.5", help="photon peak search window (the RF comb repeats every 4 ns)")
    ap.add_argument("--exclude", default="101:10,200:7", help="plane:paddle list left out of the summed check")
    a = ap.parse_args()
    lo, hi = (float(x) for x in a.window.split(","))
    f = ROOT.TFile.Open(a.btc)
    hc = f.Get("h_closure")
    print(f"closure (old recomputed - replay GoodHitTimeAvg): entries {hc.GetEntries():.0f}, mean {hc.GetMean():.4f}, "
          f"rms {hc.GetRMS():.4f}, frac |d|<0.05: {hc.Integral(hc.FindBin(-0.05), hc.FindBin(0.05)) / max(hc.Integral(), 1):.3f}")
    rows = []
    for pl in PLANES:
        cat = "bveto" if pl in ("001", "101") else "all"
        for b in range(1, 12):
            row = dict(plane=pl, paddle=b, cat=cat)
            for s in ["old", "new"]:
                r = photon_peak(f.Get(f"h_{s}_{cat}_{pl}_{b}"), lo, hi)
                row[f"{s}_mu"] = r[0] if r else np.nan
                row[f"{s}_emu"] = r[1] if r else np.nan
                row[f"{s}_sig"] = r[2] if r else np.nan
            rows.append(row)
    for s in ["old", "new"]:
        mus = np.array([r[f"{s}_mu"] for r in rows])
        med = np.nanmedian(mus)
        for r in rows:
            r[f"{s}_rel"] = r[f"{s}_mu"] - med
        print(f"{s}: median peak {med:.2f}; rms of bar peaks {np.nanstd(mus):.3f} ns; per plane:")
        for pl in PLANES:
            v = np.array([r[f"{s}_rel"] for r in rows if r["plane"] == pl])
            sg = np.array([r[f"{s}_sig"] for r in rows if r["plane"] == pl])
            print(f"   {pl}: rel " + " ".join(f"{x:+5.2f}" for x in v) + f"   | median sigma {np.nanmedian(sg):.2f}")
    # summed photon peak: back planes with front veto, and front/200 all hits, each bar shifted to its own peak
    excl = {tuple(x.split(":")) for x in a.exclude.split(",") if x}
    for cat, pls in [("bveto", ["001", "101"]), ("all", ["000", "100", "200"])]:
        for s_ in ["old", "new"]:
            tot = None
            for pl in pls:
                for b in range(1, 12):
                    if (pl, str(b)) in excl:
                        continue
                    h = f.Get(f"h_{s_}_{cat}_{pl}_{b}")
                    if tot is None:
                        tot = h.Clone(f"tot_{cat}_{s_}")
                    else:
                        tot.Add(h)
            r = photon_peak(tot, lo, hi)
            print(f"summed {cat} {'+'.join(pls)} {s_}: peak {r[0]:.2f} sigma {r[2]:.3f} ns (no per-bar realignment)" if r else f"summed {cat} {s_}: fit failed")
    with open(a.prefix + ".csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
