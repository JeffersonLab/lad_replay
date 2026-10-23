#!/usr/bin/env python3
"""Photon-peak widths from lad_photon_sigma.C outputs.

usage: lad_photon_sigma_fit.py [--csv out.csv] LABEL=file1.root[,file2.root...] ...
Each LABEL sums its files (e.g. one per HV period); its first letter (P/H) is the reference arm.
"""
import sys
import numpy as np
import ROOT
from scipy.optimize import curve_fit

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
NSET, SELS, RFS = 4, ["all", "back", "bveto", "bvetoE"], ["norf", "rf"]


def get(files, key):
    h = None
    for f in files:
        tf = ROOT.TFile.Open(f)
        x = tf.Get(key)
        if h is None:
            h = x.Clone(); h.SetDirectory(0)
        else:
            h.Add(x)
        tf.Close()
    return h


CENTRE = {"P": 1728.0, "H": 1752.0}  # rough photon peak (ns, before the global offset) to search around


def fit(h, c0, search=3.0, w=1.6):
    """Gaussian + linear background in +-w ns around the highest 5-bin-smoothed bin within c0 +- search (as
    analysis2026/analyze.py gfit, with a wider fixed window). Returns (mean, sigma, sigma error, peak counts,
    sigma from the FWHM above the fitted linear background)."""
    h = h.Clone(); h.Rebin(2)  # 0.1 ns
    c = np.array([h.GetBinCenter(b) for b in range(1, h.GetNbinsX() + 1)])
    y = np.array([h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)])
    s = np.abs(c - c0) < search
    i = int(np.argmax(np.convolve(y, np.ones(5) / 5, "same") * s)); mu = c[i]; m = np.abs(c - mu) < w
    f = lambda t, A, m0, sg, b0, b1: A * np.exp(-0.5 * ((t - m0) / sg) ** 2) + b0 + b1 * (t - mu)
    try:
        p, cv = curve_fit(f, c[m], y[m], p0=[max(y[i] - np.median(y[m]), 1), mu, w / 3, np.median(y[m]), 0],
                          sigma=np.sqrt(np.maximum(y[m], 1)), maxfev=5000)
    except Exception:
        return np.nan, np.nan, np.nan, 0, np.nan
    if not (abs(p[2]) < w and abs(p[1] - mu) < w):
        return np.nan, np.nan, np.nan, 0, np.nan
    return p[1], abs(p[2]), np.sqrt(cv[2, 2]), p[0] * abs(p[2]) * 2.5066 / (c[1] - c[0]), fwhm_sigma(c, y, p[1], mu, p[3], p[4])


def fwhm_sigma(c, y, m0, mu, b0, b1, half=2.5):
    """FWHM / 2.355 of the 3-bin-smoothed peak above the linear background b0 + b1 (t - mu), within m0 +- half."""
    s = np.abs(c - m0) < half
    x, r = c[s], np.convolve(y, np.ones(3) / 3, "same")[s] - (b0 + b1 * (c[s] - mu))
    i = int(np.argmax(r)); hm = r[i] / 2
    lo = next((k for k in range(i, 0, -1) if r[k - 1] < hm), None)
    hi = next((k for k in range(i, len(r) - 1) if r[k + 1] < hm), None)
    if lo is None or hi is None:
        return np.nan
    xl = x[lo - 1] + (hm - r[lo - 1]) / (r[lo] - r[lo - 1]) * (x[lo] - x[lo - 1])
    xh = x[hi] + (r[hi] - hm) / (r[hi] - r[hi + 1]) * (x[hi + 1] - x[hi])
    return (xh - xl) / 2.3548


args = sys.argv[1:]
csv = None
if args[:1] == ["--csv"]:
    csv, args = args[1], args[2:]
rows = []
for a in args:
    lab, fl = a.split("=", 1)
    files = fl.split(",")
    for s in range(NSET):
        for r in RFS:
            for c in SELS:
                mu, sg, se, n, sf = fit(get(files, f"h_s{s}_{r}_{c}"), CENTRE[lab[0]])
                rows.append((lab, s, r, c, mu, sg, se, n, sf))
                print(f"{lab:8s} set {s} {r:4s} {c:7s} mean {mu:8.2f}  sigma {sg:.3f} +- {se:.3f}  FWHM/2.355 {sf:.3f}"
                      f"  n {n:9.0f}")
if csv:
    with open(csv, "w") as f:
        f.write("label,set,rf,sel,mean,sigma,sigma_err,n,sigma_fwhm\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
