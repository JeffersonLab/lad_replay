#!/usr/bin/env python3
"""Back-plane (001, 101) energy calibration from same-paddle front-back pairs, two ways:

1. MIP match: fast particles (t_back - t_front ~ 1.3 ns) deposit the same energy in both 5.08 cm bars, so
   g_back = g_front * MIP_back / MIP_front, with g_front from lad_edep_fit.py (proton scale).
2. Delta-t punch-through: the back-bar light of protons vs the front-to-back flight time, fit with the
   range-energy model (lad_edep_model.py): ADC_back = g * Lnorm(T_back(dt + d)); free g and d. Protons that
   stop in the back bar spread over several ns in dt, so both branches and the cusp are constrained.

usage: lad_edep_back.py histos.root fit_<var>.csv out_prefix [--var int] [--kB 0.0126*rho] [--png]
"""
import argparse
import csv
import os
import sys

import numpy as np
from scipy.optimize import least_squares

import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lad_edep_model as M  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
C_CM_NS = 29.9792458
GAP_CM = 39.7 - M.THICK_CM  # air between the front-bar exit and the back-bar entrance (radial plane spacing 39.7 cm)
PAIRS = [("001", "000"), ("101", "100")]


def beta(T):
    g = 1 + np.asarray(T, dtype=float) / M.M_P
    return np.sqrt(1 - 1 / g**2)


class DtModel:
    """Back-bar light (proton scale) vs t_back - t_front for protons entering the front bar with T_f."""

    def __init__(self, kB, nstep=60):
        Tf = np.arange(82.0, 600.0, 0.5)
        dt = np.zeros_like(Tf)
        Tb = np.zeros_like(Tf)
        rho_air = 1.205e-3 * 0.92  # PVT-equivalent g/cm3
        for i, T in enumerate(Tf):
            # through the front bar in steps, then the air gap
            t, x, E = 0.0, 0.0, T
            dx = M.THICK_CM / nstep
            for _ in range(nstep):
                if E <= 0.05:
                    break
                t += dx / (beta(E) * C_CM_NS)
                E = float(M.t_out(E, dx * M.RHO_PVT))
            if E <= 0.05:
                dt[i], Tb[i] = np.nan, 0.0
                continue
            dxa = GAP_CM / 10
            for _ in range(10):
                t += dxa / (beta(E) * C_CM_NS)
                E = float(M.t_out(E, dxa * rho_air))
            dt[i], Tb[i] = t, E
        ok = np.isfinite(dt) & (Tb > 0.5)
        t_bar = M.T_BAR
        Ept = M.punch_through_T(t_bar)
        Lpt = float(M.light(np.array([Ept]), t_bar, kB)[0]) if kB > 0 else Ept
        L = (M.light(Tb[ok], t_bar, kB) if kB > 0 else M.edep(Tb[ok], t_bar)) * Ept / Lpt
        self.dt, self.L = dt[ok], L  # dt decreasing with T
        i = int(np.argmax(self.L))
        self.dt_cusp, self.Lmax = float(self.dt[i]), float(self.L[i])

    def __call__(self, x, g, d):
        x = np.asarray(x, dtype=float) + d
        return g * np.interp(x, self.dt[::-1], self.L[::-1], left=np.nan, right=0.0)


def flat_sub(h2, lo=8.0, hi=9.9):
    """Subtract the flat (accidental) part, estimated from a late dt sideband, per ADC bin.
    The sideband must lie inside LADlib's front-back matching window (ladhodo_matching_time_tol, 10 ns): pairs
    beyond it do not exist, so a sideband reaching past 10 ns underestimates the accidentals."""
    ax = h2.GetXaxis()
    b1, b2 = ax.FindBin(lo), ax.FindBin(hi)
    acc = h2.ProjectionY("_accdt", b1, b2)
    acc.Scale(1.0 / (b2 - b1 + 1))
    s = h2.Clone(h2.GetName() + "_s")
    for i in range(1, h2.GetNbinsX() + 1):
        for j in range(1, h2.GetNbinsY() + 1):
            s.SetBinContent(i, j, h2.GetBinContent(i, j) - acc.GetBinContent(j))
    return s


def mip(h2s, lo=0.8, hi=2.0):
    """Landau MPV of the ADC distribution of fast pairs."""
    ax = h2s.GetXaxis()
    py = h2s.ProjectionY("_mipdt", ax.FindBin(lo), ax.FindBin(hi))
    if py.Integral() < 300:
        return None
    q = np.zeros(1)
    py.GetQuantiles(1, q, np.array([0.25]))
    py.GetXaxis().SetRangeUser(q[0], py.GetXaxis().GetXmax())
    mu = py.GetBinCenter(py.GetMaximumBin())
    f = ROOT.TF1("lan", "landau", 0.75 * mu, 1.8 * mu)
    f.SetParameters(py.GetMaximum() * 3, mu, 0.1 * mu)
    r = py.Fit(f, "QNRS")
    return f.GetParameter(1) if int(r) == 0 and f.GetParameter(1) > 0 else mu


def ridge(h2s, model, g, d, x_lo, x_hi, w=0.2, nmin=60):
    ax = h2s.GetXaxis()
    pts = []
    x = x_lo
    while x < x_hi:
        c = float(model(np.array([x + w / 2]), g, d)[0])
        if np.isfinite(c) and c > 0:
            py = h2s.ProjectionY("_r", ax.FindBin(x + 1e-6), ax.FindBin(x + w - 1e-6))
            lo, hi = 0.6 * c, 1.5 * c
            k1, k2 = py.GetXaxis().FindBin(lo), py.GetXaxis().FindBin(hi)
            if py.Integral(k1, k2) >= nmin:
                py.GetXaxis().SetRange(k1, k2)
                mu = py.GetBinCenter(py.GetMaximumBin())
                s = 0.12 * mu
                f = ROOT.TF1("gr", "gaus", mu - s, mu + s)
                f.SetParameters(py.GetMaximum(), mu, s)
                r = py.Fit(f, "QNRS")
                if int(r) == 0 and lo < f.GetParameter(1) < hi:
                    pts.append((x + w / 2, f.GetParameter(1), max(f.GetParError(1), 1e-3 * mu)))
        x += w
    return pts


def fit_dt(h2s, model, g0, cusp_excl=0.25):
    g, d = g0, 0.0
    pts, res = [], None
    for _ in range(4):
        pts = ridge(h2s, model, g, d, 1.8, 7.5)
        if len(pts) < 6:
            return None, pts
        P = np.array(pts)
        use = np.abs(P[:, 0] + d - model.dt_cusp) > cusp_excl
        if use.sum() < 5:
            return None, pts
        x, y = P[use, 0], P[use, 1]
        e = np.sqrt(P[use, 2] ** 2 + (0.02 * y) ** 2)
        fun = lambda p: np.nan_to_num((model(x, p[0], p[1]) - y) / e, nan=50.0)
        r = least_squares(fun, [g, d], bounds=([0.05 * g0, -1.5], [20 * g0, 1.5]))
        g, d = r.x
        res = dict(g=g, d=d, chi2=float(np.sum(r.fun**2)), ndf=int(use.sum()) - 2, npts=len(pts))
    return res, pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histos")
    ap.add_argument("fitcsv", help="lad_edep_fit.py csv for the same variable (front gains)")
    ap.add_argument("prefix")
    ap.add_argument("--var", default="int", choices=["int", "amp"])
    ap.add_argument("--kB", type=float, default=0.0126 * M.RHO_PVT)
    ap.add_argument("--periods", default="")
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()

    front = {}
    for r in csv.DictReader(open(a.fitcsv)):
        if r.get("g") not in (None, "", "nan"):
            front[(r["period"], r["plane"], int(r["paddle"]))] = float(r["g"])
    model = DtModel(a.kB)
    print(f"dt model: back cusp at dt = {model.dt_cusp:.2f} ns, Lmax {model.Lmax:.1f} MeV")
    fin = ROOT.TFile.Open(a.histos)
    periods = [k.GetName() for k in fin.GetListOfKeys() if k.ReadObj().InheritsFrom("TDirectory")]
    if a.periods:
        periods = [p for p in periods if p in a.periods.split(",")]
    rows, plots = [], {}
    for per in periods:
        d = fin.Get(per)
        for pb, pf in PAIRS:
            for b in range(1, 12):
                hb = d.Get(f"h_dtb_{a.var}_{pb}_{b}")
                hf = d.Get(f"h_dtf_{a.var}_{pb}_{b}")
                row = dict(period=per, plane=pb, paddle=b, entries=int(hb.GetEntries()) if hb else 0)
                if not hb or hb.GetEntries() < 2000:
                    rows.append(row)
                    continue
                hbs, hfs = flat_sub(hb), flat_sub(hf)
                mb, mf = mip(hbs), mip(hfs)
                gf = front.get((per, pf, b))
                row.update(mip_back=mb or np.nan, mip_front=mf or np.nan, g_front=gf or np.nan)
                if mb and mf and gf:
                    row["g_mip"] = gf * mb / mf
                g0 = row.get("g_mip") or (gf if gf else 0.5)
                r, pts = fit_dt(hbs, model, g0)
                if r:
                    row.update(g_dt=r["g"], d_dt=r["d"], chi2ndf_dt=r["chi2"] / max(r["ndf"], 1), npts_dt=r["npts"])
                plots[(per, pb, b)] = (hbs, pts, r)
                rows.append(row)
                gm, gd = row.get("g_mip", np.nan), row.get("g_dt", np.nan)
                print(f"{per} {pb} {b:2d}: MIP-match 1/g={1/gm if gm==gm else np.nan:.3f}  dt-fit 1/g={1/gd if gd==gd else np.nan:.3f}"
                      f"  ratio {gd/gm if gm==gm and gd==gd else np.nan:.3f}")
    keys = ["period", "plane", "paddle", "entries", "mip_back", "mip_front", "g_front", "g_mip", "g_dt", "d_dt",
            "chi2ndf_dt", "npts_dt"]
    with open(f"{a.prefix}_{a.var}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    rat = [r["g_dt"] / r["g_mip"] for r in rows if r.get("g_dt") and r.get("g_mip")]
    if rat:
        print(f"dt-fit / MIP-match gain ratio: median {np.median(rat):.3f}, rms {np.std(rat):.3f}, n={len(rat)}")
    if a.png:
        c = ROOT.TCanvas("c", "", 1600, 1000)
        keep = []
        for per in periods:
            for pb, _ in PAIRS:
                c.Clear()
                c.Divide(4, 3)
                for b in range(1, 12):
                    if (per, pb, b) not in plots:
                        continue
                    c.cd(b)
                    ROOT.gPad.SetLogz()
                    hbs, pts, r = plots[(per, pb, b)]
                    hbs.SetMinimum(0.5)
                    hbs.GetXaxis().SetRangeUser(0, 9)
                    hbs.Draw("colz")
                    if pts:
                        g = ROOT.TGraph(len(pts))
                        for i, q in enumerate(pts):
                            g.SetPoint(i, q[0], q[1])
                        g.SetMarkerStyle(20)
                        g.SetMarkerSize(0.4)
                        g.Draw("P same")
                        keep.append(g)
                    if r:
                        x = np.linspace(0.5, 9, 400)
                        y = np.nan_to_num(model(x, r["g"], r["d"]))
                        gm = ROOT.TGraph(len(x), x, y.astype(float))
                        gm.SetLineColor(ROOT.kRed)
                        gm.SetLineWidth(2)
                        gm.Draw("L same")
                        keep.append(gm)
                c.Print(f"{a.prefix}_{a.var}_{per}_dt_{pb}.png")


if __name__ == "__main__":
    main()
