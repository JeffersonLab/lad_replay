#!/usr/bin/env python3
"""Fit the LAD energy calibration (proton scale) from lad_edep_histos.C output.

For every PMT-HV period, plane and paddle, the proton ridge in (1/beta, ADC) is extracted slice by
slice and fit with
    ADC(1/beta) = g * Lnorm(T(1/beta + d))
where T is the proton kinetic energy at the vertex from 1/beta, Lnorm the Birks light in the bar for
that proton (after the material in front of the bar; for back planes that includes the front bar),
normalised so that a proton that just punches through reads its deposited energy (the punch-through
energy, 80.4 MeV at normal incidence). d absorbs ToF offsets. Free per bar: g, d.
The calibration constant is 1/g (MeV per ADC unit). Works for plane 200, which has no back plane.

usage: lad_edep_fit.py histos.root out_prefix [--var int|amp] [--kB 0.0126] [--pre 0.9]
writes <out_prefix>_<var>.csv, <out_prefix>_<var>.pdf and, per period, <out_prefix>_<period>.param
(both int and amp constants once both vars are fit; see --param).
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

PLANES = ["000", "001", "100", "101", "200"]
BACK = {"001", "101"}
GAP_AIR_CM = 40.0  # front-back plane separation


def model_table(back, kB, pre, path=1.0):
    """Lnorm(T_vertex) on a grid; returns (T, L)."""
    Tv = np.arange(1.0, 600.0, 0.5)
    T = M.t_out(Tv, pre)
    if back:
        T = M.t_out(T, M.T_BAR + GAP_AIR_CM * 1.205e-3 * 0.92)
    t = path * M.T_BAR
    L = M.light(T, t, kB) if kB > 0 else M.edep(T, t)
    Ept = M.punch_through_T(t)
    Lpt = float(M.light(np.array([Ept]), t, kB)[0]) if kB > 0 else Ept
    return Tv, L * Ept / Lpt


class Model:
    def __init__(self, back, kB, pre, path=1.0):
        self.Tv, self.L = model_table(back, kB, pre, path)
        self.ib_grid = M.invbeta_from_t(self.Tv)  # decreasing with T
        # cusp: T_vertex of the maximum
        i = int(np.argmax(self.L))
        self.ib_cusp = float(self.ib_grid[i])
        self.Lmax = float(self.L[i])

    def __call__(self, ib, g, d):
        ib = np.asarray(ib, dtype=float) + d
        # interpolate in ib (grid decreasing)
        return g * np.interp(ib, self.ib_grid[::-1], self.L[::-1], left=self.L[-1], right=0.0)


ACC_WIN = (-1.8, -0.3)  # 1/beta window before the photon flash: accidentals only


def acc_spectrum(h2):
    """ADC spectrum of accidentals per 1/beta bin."""
    ax = h2.GetXaxis()
    b1, b2 = ax.FindBin(ACC_WIN[0]), ax.FindBin(ACC_WIN[1])
    acc = h2.ProjectionY("_acc", b1, b2)
    acc.Scale(1.0 / (b2 - b1 + 1))
    acc.SetDirectory(0)
    return acc


def slice_peak(h2, ib_lo, ib_hi, center, frac_lo=0.6, frac_hi=1.5, nmin=80, acc=None):
    """Peak of the (accidental-subtracted) ADC distribution in [ib_lo, ib_hi], searched in
    [frac_lo, frac_hi] * center."""
    ax = h2.GetXaxis()
    b1, b2 = ax.FindBin(ib_lo + 1e-9), ax.FindBin(ib_hi - 1e-9)
    py = h2.ProjectionY("_py", b1, b2)
    if acc is not None:
        py.Add(acc, -(b2 - b1 + 1))
    lo, hi = frac_lo * center, frac_hi * center
    ya = py.GetXaxis()
    k1, k2 = ya.FindBin(lo), ya.FindBin(hi)
    if py.Integral(k1, k2) < nmin:
        return None
    py.GetXaxis().SetRange(k1, k2)
    mu = py.GetBinCenter(py.GetMaximumBin())
    s = 0.12 * mu
    f = ROOT.TF1("g", "gaus", mu - s, mu + s)
    res = None
    for _ in range(3):
        f.SetRange(mu - 1.0 * s, mu + 1.0 * s)
        f.SetParameters(py.GetMaximum(), mu, s)
        r = py.Fit(f, "QNRS")
        if int(r) != 0:
            break
        m2, s2, e2 = f.GetParameter(1), abs(f.GetParameter(2)), f.GetParError(1)
        if not (lo < m2 < hi) or s2 <= 0 or s2 > mu:
            break
        mu, s, res = m2, s2, (m2, max(e2, 1e-3 * m2), s2, py.Integral(k1, k2))
    return res


def fit_bar(h2, model, ib_range, cusp_excl, g0, acc=None):
    """Iteratively extract the ridge around the model and fit (g, d)."""
    g, d = g0, 0.0
    pts = []
    for it in range(4):
        pts = []
        ib = ib_range[0]
        while ib < ib_range[1]:
            w = 0.04 if ib < 2.5 else 0.08
            c = float(model(np.array([ib + w / 2]), g, d)[0])
            if c > 0:
                r = slice_peak(h2, ib, ib + w, c, acc=acc)
                if r is not None:
                    pts.append((ib + w / 2, r[0], r[1], r[2], r[3]))
            ib += w
        if len(pts) < 6:
            return None, pts
        P = np.array(pts)
        use = np.abs(P[:, 0] + d - model.ib_cusp) > cusp_excl
        if use.sum() < 5:
            return None, pts
        x, y = P[use, 0], P[use, 1]
        # errors: peak error plus a 2% model systematic so a few slices with tiny errors don't dominate
        e = np.sqrt(P[use, 2] ** 2 + (0.02 * y) ** 2)
        res = least_squares(lambda p: (model(x, p[0], p[1]) - y) / e, [g, d], bounds=([0.01 * g0, -0.5], [100 * g0, 0.5]))
        g, d = res.x
        chi2 = float(np.sum(res.fun**2))
    return dict(g=g, d=d, chi2=chi2, ndf=int(use.sum()) - 2, npts=len(pts)), pts


def scan_gain(h2, model, ib_range, acc, frac=0.12):
    """Initial gain: the g whose model curve (+-frac) collects the most accidental-subtracted counts."""
    ax, ay = h2.GetXaxis(), h2.GetYaxis()
    sub = np.zeros((h2.GetNbinsX(), h2.GetNbinsY()))
    accv = np.array([acc.GetBinContent(j + 1) for j in range(h2.GetNbinsY())])
    for i in range(h2.GetNbinsX()):
        for j in range(h2.GetNbinsY()):
            sub[i, j] = h2.GetBinContent(i + 1, j + 1) - accv[j]
    xc = np.array([ax.GetBinCenter(i + 1) for i in range(h2.GetNbinsX())])
    yc = np.array([ay.GetBinCenter(j + 1) for j in range(h2.GetNbinsY())])
    xs = (xc > ib_range[0]) & (xc < ib_range[1])
    best, gbest = -1e30, 1.0
    for g in np.exp(np.linspace(np.log(0.03), np.log(30), 160)):
        m = model(xc[xs], g, 0.0)
        ok = m > 3 * (yc[1] - yc[0])
        tot = 0.0
        for k, i in enumerate(np.where(xs)[0]):
            if not ok[k] or m[k] > yc[-1]:
                continue
            sel = np.abs(yc - m[k]) < frac * m[k]
            tot += sub[i, sel].sum()
        if tot > best:
            best, gbest = tot, g
    return gbest


class DEEModel:
    """Front-back Delta-E - E: back-bar light vs front-bar light (both proton-scale MeV) for protons
    that cross the front bar and stop in the back bar."""

    def __init__(self, kB):
        Tf = np.arange(70.0, 200.0, 0.25)  # kinetic energy entering the front bar
        t = M.T_BAR
        Ept = M.punch_through_T(t)
        Lpt = float(M.light(np.array([Ept]), t, kB)[0]) if kB > 0 else Ept
        norm = Ept / Lpt
        Ef = (M.light(Tf, t, kB) if kB > 0 else M.edep(Tf, t)) * norm
        Tb = M.t_out(M.t_out(Tf, t), GAP_AIR_CM * 1.205e-3 * 0.92)
        Eb = (M.light(Tb, t, kB) if kB > 0 else M.edep(Tb, t)) * norm
        # stopping-in-back branch: Tb between ~0 and the back punch-through, i.e. up to the Eb maximum
        imax = int(np.argmax(Eb))
        sel = (Tb > 2.0) & (np.arange(len(Tf)) <= imax)
        self.Ef, self.Eb = Ef[sel], Eb[sel]  # Ef decreasing, Eb increasing along the branch
        self.Ef_cusp, self.Eb_max = float(Ef[imax]), float(Eb[imax])

    def __call__(self, ef):
        return np.interp(ef, self.Ef[::-1], self.Eb[::-1], left=np.nan, right=np.nan)


def fit_dee(h2, g_front, dee, g0, nmin=60):
    """Back gain from the Delta-E - E ridge, given the front gain (ADC per MeV)."""
    ax = h2.GetXaxis()
    gb = g0
    pts = []
    for it in range(4):
        pts = []
        ef_lo, ef_hi = dee.Ef_cusp + 3, min(dee.Ef.max() - 3, 75.0)
        for ef in np.arange(ef_lo, ef_hi, 2.0):
            a1, a2 = ef * g_front, (ef + 2.0) * g_front
            c = float(dee(np.array([ef + 1.0]))[0]) * gb
            if not np.isfinite(c) or c <= 0:
                continue
            py = h2.ProjectionY("_dee", ax.FindBin(a1), ax.FindBin(a2))
            r = None
            if py.GetEntries() >= nmin:
                lo, hi = 0.6 * c, 1.5 * c
                k1, k2 = py.GetXaxis().FindBin(lo), py.GetXaxis().FindBin(hi)
                if py.Integral(k1, k2) >= nmin:
                    py.GetXaxis().SetRange(k1, k2)
                    mu = py.GetBinCenter(py.GetMaximumBin())
                    sg = 0.15 * mu
                    f = ROOT.TF1("gd", "gaus", mu - sg, mu + sg)
                    f.SetParameters(py.GetMaximum(), mu, sg)
                    rr = py.Fit(f, "QNRS")
                    if int(rr) == 0 and lo < f.GetParameter(1) < hi:
                        r = (f.GetParameter(1), max(f.GetParError(1), 1e-3 * mu))
            if r:
                pts.append((ef + 1.0, r[0], r[1]))
        if len(pts) < 4:
            return None, pts
        P = np.array(pts)
        pred = dee(P[:, 0])
        e = np.sqrt(P[:, 2] ** 2 + (0.02 * P[:, 1]) ** 2)
        w = 1 / e**2
        gb = float(np.sum(w * pred * P[:, 1]) / np.sum(w * pred * pred))  # linear least squares
        chi2 = float(np.sum(((gb * pred - P[:, 1]) / e) ** 2))
    return dict(g=gb, chi2=chi2, ndf=len(pts) - 1, npts=len(pts)), pts


def mip_peak(h2, ib_lo=0.95, ib_hi=1.12, acc=None):
    ax = h2.GetXaxis()
    b1, b2 = ax.FindBin(ib_lo), ax.FindBin(ib_hi)
    py = h2.ProjectionY("_mip", b1, b2)
    if acc is not None:
        py.Add(acc, -(b2 - b1 + 1))
    if py.GetEntries() < 200:
        return None
    # skip the threshold region: start the search above the 30% quantile
    q = np.zeros(1)
    py.GetQuantiles(1, q, np.array([0.3]))
    py.GetXaxis().SetRangeUser(q[0], py.GetXaxis().GetXmax())
    mu = py.GetBinCenter(py.GetMaximumBin())
    f = ROOT.TF1("lan", "landau", 0.7 * mu, 2.0 * mu)
    f.SetParameters(py.GetMaximum(), mu, 0.1 * mu)
    r = py.Fit(f, "QNRS")
    if int(r) != 0:
        return mu
    return f.GetParameter(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histos")
    ap.add_argument("prefix")
    ap.add_argument("--var", default="int", choices=["int", "amp"])
    ap.add_argument("--kB", type=float, default=0.0126 * M.RHO_PVT, help="Birks kB in g/(cm2 MeV); 0 = no quenching")
    ap.add_argument("--pre", type=float, default=0.9, help="material before the (front) bar, g/cm2 PVT-equivalent")
    ap.add_argument("--periods", default="", help="comma list of period dirs (default all)")
    ap.add_argument("--cusp-excl", type=float, default=0.08)
    ap.add_argument("--png", action="store_true", help="also write one PNG per page")
    ap.add_argument("--cat", default="trk", choices=["all", "pair", "trk"], help="hit category (lad_edep_histos.C)")
    ap.add_argument("--cat-back", default="pair", choices=["all", "pair", "trk"], help="hit category for back planes")
    a = ap.parse_args()

    fin = ROOT.TFile.Open(a.histos)
    periods = [k.GetName() for k in fin.GetListOfKeys() if k.ReadObj().InheritsFrom("TDirectory")]
    if a.periods:
        periods = [p for p in periods if p in a.periods.split(",")]
    models = {False: Model(False, a.kB, a.pre), True: Model(True, a.kB, a.pre)}
    print("cusp 1/beta: front %.3f (Lmax %.1f), back %.3f (Lmax %.1f)" % (
        models[False].ib_cusp, models[False].Lmax, models[True].ib_cusp, models[True].Lmax))

    rows, allpts = [], {}
    for per in periods:
        d = fin.Get(per)
        for pl in PLANES:
            back = pl in BACK
            mdl = models[back]
            ib_range = (1.25, 3.6) if not back else (1.25, 2.6)
            for b in range(1, 12):
                cat = a.cat_back if back else a.cat
                h = d.Get(f"h_{a.var}_{cat}_{pl}_{b}")
                if not h or h.GetEntries() < 2000:
                    rows.append(dict(period=per, plane=pl, paddle=b, entries=int(h.GetEntries()) if h else 0, ok=0))
                    continue
                acc = acc_spectrum(h)
                g0 = scan_gain(h, mdl, ib_range, acc)
                r, pts = fit_bar(h, mdl, ib_range, a.cusp_excl, g0, acc)
                allpts[(per, pl, b)] = (pts, r)
                mip = mip_peak(h, acc=acc)
                row = dict(period=per, plane=pl, paddle=b, cat=cat, entries=int(h.GetEntries()), ok=int(r is not None), g0=g0)
                if r:
                    row.update(g=r["g"], d=r["d"], chi2ndf=r["chi2"] / max(r["ndf"], 1), npts=r["npts"],
                               MeV_per_unit=1.0 / r["g"], adc_pt=r["g"] * mdl.Lmax,
                               mip_adc=mip if mip else float("nan"),
                               mip_MeV=(mip / r["g"]) if mip else float("nan"))
                rows.append(row)
                if r:
                    print(f"{per} {pl} {b:2d}: g={r['g']:.4f} ({a.var}/MeV)  1/g={1/r['g']:.3f}  d={r['d']:+.3f}  "
                          f"chi2/n={r['chi2']/max(r['ndf'],1):.1f}  mip={row['mip_MeV']:.1f} MeV")
                else:
                    print(f"{per} {pl} {b:2d}: fit failed ({len(pts)} pts)")

    # back planes: Delta-E - E against the same paddle of the front plane (primary method for 001/101)
    dee = DEEModel(a.kB)
    print("Delta-E - E: front light at the back-bar cusp %.1f MeV, back max %.1f MeV" % (dee.Ef_cusp, dee.Eb_max))
    byk = {(r["period"], r["plane"], r["paddle"]): r for r in rows}
    deepts = {}
    for per in periods:
        d = fin.Get(per)
        for pl, pf in [("001", "000"), ("101", "100")]:
            for b in range(1, 12):
                rf = byk.get((per, pf, b))
                rb = byk.get((per, pl, b))
                h = d.Get(f"h_dee_{a.var}_{pl}_{b}")
                if rb is None or not h or h.GetEntries() < 500 or not rf or not rf.get("g"):
                    continue
                g0 = rb["g"] if rb.get("g") else rf["g"]
                r, pts = fit_dee(h, rf["g"], dee, g0)
                deepts[(per, pl, b)] = (pts, r, rf["g"])
                if r:
                    rb.update(dee_g=r["g"], dee_chi2ndf=r["chi2"] / max(r["ndf"], 1), dee_npts=r["npts"])
                    print(f"{per} {pl} {b:2d}: dE-E g={r['g']:.4f} (ToF g={rb.get('g', float('nan')):.4f})  chi2/n={r['chi2']/max(r['ndf'],1):.1f}")
                else:
                    print(f"{per} {pl} {b:2d}: dE-E failed ({len(pts)} pts)")
    for r in rows:
        if r["plane"] in BACK and r.get("dee_g"):
            r["g_final"], r["method"] = r["dee_g"], "dEE"
        elif r.get("g"):
            r["g_final"], r["method"] = r["g"], "tof"
        if r.get("g_final"):
            r["MeV_per_unit_final"] = 1.0 / r["g_final"]

    keys = sorted({k for r in rows for k in r.keys()}, key=lambda k: list(rows[0].keys()).index(k) if k in rows[0] else 99)
    with open(f"{a.prefix}_{a.var}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["period", "plane", "paddle", "cat", "entries", "ok", "g0", "g", "d", "chi2ndf", "npts",
                                          "MeV_per_unit", "adc_pt", "mip_adc", "mip_MeV",
                                          "dee_g", "dee_chi2ndf", "dee_npts", "g_final", "method", "MeV_per_unit_final"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    # plots
    c = ROOT.TCanvas("c", "", 1600, 1000)
    pdf = f"{a.prefix}_{a.var}.pdf"
    c.Print(pdf + "[")
    keep = []
    for per in periods:
        d = fin.Get(per)
        for pl in PLANES:
            mdl = models[pl in BACK]
            c.Clear()
            c.Divide(4, 3)
            for b in range(1, 12):
                c.cd(b)
                ROOT.gPad.SetLogz()
                cat = a.cat_back if pl in BACK else a.cat
                h = d.Get(f"h_{a.var}_{cat}_{pl}_{b}")
                if not h:
                    continue
                acc = acc_spectrum(h)
                h = h.Clone(h.GetName() + "_sub")
                for ix in range(1, h.GetNbinsX() + 1):
                    for iy in range(1, h.GetNbinsY() + 1):
                        h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) - acc.GetBinContent(iy))
                h.SetMinimum(0)
                keep.append(h)
                h.GetXaxis().SetRangeUser(0.9, 4.0)
                h.Draw("colz")
                pts, r = allpts.get((per, pl, b), ([], None))
                if pts:
                    g = ROOT.TGraph(len(pts))
                    for i, p in enumerate(pts):
                        g.SetPoint(i, p[0], p[1])
                    g.SetMarkerStyle(20)
                    g.SetMarkerSize(0.4)
                    g.Draw("P same")
                    keep.append(g)
                if r:
                    x = np.linspace(1.0, 4.0, 300)
                    y = mdl(x, r["g"], r["d"])
                    gm = ROOT.TGraph(len(x), x.astype(float), y.astype(float))
                    gm.SetLineColor(ROOT.kRed)
                    gm.SetLineWidth(2)
                    gm.Draw("L same")
                    keep.append(gm)
                    t = ROOT.TLatex()
                    t.SetNDC()
                    t.SetTextSize(0.05)
                    t.DrawLatex(0.4, 0.85, f"1/g={1/r['g']:.3f} d={r['d']:+.2f}")
                    keep.append(t)
            c.cd(12)
            t = ROOT.TLatex()
            t.SetTextSize(0.08)
            t.DrawLatexNDC(0.05, 0.5, f"{per} plane {pl} ({a.var})")
            keep.append(t)
            c.Print(pdf)
            if a.png:
                c.Print(f"{a.prefix}_{a.var}_{per}_{pl}.png")
    for per in periods:
        d = fin.Get(per)
        for pl in ["001", "101"]:
            c.Clear()
            c.Divide(4, 3)
            for b in range(1, 12):
                c.cd(b)
                ROOT.gPad.SetLogz()
                h = d.Get(f"h_dee_{a.var}_{pl}_{b}")
                if not h:
                    continue
                h.Draw("colz")
                pts, r, gf = deepts.get((per, pl, b), ([], None, None))
                if pts:
                    g = ROOT.TGraph(len(pts))
                    for i, q in enumerate(pts):
                        g.SetPoint(i, q[0] * gf, q[1])
                    g.SetMarkerStyle(20)
                    g.SetMarkerSize(0.4)
                    g.Draw("P same")
                    keep.append(g)
                if r:
                    ef = dee.Ef
                    gm = ROOT.TGraph(len(ef), (ef * gf).astype(float), (dee.Eb * r["g"]).astype(float))
                    gm.SetLineColor(ROOT.kRed)
                    gm.SetLineWidth(2)
                    gm.Draw("L same")
                    keep.append(gm)
            c.cd(12)
            t = ROOT.TLatex()
            t.SetTextSize(0.08)
            t.DrawLatexNDC(0.05, 0.5, f"{per} dE-E back {pl} ({a.var})")
            keep.append(t)
            c.Print(pdf)
            if a.png:
                c.Print(f"{a.prefix}_{a.var}_{per}_dee_{pl}.png")
    c.Print(pdf + "]")


if __name__ == "__main__":
    main()
