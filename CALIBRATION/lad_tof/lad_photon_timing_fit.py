#!/usr/bin/env python3
"""Fit the LAD photon peak from lad_photon_timing.C output and derive lglobal_time_offset.

usage: python3 lad_photon_timing_fit.py <in.root> <outdir> [--ecut 5,10,20]

The ToF offsets the histograms were made with are read from the TNamed
"tof_offsets" that lad_photon_timing.C writes (files without it predate the
calibration and used the shared -1710 ns); --offset-p/--offset-h override.
Suggested offset = offset used - fitted photon tof-L/c.

For every histogram the accidental background is removed with a template made
by folding the region before the photon peak ([mu0-SB_LO, mu0-SB_HI]) on the RF
period, then a Gaussian is fit iteratively to the peak core ([mu-1.5, mu+1] sigma; slower particles sit on the late side).
The FWHM is also measured directly on the background-subtracted histogram.

Outputs in <outdir>:
  photon_timing_plots.root   canvases (per bar / per plane / total, tof and tof-L/c, all cuts)
  photon_timing_plots.pdf    the same, one canvas per page
  fits.csv                   one row per (spec, cat, var, plane, bar) fit
  summary.txt                widths, peak positions and the offset recommendation
"""
import argparse
import csv
import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)
ROOT.TH1.AddDirectory(False)
KEEP = []  # keep drawn objects alive until the canvases are written

RF_PERIOD = 4.00801
SB_LO, SB_HI = 45.0, 8.0  # sideband = [mu0-45, mu0-8] ns (before the photon peak)
FIT_NSIG_LO, FIT_NSIG_HI = 1.5, 1.0  # asymmetric core window: slower particles sit on the late side
C_LIGHT = 29.9792458
PLANES = ["000", "001", "100", "101", "200"]
ZPOS = [618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935]
THETA = [2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772]
SPECS = ["P", "H"]
# run periods = row boundaries of l_rf_offset in lladkine.param
PERIODS = [22013, 22560, 22661, 22745, 23597, 24000]
SPEC_NAME = {"P": "SHMS", "H": "HMS"}

COL = {k: ROOT.TColor.GetColor(v) for k, v in {
    "all": "#2a78d6", "back": "#2a78d6", "back_fv": "#eb6834",
    "fv_lo": "#1baf7a", "fv_hi": "#4a3aa7"}.items()}
LABEL = {"all": "all hits", "back": "back-plane hits", "back_fv": "back + front veto",
         "fv_lo": "back + veto, E<{e} MeV", "fv_hi": "back + veto, E>{e} MeV"}


def expected_Lc(p, b, y=0.0, z=0.0):
    x0, z0 = 110.0 - 22.0 * b, ZPOS[p]
    c, s = math.cos(THETA[p]), math.sin(THETA[p])
    X, Z = x0 * c + z0 * s, -x0 * s + z0 * c - z
    return math.sqrt(X * X + y * y + Z * Z) / C_LIGHT


def bgsub(h, mu0, period=RF_PERIOD):
    """Periodic background template from the pre-peak sideband, subtracted everywhere."""
    lo, hi = mu0 - SB_LO, mu0 - SB_HI
    nph = max(1, int(round(period / h.GetBinWidth(1))))
    tmpl, cnt = [0.0] * nph, [0] * nph
    for i in range(1, h.GetNbinsX() + 1):
        x = h.GetBinCenter(i)
        if lo <= x < hi:
            k = int(((x % period) / period) * nph) % nph
            tmpl[k] += h.GetBinContent(i)
            cnt[k] += 1
    tmpl = [t / c if c else 0.0 for t, c in zip(tmpl, cnt)]
    out = h.Clone(h.GetName() + "_bgsub")
    for i in range(1, h.GetNbinsX() + 1):
        k = int(((h.GetBinCenter(i) % period) / period) * nph) % nph
        out.SetBinContent(i, h.GetBinContent(i) - tmpl[k])
        out.SetBinError(i, math.sqrt(max(h.GetBinContent(i), 0.0) + tmpl[k] / max(1, cnt[k])))
    bg_level = sum(tmpl) / nph
    return out, bg_level


def fwhm(h, mu, win=4.0):
    """Direct FWHM on a (bg-subtracted) histogram around mu, linear interpolation."""
    b0, b1 = h.FindBin(mu - win), h.FindBin(mu + win)
    ib = max(range(b0, b1 + 1), key=lambda i: h.GetBinContent(i))
    half = h.GetBinContent(ib) / 2.0
    if half <= 0:
        return float("nan")
    i = ib
    while i > 1 and h.GetBinContent(i) > half:
        i -= 1
    xl = h.GetBinCenter(i) + (half - h.GetBinContent(i)) / max(h.GetBinContent(i + 1) - h.GetBinContent(i), 1e-9) * h.GetBinWidth(1)
    j = ib
    while j < h.GetNbinsX() and h.GetBinContent(j) > half:
        j += 1
    xr = h.GetBinCenter(j) - (half - h.GetBinContent(j)) / max(h.GetBinContent(j - 1) - h.GetBinContent(j), 1e-9) * h.GetBinWidth(1)
    return xr - xl


def fit_peak(h, mu0, search=5.0, min_counts=50):
    """Return dict(mu, emu, sigma, esigma, fwhm, nsig, sb) and the bg-subtracted clone."""
    res = dict(mu=float("nan"), emu=float("nan"), sigma=float("nan"), esigma=float("nan"),
               fwhm=float("nan"), nsig=0.0, sb=float("nan"), ok=False)
    if h.GetEntries() < min_counts:
        return res, None
    hs, bg = bgsub(h, mu0)
    # coarse peak in mu0 +- search on a lightly smoothed copy
    hc = hs.Clone(hs.GetName() + "_c")
    hc.GetXaxis().SetRangeUser(mu0 - search, mu0 + search)
    if hc.Integral() <= 0:
        return res, hs
    reb = max(1, int(round(0.3 / hc.GetBinWidth(1))))
    hr = hc.Clone(hc.GetName() + "_r")
    hr.Rebin(reb)
    hr.GetXaxis().SetRangeUser(mu0 - search, mu0 + search)
    mu, sig = hr.GetBinCenter(hr.GetMaximumBin()), 1.0
    f = ROOT.TF1("g_" + h.GetName(), "gaus", mu - 2, mu + 2)
    for _ in range(6):
        f.SetRange(mu - FIT_NSIG_LO * sig, mu + FIT_NSIG_HI * sig)
        f.SetParameters(max(hs.GetBinContent(hs.FindBin(mu)), 1.0), mu, sig)
        r = hs.Fit(f, "QNRS0")
        if int(r) != 0 or f.GetParameter(2) <= 0:
            break
        mu_new, sig_new = f.GetParameter(1), abs(f.GetParameter(2))
        if abs(mu_new - mu) > search or sig_new > 5 or sig_new < 0.05:
            break
        mu, sig = mu_new, max(sig_new, 2 * hs.GetBinWidth(1))
        res.update(mu=mu, emu=f.GetParError(1), sigma=sig, esigma=f.GetParError(2), ok=True)
    if res["ok"]:
        b0, b1 = hs.FindBin(mu - 2 * sig), hs.FindBin(mu + 2 * sig)
        nsig = sum(hs.GetBinContent(i) for i in range(b0, b1 + 1))
        nbg = bg * (b1 - b0 + 1)
        res.update(nsig=nsig, sb=nsig / nbg if nbg > 0 else float("inf"), fwhm=fwhm(hs, mu))
        f.SetRange(mu - FIT_NSIG_LO * sig, mu + FIT_NSIG_HI * sig)
        f.SetLineColor(ROOT.kBlack)
        f.SetLineWidth(2)
        hs.GetListOfFunctions().Add(f.Clone())
    return res, hs


def get(fin, path):
    o = fin.Get(path)
    if not o:
        raise KeyError(path)
    return o


def proj_e(h2, emin=None, emax=None, name=None):
    """Project tofc-vs-edep 2D onto tofc for emin <= edep < emax (MeV)."""
    ay = h2.GetYaxis()
    b0 = 1 if emin is None else ay.FindBin(emin + 1e-6)
    b1 = ay.GetNbins() if emax is None else ay.FindBin(emax - 1e-6)
    return h2.ProjectionX(name or (h2.GetName() + "_e%s_%s" % (emin, emax)), b0, b1)


def style(h, key, rebin=1):
    h = h.Clone()
    if rebin > 1:
        h.Rebin(rebin)
    h.SetStats(0)
    h.SetLineColor(COL[key])
    h.SetLineWidth(2)
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outdir")
    ap.add_argument("--offset-p", type=float, default=None, help="SHMS ToF offset the input was made with (ns)")
    ap.add_argument("--offset-h", type=float, default=None, help="HMS ToF offset the input was made with (ns)")
    ap.add_argument("--ecut", default="2,5,10,15,20", help="back-plane edep cuts to scan (MeV)")
    ap.add_argument("--eplot", type=float, default=10.0, help="edep split (MeV) used in the overlay plots")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fin = ROOT.TFile.Open(a.infile)
    off = {"P": -1710.0, "H": -1710.0}  # files without "tof_offsets" predate the per-spectrometer calibration
    tn = fin.Get("tof_offsets")
    if tn:
        off.update({k: float(v) for k, v in (kv.split("=") for kv in tn.GetTitle().split(","))})
    if a.offset_p is not None:
        off["P"] = a.offset_p
    if a.offset_h is not None:
        off["H"] = a.offset_h
    ecuts = [float(x) for x in a.ecut.split(",")]

    rows = []
    summary = []

    def record(sp, cat, var, plane, bar, res):
        rows.append(dict(spec=sp, cat=cat, var=var, plane=plane, bar=bar,
                         **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()}))

    # ---- coarse photon-peak position per spec and var (from all-hits total) ----
    mu0 = {}
    for sp in SPECS:
        for var in ["tofc", "tofcrf"]:
            h = get(fin, f"{sp}/all/{var}_total").Clone()
            h.GetXaxis().SetRangeUser(-20, 80)
            h.Rebin(5)
            h.GetXaxis().SetRangeUser(-20, 80)
            mu0[(sp, var)] = h.GetBinCenter(h.GetMaximumBin())
        mu0[(sp, "tof")] = mu0[(sp, "tofc")] + 20.0  # tof window centre is refined per bar below

    # ---- totals / planes / bars for 1D cats ----
    fits = {}
    for sp in SPECS:
        for cat in ["all", "back", "back_fv"]:
            for var in ["tofc", "tofcrf"]:
                m0 = mu0[(sp, var)]
                r, _ = fit_peak(get(fin, f"{sp}/{cat}/{var}_total"), m0)
                fits[(sp, cat, var, "total", -1)] = r
                record(sp, cat, var, "total", -1, r)
                for p in range(5):
                    if cat != "all" and p not in (1, 3):
                        continue
                    r, _ = fit_peak(get(fin, f"{sp}/{cat}/{var}_p{p}_sum"), m0)
                    fits[(sp, cat, var, PLANES[p], -1)] = r
                    record(sp, cat, var, PLANES[p], -1, r)
                    for b in range(11):
                        r, _ = fit_peak(get(fin, f"{sp}/{cat}/{var}_p{p}_b{b}"), m0, min_counts=200)
                        fits[(sp, cat, var, PLANES[p], b)] = r
                        record(sp, cat, var, PLANES[p], b, r)
            # raw tof per bar: peak should sit at L/c once calibrated
            for p in range(5):
                if cat != "all" and p not in (1, 3):
                    continue
                for b in range(11):
                    m0 = mu0[(sp, "tofc")] + expected_Lc(p, b)
                    r, _ = fit_peak(get(fin, f"{sp}/{cat}/tof_p{p}_b{b}"), m0, min_counts=200)
                    r["expected_Lc"] = expected_Lc(p, b)
                    fits[(sp, cat, "tof", PLANES[p], b)] = r
                    record(sp, cat, "tof", PLANES[p], b, r)

    # ---- energy-cut scan on back-plane hits (with and without front veto) ----
    escan = {}
    for sp in SPECS:
        for cat in ["back", "back_fv"]:
            for var in ["tofc", "tofcrf"]:
                h2 = get(fin, f"{sp}/{cat}/{var}_vs_edep_total")
                m0 = mu0[(sp, var)]
                for e in ecuts:
                    for side, (lo, hi) in (("lo", (None, e)), ("hi", (e, None))):
                        r, _ = fit_peak(proj_e(h2, lo, hi), m0)
                        escan[(sp, cat, var, side, e)] = r
                        record(sp, f"{cat}_E{side}{e:g}", var, "total", -1, r)
    eplot = a.eplot

    with open(os.path.join(a.outdir, "fits.csv"), "w", newline="") as fo:
        keys = ["spec", "cat", "var", "plane", "bar", "mu", "emu", "sigma", "esigma", "fwhm", "nsig", "sb", "ok", "expected_Lc"]
        w = csv.DictWriter(fo, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- per-run peak position ----
    runfits = {}
    perfits = {}
    aligned = {}
    for sp in SPECS:
        for var in ["tofc", "tofcrf"]:
            h2 = get(fin, f"{sp}/run/{var}_vs_run_back_fv")
            pts = []
            for ib in range(1, h2.GetNbinsX() + 1):
                hp = h2.ProjectionY(f"_r{sp}{var}{ib}", ib, ib)
                if hp.GetEntries() < 1500:
                    continue
                r, _ = fit_peak(hp, mu0[(sp, var)], min_counts=1500)
                if r["ok"] and r["emu"] < 0.4:
                    pts.append((h2.GetXaxis().GetBinLowEdge(ib), r["mu"], r["emu"], r["sigma"]))
            runfits[(sp, var)] = pts
            for lo, hi in zip(PERIODS[:-1], PERIODS[1:]):
                hp = h2.ProjectionY(f"_per{sp}{var}{lo}", h2.GetXaxis().FindBin(lo), h2.GetXaxis().FindBin(hi - 1))
                r, _ = fit_peak(hp, mu0[(sp, var)], min_counts=500)
                r["n"] = hp.GetEntries()
                perfits[(sp, var, lo)] = r
            # Effect of the run-period shift on the combined width (back + veto):
            # all runs as-is / without runs 22560-22660 / each period moved onto the
            # peak of the largest period (whole 0.2 ns bins).
            hall = hexcl = hal = None
            ref = max((perfits[(sp, var, lo)] for lo in PERIODS[:-1] if perfits[(sp, var, lo)]["ok"]),
                      key=lambda r: r["n"], default=None)
            for lo, hi in zip(PERIODS[:-1], PERIODS[1:]):
                hp = h2.ProjectionY(f"_pa{sp}{var}{lo}", h2.GetXaxis().FindBin(lo), h2.GetXaxis().FindBin(hi - 1))
                r = perfits[(sp, var, lo)]
                if hall is None:
                    hall = hp.Clone(f"_pall{sp}{var}")
                else:
                    hall.Add(hp)
                if lo >= 22745:
                    if hexcl is None:
                        hexcl = hp.Clone(f"_pexcl{sp}{var}")
                    else:
                        hexcl.Add(hp)
                hs = hp.Clone(f"_psh{sp}{var}{lo}")
                hs.Reset()
                k = int(round((r["mu"] - ref["mu"]) / hp.GetBinWidth(1))) if (r["ok"] and ref) else 0
                for b in range(1, hp.GetNbinsX() + 1):
                    if 1 <= b - k <= hp.GetNbinsX():
                        hs.SetBinContent(b - k, hp.GetBinContent(b))
                if hal is None:
                    hal = hs
                else:
                    hal.Add(hs)
            for tag, h in (("all", hall), ("excl", hexcl), ("aligned", hal)):
                aligned[(sp, var, tag)] = fit_peak(h, mu0[(sp, var)])[0] if h else dict(ok=False)

    # =================== plots ===================
    pdf = os.path.join(a.outdir, "photon_timing_plots.pdf")
    fout = ROOT.TFile(os.path.join(a.outdir, "photon_timing_plots.root"), "RECREATE")
    canv_list = []

    def save(c, d):
        dd = fout
        for seg in d.split("/"):
            dd = dd.GetDirectory(seg) or dd.mkdir(seg)
        dd.cd()
        c.Write()
        canv_list.append(c)
        os.makedirs(os.path.join(a.outdir, "png", d), exist_ok=True)
        c.SaveAs(os.path.join(a.outdir, "png", d, c.GetName() + ".png"))

    def overlay(pad, hists, title, xr, logy=False, fitres=None):
        pad.cd()
        if logy:
            pad.SetLogy()
        ymax = max(h.GetMaximum() for _, h in hists) * (3 if logy else 1.15)
        leg = ROOT.TLegend(0.52, 0.62, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        for i, (lab, h) in enumerate(hists):
            h.SetTitle(title)
            h.GetXaxis().SetRangeUser(*xr)
            h.SetMaximum(ymax)
            h.SetMinimum(0.5 if logy else 0.0)
            h.Draw("HIST" if i == 0 else "HIST SAME")
            leg.AddEntry(h, lab, "l")
        if fitres:
            for lab, r in fitres:
                if r and r["ok"]:
                    leg.AddEntry(ROOT.nullptr, "%s: #mu=%.2f #sigma=%.2f ns" % (lab, r["mu"], r["sigma"]), "")
        leg.Draw()
        KEEP.append((hists, leg))

    def cat_hists(sp, var, key_suffix, e=eplot, rebin=2):
        """(label, hist) list for all, back, back_fv, back_fv E<e, back_fv E>e."""
        if key_suffix == "total":
            n1, n2 = f"{var}_total", f"{var}_vs_edep_total"
        elif key_suffix.startswith("p") and "_b" not in key_suffix:
            n1, n2 = f"{var}_{key_suffix}_sum", f"{var}_vs_edep_{key_suffix}_sum"
        else:
            n1, n2 = f"{var}_{key_suffix}", f"{var}_vs_edep_{key_suffix}"
        out = [(LABEL["all"], style(get(fin, f"{sp}/all/{n1}"), "all", rebin))]
        p = int(key_suffix[1]) if key_suffix != "total" else -1
        if p in (-1, 1, 3):
            out.append((LABEL["back_fv"], style(get(fin, f"{sp}/back_fv/{n1}"), "back_fv", rebin)))
            if var != "tof":
                h2 = get(fin, f"{sp}/back_fv/{n2}")
                rb = max(1, rebin // 2)  # 2D x-bins are 0.2 ns
                out.append((LABEL["fv_lo"].format(e=f"{e:g}"), style(proj_e(h2, None, e), "fv_lo", rb)))
                out.append((LABEL["fv_hi"].format(e=f"{e:g}"), style(proj_e(h2, e, None), "fv_hi", rb)))
        return out

    for sp in SPECS:
        sname = f"{sp} ({SPEC_NAME[sp]} vertex time)"
        for var, xr_fn in (("tof", lambda m: (m - 35, m + 35)), ("tofc", lambda m: (m - 35, m + 35)),
                           ("tofcrf", lambda m: (m - 35, m + 35))):
            m = mu0[(sp, "tofc" if var == "tof" else var)] + (20 if var == "tof" else 0)
            xt = "tof (ns)" if var == "tof" else ("tof - L/c (ns)" if var == "tofc" else "tof_{RF} - L/c (ns)")
            # summary: 5 planes + total
            for logy in (False, True):
                c = ROOT.TCanvas(f"c_{sp}_{var}_summary{'_log' if logy else ''}", f"{sp} {var} summary", 1800, 1100)
                c.Divide(3, 2)
                for p in range(5):
                    hs = cat_hists(sp, var, f"p{p}")
                    fr = None
                    if var != "tof":
                        fr = [("all", fits.get((sp, "all", var, PLANES[p], -1)))]
                        if p in (1, 3):
                            fr.append(("veto", fits.get((sp, "back_fv", var, PLANES[p], -1))))
                    overlay(c.cd(p + 1), hs, f"{sname} plane {PLANES[p]};{xt};hits / 0.2 ns", xr_fn(m), logy, fr)
                hs = cat_hists(sp, var, "total")
                fr = None if var == "tof" else [("all", fits[(sp, "all", var, "total", -1)]),
                                                 ("veto", fits[(sp, "back_fv", var, "total", -1)])]
                overlay(c.cd(6), hs, f"{sname} all planes;{xt};hits / 0.2 ns", xr_fn(m), logy, fr)
                save(c, f"{sp}/{var}")
            # per bar
            for p in range(5):
                c = ROOT.TCanvas(f"c_{sp}_{var}_pl{PLANES[p]}", f"{sp} {var} plane {PLANES[p]} per bar", 1800, 1100)
                c.Divide(4, 3)
                for b in range(11):
                    mb = m + (expected_Lc(p, b) - 20 if var == "tof" else 0)
                    hs = cat_hists(sp, var, f"p{p}_b{b}", rebin=4)
                    fr = None
                    if var != "tof":
                        fr = [("all", fits.get((sp, "all", var, PLANES[p], b)))]
                        if p in (1, 3):
                            fr.append(("veto", fits.get((sp, "back_fv", var, PLANES[p], b))))
                    overlay(c.cd(b + 1), hs, f"{sp} plane {PLANES[p]} bar {b};{xt};hits / 0.4 ns", (mb - 30, mb + 35), False, fr)
                save(c, f"{sp}/{var}")

        # zoomed, background-subtracted photon peak with fits (total), tofc and tofcrf
        c = ROOT.TCanvas(f"c_{sp}_photon_zoom", f"{sp} photon peak (bg-subtracted)", 1800, 1100)
        c.Divide(3, 2)
        ipad = 1
        for var in ("tofc", "tofcrf"):
            m0 = mu0[(sp, var)]
            for cat, src in (("all", get(fin, f"{sp}/all/{var}_total")),
                             ("back_fv", get(fin, f"{sp}/back_fv/{var}_total")),
                             ("fv_lo", proj_e(get(fin, f"{sp}/back_fv/{var}_vs_edep_total"), None, eplot))):
                r, hs = fit_peak(src, m0)
                pad = c.cd(ipad)
                ipad += 1
                if hs is None:
                    continue
                hs = style(hs, cat)
                hs.SetStats(0)
                lab = LABEL[cat].format(e=f"{eplot:g}")
                hs.SetTitle(f"{sp} {lab} [{var}], bg-subtracted;{'tof - L/c' if var == 'tofc' else 'tof_{RF} - L/c'} (ns);hits")
                hs.GetXaxis().SetRangeUser(m0 - 8, m0 + 12)
                hs.Draw("HIST")
                for fobj in hs.GetListOfFunctions():
                    fobj.Draw("SAME")
                t = ROOT.TLatex()
                t.SetNDC()
                t.SetTextSize(0.045)
                if r["ok"]:
                    t.DrawLatex(0.55, 0.82, "#mu = %.3f #pm %.3f ns" % (r["mu"], r["emu"]))
                    t.DrawLatex(0.55, 0.76, "#sigma = %.3f ns" % r["sigma"])
                    t.DrawLatex(0.55, 0.70, "FWHM = %.2f ns" % r["fwhm"])
                    t.DrawLatex(0.55, 0.64, "S/B(#pm2#sigma) = %.2f" % r["sb"])
                KEEP.append((hs, t))
        save(c, f"{sp}")

        # tofc vs edep 2D (back, back_fv)
        c = ROOT.TCanvas(f"c_{sp}_tofc_vs_edep", f"{sp} tof-L/c vs edep", 1800, 700)
        c.Divide(2, 1)
        for i, cat in enumerate(("back", "back_fv")):
            c.cd(i + 1)
            ROOT.gPad.SetLogz()
            ROOT.gPad.SetRightMargin(0.13)
            h2 = get(fin, f"{sp}/{cat}/tofc_vs_edep_total").Clone()
            h2.SetTitle(f"{sp} {LABEL[cat]};tof - L/c (ns);edep (MeV)")
            h2.GetYaxis().SetRangeUser(0, 60)
            h2.Draw("COLZ")
            KEEP.append(h2)
        save(c, f"{sp}")

    # P vs H overlay (all hits and veto, tofc), shifted by fit
    c = ROOT.TCanvas("c_PH_overlay", "P vs H photon peak", 1800, 700)
    c.Divide(2, 1)
    for i, cat in enumerate(("all", "back_fv")):
        pad = c.cd(i + 1)
        hp = get(fin, f"P/{cat}/tofc_total").Clone("hp_" + cat)
        hh = get(fin, f"H/{cat}/tofc_total").Clone("hh_" + cat)
        for h, col in ((hp, "#2a78d6"), (hh, "#eb6834")):
            h.Rebin(2)
            h.SetStats(0)
            h.SetMinimum(0)
            h.SetLineColor(ROOT.TColor.GetColor(col))
            h.SetLineWidth(2)
            h.SetTitle(f"{LABEL[cat]}: P (SHMS) vs H (HMS);tof - L/c (ns);hits / 0.2 ns")
            h.GetXaxis().SetRangeUser(-40, 60)
        hp.SetMaximum(1.15 * max(hp.GetMaximum(), hh.GetMaximum()))
        hp.Draw("HIST")
        hh.Draw("HIST SAME")
        leg = ROOT.TLegend(0.6, 0.75, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.AddEntry(hp, "P (SHMS)", "l")
        leg.AddEntry(hh, "H (HMS)", "l")
        leg.Draw()
        KEEP.append((hp, hh, leg))
    save(c, "PH")

    # width vs category / ecut (bar chart-like graph per spec)
    c = ROOT.TCanvas("c_escan", "photon sigma vs edep cut", 1800, 700)
    c.Divide(2, 1)
    for i, var in enumerate(("tofc", "tofcrf")):
        pad = c.cd(i + 1)
        mg = ROOT.TMultiGraph()
        leg = ROOT.TLegend(0.55, 0.7, 0.89, 0.89)
        leg.SetBorderSize(0)
        for sp, col in (("P", "#2a78d6"), ("H", "#eb6834")):
            for side, mk in (("lo", 20), ("hi", 24)):
                g = ROOT.TGraphErrors()
                for e in ecuts:
                    r = escan[(sp, "back_fv", var, side, e)]
                    if r["ok"]:
                        n = g.GetN()
                        g.SetPoint(n, e, r["sigma"])
                        g.SetPointError(n, 0, r["esigma"])
                g.SetMarkerStyle(mk)
                g.SetMarkerSize(1.4)
                g.SetMarkerColor(ROOT.TColor.GetColor(col))
                g.SetLineColor(ROOT.TColor.GetColor(col))
                mg.Add(g, "PL")
                leg.AddEntry(g, f"{sp}, back+veto, E {'<' if side == 'lo' else '>'} cut", "pl")
        mg.SetTitle(f"photon peak #sigma vs edep cut [{var}];edep cut (MeV);#sigma (ns)")
        mg.Draw("A")
        leg.Draw()
        KEEP.append((mg, leg))
    save(c, "PH")

    # peak position vs run
    c = ROOT.TCanvas("c_run", "photon peak vs run", 1800, 700)
    c.Divide(2, 1)
    for i, var in enumerate(("tofc", "tofcrf")):
        pad = c.cd(i + 1)
        mg = ROOT.TMultiGraph()
        leg = ROOT.TLegend(0.12, 0.75, 0.4, 0.89)
        leg.SetBorderSize(0)
        for sp, col in (("P", "#2a78d6"), ("H", "#eb6834")):
            g = ROOT.TGraphErrors()
            for run, mu, emu, _ in runfits[(sp, var)]:
                n = g.GetN()
                g.SetPoint(n, run, mu)
                g.SetPointError(n, 0, emu)
            g.SetMarkerStyle(20)
            g.SetMarkerSize(0.8)
            g.SetMarkerColor(ROOT.TColor.GetColor(col))
            g.SetLineColor(ROOT.TColor.GetColor(col))
            mg.Add(g, "P")
            leg.AddEntry(g, f"{sp} ({SPEC_NAME[sp]})", "p")
        mg.SetTitle(f"photon peak position per run (back + front veto) [{var}];run;#mu (ns)")
        mg.Draw("A")
        leg.Draw()
        KEEP.append((mg, leg))
    save(c, "PH")

    # per-bar mu (tofc) map: plane x bar, P and H, all hits
    c = ROOT.TCanvas("c_bar_mu", "photon peak per bar", 1800, 700)
    c.Divide(2, 1)
    for i, sp in enumerate(SPECS):
        pad = c.cd(i + 1)
        pad.SetRightMargin(0.15)
        h = ROOT.TH2D(f"mu_{sp}", f"{sp}: photon #mu (tof-L/c, all hits) - plane total;bar;plane", 11, -0.5, 10.5, 5, -0.5, 4.5)
        for p in range(5):
            h.GetYaxis().SetBinLabel(p + 1, PLANES[p])
            ref = fits[(sp, "all", "tofc", "total", -1)]["mu"]
            for b in range(11):
                r = fits[(sp, "all", "tofc", PLANES[p], b)]
                if r["ok"]:
                    h.SetBinContent(b + 1, p + 1, r["mu"] - ref)
        h.SetMinimum(-1.5)
        h.SetMaximum(1.5)
        h.Draw("COLZ TEXT")
        h.SetMarkerSize(1.2)
        ROOT.gStyle.SetPaintTextFormat(".2f")
        KEEP.append(h)
    save(c, "PH")

    for i, cv in enumerate(canv_list):
        cv.Print(pdf + ("(" if i == 0 else (")" if i == len(canv_list) - 1 else "")), "pdf")
    fout.Close()

    # =================== summary ===================
    L = summary.append
    L(f"Input: {a.infile}")
    L(f"ToF offsets of the input: P (SHMS) = {off['P']} ns, H (HMS) = {off['H']} ns")
    L(f"Background: RF-periodic template ({RF_PERIOD} ns) from [mu0-{SB_LO}, mu0-{SB_HI}] ns; Gaussian core fit [mu-{FIT_NSIG_LO}, mu+{FIT_NSIG_HI}] sigma")
    L("")
    L("Photon peak, all-planes total (bars 1,9 of 100/101 excluded):")
    L(f"{'spec':4} {'var':7} {'category':28} {'mu':>8} {'+-':>6} {'sigma':>7} {'FWHM':>6} {'Nsig':>9} {'S/B':>6}")

    def line(sp, var, lab, r):
        L(f"{sp:4} {var:7} {lab:28} {r['mu']:8.3f} {r['emu']:6.3f} {r['sigma']:7.3f} {r['fwhm']:6.2f} {r['nsig']:9.0f} {r['sb']:6.2f}")

    for var in ("tofc", "tofcrf"):
        for sp in SPECS:
            for cat in ("all", "back", "back_fv"):
                line(sp, var, LABEL[cat], fits[(sp, cat, var, "total", -1)])
            for e in ecuts:
                line(sp, var, f"back+veto E<{e:g}", escan[(sp, "back_fv", var, "lo", e)])
                line(sp, var, f"back+veto E>{e:g}", escan[(sp, "back_fv", var, "hi", e)])
            for e in ecuts:
                line(sp, var, f"back (no veto) E<{e:g}", escan[(sp, "back", var, "lo", e)])
        L("")
    L(f"Energy split used in overlay plots: {eplot:g} MeV")
    L("")
    L("Per plane (all hits, tofc):")
    for sp in SPECS:
        L("  " + sp + ": " + "  ".join(f"{PLANES[p]}: {fits[(sp, 'all', 'tofc', PLANES[p], -1)]['mu']:.2f}/{fits[(sp, 'all', 'tofc', PLANES[p], -1)]['sigma']:.2f}" for p in range(5)) + "   (mu/sigma ns)")
    L("")
    L("Offset determination (photon tof - L/c should be 0):")
    for var in ("tofc", "tofcrf"):
        L(f"  [{var}]")
        mus = {}
        for sp in SPECS:
            r = fits[(sp, "back_fv", var, "total", -1)]
            ra = fits[(sp, "all", var, "total", -1)]
            mus[sp] = r["mu"]
            L(f"    {sp}: mu(back+veto) = {r['mu']:.3f} +- {r['emu']:.3f} ns, mu(all) = {ra['mu']:.3f} ns"
              f"  -> lglobal_time_offset_{'shms' if sp == 'P' else 'hms'} = {off[sp] - r['mu']:.2f} ns")
        L(f"    H - P = {mus['H'] - mus['P']:.3f} ns")
    L("")
    L("Per run period (back + front veto), periods = l_rf_offset rows:")
    for var in ("tofc", "tofcrf"):
        L(f"  [{var}]")
        for lo, hi in zip(PERIODS[:-1], PERIODS[1:]):
            rp, rh = perfits[("P", var, lo)], perfits[("H", var, lo)]
            if not (rp["ok"] or rh["ok"]):
                continue
            L(f"    runs {lo}-{hi - 1}: P mu={rp['mu']:.3f}+-{rp['emu']:.3f} sig={rp['sigma']:.3f} (N={rp['n']:.0f})  "
              f"H mu={rh['mu']:.3f}+-{rh['emu']:.3f} sig={rh['sigma']:.3f} (N={rh['n']:.0f})  H-P={rh['mu'] - rp['mu']:.3f}  "
              f"-> offset_P={off['P'] - rp['mu']:.2f} offset_H={off['H'] - rh['mu']:.2f}")
    L("")
    L("Run-period shift and the combined width (back + front veto, from the per-run histograms, 0.2 ns bins):")
    for var in ("tofc", "tofcrf"):
        for sp in SPECS:
            parts = []
            for tag, lab in (("all", "all runs"), ("excl", "without 22560-22660"), ("aligned", "periods aligned")):
                r = aligned[(sp, var, tag)]
                parts.append(f"{lab}: sig={r['sigma']:.3f} FWHM={r['fwhm']:.2f}" if r.get("ok") else f"{lab}: no fit")
            L(f"  [{var}] {sp}: " + " | ".join(parts))
    L("")
    L("Run dependence of the photon peak (back + front veto, tofc): runs with >=1500 hits")
    for sp in SPECS:
        pts = runfits[(sp, "tofc")]
        if pts:
            mus = [p[1] for p in pts]
            mean = sum(mus) / len(mus)
            rms = math.sqrt(sum((x - mean) ** 2 for x in mus) / len(mus))
            L(f"  {sp}: {len(pts)} runs, mean {mean:.3f}, rms {rms:.3f}, min {min(mus):.2f} (run {pts[mus.index(min(mus))][0]:.0f}), max {max(mus):.2f} (run {pts[mus.index(max(mus))][0]:.0f})")
    L("")
    L("Expected photon tof L/c (target centre, ypos=0), ns:")
    for p in range(5):
        L(f"  plane {PLANES[p]}: " + " ".join(f"{expected_Lc(p, b):.2f}" for b in range(11)))
    txt = "\n".join(summary)
    with open(os.path.join(a.outdir, "summary.txt"), "w") as fo:
        fo.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    sys.exit(main())
