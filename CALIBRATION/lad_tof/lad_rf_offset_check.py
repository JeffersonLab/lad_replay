#!/usr/bin/env python3
"""Check an RF offset table (l_rf_offset) against the LAD photon flash.

usage: python3 lad_rf_offset_check.py <lad_rf_offset_check.C output.root> <outdir>
           [--user rf_offset_results.root] [--table proposed.param] [--group-min 150000]

Input: lad_rf_offset_check.C, which stores per run and arm the RF phase
u = fmod(t_b - RF, T) of every event and (tof - L/c, u) of every photon candidate.
For an offset table the RF-corrected photon time is (tof - L/c) + remainder(u + rf_offset, T),
so any table is evaluated without reprocessing the data.

What is tested. Write t_vertex = t_bunch + d + jitter (d: spectrometer timing) and
RF = t_bunch - psi (psi: RF signal phase). The measured phase is phi = d + psi (mod T).
  - The photon peak WITHOUT the RF correction (mu_noRF) moves with -d and not with psi.
  - With a table that follows phi, the RF-corrected peak equals mu_noRF plus a constant
    (circular mean - mode of u); a table error e shifts it by e.
So phase changes that are RF (psi) leave mu_noRF flat, and the table must follow them;
phase changes that are spectrometer timing (d) move mu_noRF by -1 ns per ns of phase.
The drift test fits mu_noRF = a_period + b * (phi - phi_table) over run groups: b = 0
means the phase drift is in the RF, b = -1 means it is in the spectrometer time.

Per run: phase (get_RF_offset_fit.C method: Gaussian on the tiled phase histogram,
+-0.20/0.25/0.30 ns around the maximum, averaged; rf_offset = -phase). Photon peaks are
fitted with the width fixed to the summed sample (only the position is free), on the
accidental-background-subtracted histogram, per run and per run group (consecutive runs
in one table period with >= --group-min photon-candidate hits).

Outputs in <outdir>: runs.csv, groups.csv, summary.txt, plots.pdf (+ png).
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lad_photon_timing_fit import bgsub, fit_peak  # noqa: E402

ROOT.gErrorIgnoreLevel = ROOT.kFatal  # after the import (it sets its own level): rebin/empty-fit chatter

T = 4.00801
SPECS = ("P", "H")
COL = {"noRF": "#52514e", "cur": "#e34948", "run": "#2a78d6", "new": "#1baf7a"}
LAB = {"noRF": "no RF", "cur": "RF, current table", "run": "RF, per-run offset", "new": "RF, proposed table"}
DEFAULT_PARAM = os.path.join(HERE, "..", "..", "PARAM", "LAD", "LADKINE", "lladkine.param")


def rem(x):
    return np.remainder(np.asarray(x) + T / 2, T) - T / 2


def read_table(path):
    """l_rf_offset rows (run_start, SHMS, HMS) from a param file."""
    nums, inside = [], False
    for line in open(path):
        line = line.split(";")[0].strip()
        if not line:
            continue
        if "=" in line:
            if inside:
                break
            key, line = (x.strip() for x in line.split("=", 1))
            if key != "l_rf_offset":
                continue
            inside = True
        if inside:
            nums += [float(x) for x in line.replace(",", " ").split()]
    if not nums:
        raise ValueError(f"no l_rf_offset table in {path}")
    return [(int(nums[k]), nums[k + 1], nums[k + 2]) for k in range(0, len(nums) - 2, 3)]


def table_lookup(table, run, spec):
    v = table[0]
    for row in table:
        if run >= row[0]:
            v = row
    return v[1] if spec == "P" else v[2]


def phase_fit(h):
    """get_RF_offset_fit.C: maximum, Gaussian fits on the 3x tiled histogram, mean of 3 ranges."""
    n = h.GetNbinsX()
    ext = ROOT.TH1D(h.GetName() + "_ext", "", 3 * n, -T, 2 * T)
    for i in range(1, n + 1):
        for k in range(3):
            ext.SetBinContent(i + k * n, h.GetBinContent(i))
    cont = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    sm = np.convolve(np.concatenate([cont[-5:], cont, cont[:5]]), np.ones(11) / 11, "valid")
    xmax = h.GetBinCenter(int(np.argmax(sm)) + 1)
    means = []
    f = ROOT.TF1("pf", "gaus", -T, 2 * T)
    for w in (0.20, 0.25, 0.30):
        f.SetParameters(sm.max(), xmax, 0.3)
        if int(ext.Fit(f, "QNS0", "", xmax - w, xmax + w)) == 0:
            means.append(f.GetParameter(1) % T)
    if not means:
        return float("nan")
    m = np.array(means)
    return float((m[0] + rem(m - m[0]).mean()) % T)


def circ(h):
    n = h.GetNbinsX()
    c = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    x = np.array([h.GetBinCenter(i) for i in range(1, n + 1)])
    if c.sum() <= 0:
        return float("nan"), float("nan")
    a = 2 * np.pi * x / T
    C, S = (c * np.cos(a)).sum() / c.sum(), (c * np.sin(a)).sum() / c.sum()
    R = math.hypot(C, S)
    width = T / 2 / np.pi * math.sqrt(-2 * math.log(R)) if 0 < R < 1 else float("nan")
    return float(np.mod(math.atan2(S, C) * T / 2 / np.pi, T)), width


def th2_arrays(h2):
    nx, ny = h2.GetNbinsX(), h2.GetNbinsY()
    a = np.array([[h2.GetBinContent(i, j) for j in range(1, ny + 1)] for i in range(1, nx + 1)])
    xc = np.array([h2.GetXaxis().GetBinCenter(i) for i in range(1, nx + 1)])
    yc = np.array([h2.GetYaxis().GetBinCenter(j) for j in range(1, ny + 1)])
    return a, xc, yc, h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax()


def shifted_counts(arr, off):
    """(tof - L/c) + remainder(u + off) histogram contents; off=None -> no RF correction."""
    a, xc, yc, lo, hi = arr
    if off is None:
        return a.sum(axis=1)
    x = (xc[:, None] + rem(yc + off)[None, :]).ravel()
    return np.histogram(x, bins=len(xc), range=(lo, hi), weights=a.ravel())[0]


def to_hist(c, lo, hi, name):
    h = ROOT.TH1D(name, "", len(c), lo, hi)
    for i, v in enumerate(c):
        h.SetBinContent(i + 1, v)
        h.SetBinError(i + 1, math.sqrt(max(v, 0.)))
    h.SetEntries(float(np.sum(c)))
    return h


def fit_fixed(h, sigma, mu0=0.0, window=3.0):
    """Peak position with the width fixed: bg-subtracted, seeded at the local maximum."""
    hs, _ = bgsub(h, mu0)
    hc = hs.Clone(hs.GetName() + "_c")
    hc.Rebin(5)
    hc.GetXaxis().SetRangeUser(mu0 - window, mu0 + window)
    mu = hc.GetBinCenter(hc.GetMaximumBin())
    g = ROOT.TF1("gfix", "gaus", -50, 50)
    for _ in range(4):
        g.SetRange(mu - 1.5 * sigma, mu + 1.0 * sigma)
        g.SetParameters(max(hs.GetBinContent(hs.FindBin(mu)), 1.), mu, sigma)
        g.FixParameter(2, sigma)
        if int(hs.Fit(g, "QNRS0")) != 0:
            return float("nan"), float("nan")
        mu = g.GetParameter(1)
        if abs(mu - mu0) > window + 1:
            return float("nan"), float("nan")
    return mu, g.GetParError(1)


def read_user(path):
    """Per-run phases from get_RF_offset_fit.C output: {spec: {run: value}}."""
    out = {s: {} for s in SPECS}
    if not path or not os.path.exists(path):
        return out
    f = ROOT.TFile.Open(path)
    for s in SPECS:
        g = [o for o in f.Get("c_offsets_" + s).GetListOfPrimitives() if o.InheritsFrom("TGraph")][0]
        for i in range(g.GetN()):
            out[s][int(round(g.GetPointX(i)))] = g.GetPointY(i)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outdir")
    ap.add_argument("--user", default=None, help="get_RF_offset_fit.C output with per-run phases (for comparison)")
    ap.add_argument("--current", default=DEFAULT_PARAM, help="param file with the current l_rf_offset table")
    ap.add_argument("--table", default=None, help="param file with a proposed l_rf_offset table")
    ap.add_argument("--cat", default="bfv", choices=("bfv", "all"), help="photon sample")
    ap.add_argument("--min-events", type=int, default=500, help="min events with t_vertex for a run's phase")
    ap.add_argument("--group-min", type=float, default=150000, help="photon-candidate hits per run group")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fin = ROOT.TFile.Open(a.infile)
    user = read_user(a.user)
    tab_cur = read_table(a.current)
    tab_new = read_table(a.table) if a.table else None
    tabs = ("noRF", "cur", "run") + (("new",) if tab_new else ())
    runs = sorted(int(k.GetName()[3:]) for k in fin.GetListOfKeys() if k.GetName().startswith("run"))

    # ---- per run: phase and photon spectra for every table ----
    R = {s: {} for s in SPECS}
    for r in runs:
        d = fin.Get(f"run{r}")
        for s in SPECS:
            hu = d.Get(f"{s}_u_new")
            if not hu or hu.GetEntries() < a.min_events:
                continue
            ph = phase_fit(hu)
            cm, width = circ(hu)
            h2 = d.Get(f"{s}_{a.cat}")
            arr = th2_arrays(h2)
            offs = {"noRF": None, "cur": table_lookup(tab_cur, r, s), "run": -ph}
            if tab_new:
                offs["new"] = table_lookup(tab_new, r, s)
            R[s][r] = dict(run=r, spec=s, n_ev=int(hu.GetEntries()), n_ph=int(h2.GetEntries()), phase=ph,
                           phase_old=phase_fit(d.Get(f"{s}_u_old")), circ_minus_mode=float(rem(cm - ph)),
                           width=width, user_phase=user[s].get(r, float("nan")), table_cur=offs["cur"],
                           dphi=float(rem(ph + offs["cur"])), lo=arr[3], hi=arr[4],
                           counts={t: shifted_counts(arr, o) for t, o in offs.items()})
        print(f"\r[lad_rf_offset_check] run {r}", end="", file=sys.stderr)
    print(file=sys.stderr)

    # ---- summed peaks: reference widths ----
    lines = []
    L = lines.append
    L(f"input {a.infile}, photon sample '{a.cat}', {len(runs)} runs; current table {a.current}"
      + (f"; proposed {a.table}" if a.table else ""))
    L("")
    L("Summed photon peak (all runs), free Gaussian core fit:")
    sig_ref = {}
    for s in SPECS:
        rr = list(R[s].values())
        if not rr:
            continue
        for t in tabs:
            h = to_hist(sum(x["counts"][t] for x in rr), rr[0]["lo"], rr[0]["hi"], f"sum_{s}_{t}")
            fr, _ = fit_peak(h, 0.0, search=4.0)
            sig_ref[(s, t)] = fr["sigma"]
            L(f"  {s} {t:5s}: mu {fr['mu']:7.3f} +- {fr['emu']:.3f}  sigma {fr['sigma']:.3f}  FWHM {fr['fwhm']:.2f}  Nsig {fr['nsig']:.0f}")
    L("")

    # ---- per run and per group fixed-width fits ----
    def fit_counts(s, c, lo, hi, t, name):
        return fit_fixed(to_hist(c, lo, hi, name), sig_ref[(s, t)])

    for s in SPECS:
        for x in R[s].values():
            for t in tabs:
                x[f"mu_{t}"], x[f"emu_{t}"] = fit_counts(s, x["counts"][t], x["lo"], x["hi"], t, f"r_{s}_{x['run']}_{t}")
    bounds = [row[0] for row in tab_cur] + [10 ** 9]
    groups = []
    for s in SPECS:
        keys = sorted(R[s])
        for b0, b1 in zip(bounds[:-1], bounds[1:]):
            rs = [r for r in keys if b0 <= r < b1]
            cur, acc, gl = [], 0, []
            for r in rs:
                cur.append(r)
                acc += R[s][r]["n_ph"]
                if acc >= a.group_min:
                    gl.append(cur)
                    cur, acc = [], 0
            if cur:
                if gl and acc < a.group_min / 2:
                    gl[-1] += cur
                else:
                    gl.append(cur)
            for g in gl:
                xs = [R[s][r] for r in g]
                w = np.array([x["n_ev"] for x in xs], float)
                ph = np.array([x["phase"] for x in xs])
                grp = dict(spec=s, first=g[0], last=g[-1], nrun=len(g), n_ph=int(sum(x["n_ph"] for x in xs)),
                           phase=float((ph[0] + np.sum(w * rem(ph - ph[0])) / w.sum()) % T),
                           table_cur=xs[0]["table_cur"])
                grp["dphi"] = float(rem(grp["phase"] + grp["table_cur"]))
                for t in tabs:
                    grp[f"mu_{t}"], grp[f"emu_{t}"] = fit_counts(s, sum(x["counts"][t] for x in xs), xs[0]["lo"],
                                                                 xs[0]["hi"], t, f"g_{s}_{g[0]}_{t}")
                groups.append(grp)

    rkeys = ["run", "spec", "n_ev", "n_ph", "phase", "phase_old", "user_phase", "circ_minus_mode", "width", "table_cur",
             "dphi"] + [f"{p}_{t}" for t in tabs for p in ("mu", "emu")]
    with open(os.path.join(a.outdir, "runs.csv"), "w", newline="") as fo:
        wr = csv.DictWriter(fo, fieldnames=rkeys, extrasaction="ignore")
        wr.writeheader()
        for s in SPECS:
            for x in R[s].values():
                wr.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in x.items() if k in rkeys})
    gkeys = ["spec", "first", "last", "nrun", "n_ph", "phase", "table_cur", "dphi"] + [f"{p}_{t}" for t in tabs for p in ("mu", "emu")]
    with open(os.path.join(a.outdir, "groups.csv"), "w", newline="") as fo:
        wr = csv.DictWriter(fo, fieldnames=gkeys)
        wr.writeheader()
        for g in groups:
            wr.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in g.items()})

    # ---- consistency with get_RF_offset_fit and the drift test ----
    for s in SPECS:
        xs = list(R[s].values())
        dd = np.array([rem(x["phase"] - x["user_phase"]) for x in xs])
        bad = [(x["run"], round(float(v), 2)) for x, v in zip(xs, dd) if np.isfinite(v) and abs(v) > 0.05]
        L(f"{s}: phase - get_RF_offset_fit phase over {np.isfinite(dd).sum()} runs: median {np.nanmedian(dd) if np.isfinite(dd).any() else float('nan'):+.3f} ns; |diff| > 0.05 ns: {bad}")
        L(f"   phase width (circular, Gaussian-equivalent) {min(x['width'] for x in xs):.2f}-{max(x['width'] for x in xs):.2f} ns; "
          f"circular mean - mode {np.median([x['circ_minus_mode'] for x in xs]):+.3f} ns; z-corrected - old phase max {max(abs(rem(x['phase'] - x['phase_old'])) for x in xs):.3f} ns")
        gs = [g for g in groups if g["spec"] == s and np.isfinite(g["mu_noRF"]) and g["emu_noRF"] < 0.3]
        per = sorted({next(i for i, b in enumerate(bounds) if b > g["first"]) for g in gs})
        if len(gs) > len(per) + 1:
            X = np.array([[1.0 if next(i for i, b in enumerate(bounds) if b > g["first"]) == p else 0.0 for p in per] + [g["dphi"]] for g in gs])
            y = np.array([g["mu_noRF"] for g in gs])
            W = 1 / np.array([g["emu_noRF"] for g in gs]) ** 2
            A = X.T @ (W[:, None] * X)
            beta = np.linalg.solve(A, X.T @ (W * y))
            cov = np.linalg.inv(A)
            chi2 = float(np.sum(W * (y - X @ beta) ** 2))
            L(f"   drift test: d(mu_noRF)/d(phase - table) = {beta[-1]:+.2f} +- {math.sqrt(cov[-1, -1]):.2f} "
              f"(RF drift -> 0, spectrometer drift -> -1), phase range {min(g['dphi'] for g in gs):+.3f}..{max(g['dphi'] for g in gs):+.3f} ns, chi2/ndf {chi2:.1f}/{len(gs) - len(beta)}")
            for p, b, e in zip(per, beta, np.sqrt(np.diag(cov))):
                L(f"     table period from run {bounds[p - 1]}: mu_noRF {b:+.3f} +- {e:.3f} ns")
    L("")
    L("Run groups (fixed-width fits; mu in ns):")
    L(f"{'sp':2} {'runs':>11} {'nrun':>4} {'n_ph':>8} {'phase':>6} {'-table':>6} | " + "  ".join(f"{t:>13}" for t in tabs))
    for g in groups:
        L(f"{g['spec']:2} {g['first']:5d}-{g['last']:5d} {g['nrun']:4d} {g['n_ph']:8d} {g['phase']:6.3f} {(-g['table_cur']) % T:6.3f} | "
          + "  ".join(f"{g[f'mu_{t}']:6.3f}+-{g[f'emu_{t}']:.3f}" for t in tabs))
    with open(os.path.join(a.outdir, "summary.txt"), "w") as fo:
        fo.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # ---- plots ----
    pdf = os.path.join(a.outdir, "plots.pdf")
    cv, keep = [], []
    for s in SPECS:
        xs = list(R[s].values())
        if not xs:
            continue
        rlo, rhi = min(x["run"] for x in xs) - 20, max(x["run"] for x in xs) + 20
        c = ROOT.TCanvas(f"c_{s}", s, 1600, 1200)
        c.Divide(1, 3)
        # 1: phase per run, with the tables as steps
        c.cd(1)
        mg, leg = ROOT.TMultiGraph(), ROOT.TLegend(0.12, 0.62, 0.42, 0.88)
        for lab, key, mk, col in (("phase (this check)", "phase", 20, "#2a78d6"),
                                  ("phase (get_RF_offset_fit)", "user_phase", 24, "#000000")):
            g = ROOT.TGraph()
            for x in xs:
                if np.isfinite(x[key]):
                    g.SetPoint(g.GetN(), x["run"], x[key] % T)
            if g.GetN():
                g.SetMarkerStyle(mk)
                g.SetMarkerColor(ROOT.TColor.GetColor(col))
                mg.Add(g, "P")
                leg.AddEntry(g, lab, "p")
        for lab, tab, col in (("-table (current)", tab_cur, COL["cur"]), ("-table (proposed)", tab_new, COL["new"])):
            if not tab:
                continue
            g = ROOT.TGraph()
            rows = [row for row in tab if row[0] <= rhi]
            for i, row in enumerate(rows):
                x0 = max(row[0], rlo)
                x1 = rows[i + 1][0] if i + 1 < len(rows) else rhi
                if x1 <= rlo:
                    continue
                v = (-(row[1] if s == "P" else row[2])) % T
                g.SetPoint(g.GetN(), x0, v)
                g.SetPoint(g.GetN(), x1, v)
            g.SetLineColor(ROOT.TColor.GetColor(col))
            g.SetLineWidth(2)
            mg.Add(g, "L")
            leg.AddEntry(g, lab, "l")
        mg.SetTitle(f"{s}: RF phase per run;run;fmod(t_{{b}} - RF, T) (ns)")
        mg.Draw("A")
        mg.GetXaxis().SetLimits(rlo, rhi)
        leg.Draw()
        keep += [mg, leg]
        # 2: photon peak per run group
        c.cd(2)
        mg, leg = ROOT.TMultiGraph(), ROOT.TLegend(0.12, 0.62, 0.42, 0.88)
        for k, t in enumerate(tabs):
            g = ROOT.TGraphErrors()
            for gr in groups:
                if gr["spec"] == s and np.isfinite(gr[f"mu_{t}"]) and gr[f"emu_{t}"] < 0.5:
                    n = g.GetN()
                    g.SetPoint(n, 0.5 * (gr["first"] + gr["last"]) + 3 * (k - 1), gr[f"mu_{t}"])
                    g.SetPointError(n, 0.5 * (gr["last"] - gr["first"]), gr[f"emu_{t}"])
            g.SetMarkerStyle(20 + k)
            g.SetMarkerColor(ROOT.TColor.GetColor(COL[t]))
            g.SetLineColor(ROOT.TColor.GetColor(COL[t]))
            mg.Add(g, "P")
            leg.AddEntry(g, LAB[t], "p")
        mg.SetTitle(f"{s}: photon peak per run group (width fixed);run;#mu (tof - L/c) (ns)")
        mg.Draw("A")
        mg.GetXaxis().SetLimits(rlo, rhi)
        leg.Draw()
        keep += [mg, leg]
        # 3: drift test
        c.cd(3)
        g = ROOT.TGraphErrors()
        for gr in groups:
            if gr["spec"] == s and np.isfinite(gr["mu_noRF"]) and gr["emu_noRF"] < 0.3:
                n = g.GetN()
                g.SetPoint(n, gr["dphi"], gr["mu_noRF"])
                g.SetPointError(n, 0, gr["emu_noRF"])
        g.SetMarkerStyle(20)
        g.SetTitle(f"{s}: photon peak without RF vs RF phase - table (run groups);phase - (-table) (ns);#mu_{{no RF}} (ns)")
        g.Draw("AP")
        keep.append(g)
        cv.append(c)
    for i, c in enumerate(cv):
        c.Print(pdf + ("(" if i == 0 and len(cv) > 1 else ")" if i == len(cv) - 1 and len(cv) > 1 else ""), "pdf")
        c.SaveAs(os.path.join(a.outdir, c.GetName() + ".png"))


if __name__ == "__main__":
    main()
