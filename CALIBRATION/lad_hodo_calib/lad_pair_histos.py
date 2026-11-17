#!/usr/bin/env python3
"""Same-paddle front-back pairs (000/001, 100/101) recomputed from the per-PMT Good arrays with two LAD parameter
sets (set 0 "old", set 1 "new"; params file from lad_bar_timing_params.py), with no pairing cut: back - front y and
time, the back-bar pulse integral against dt for |dy| < 20 cm, and the raw per-PMT occupancy of plane 101. For set 1
also h_{dtb,dtf}_{int,amp}_<back plane>_<paddle>: dt against the back/front bar sqrt(Q_top Q_btm) for |dy| < 20 cm
and |dt| < 10 ns, in the lad_edep_histos.C binning (lad_pairs_to_root.py writes them for lad_edep_back.py).
Unlike LADlib's good-hit pairing (same paddle, |dt| < 10 ns, |dy| < 20 cm with the replay's constants) this does not
depend on the constants the file was replayed with. Input: a replay with the Good* per-PMT arrays (production DEF).

usage: lad_pair_histos.py replay.root params.txt out.npz      (sum the outputs with numpy; lad_bar_from_pairs.py)
"""
import uproot, numpy as np, awkward as ak, sys
fn, parf, out = sys.argv[1], sys.argv[2], sys.argv[3]
PL = ["000", "001", "100", "101"]
par = {}
thr = {}
for l in open(parf):
    w = l.split()
    if w[1] == "thr":
        thr[int(w[0])] = float(w[2]); continue
    s, pl, side, pad = int(w[0]), w[1], w[2], int(w[3])
    par[(s, pl, side, pad)] = tuple(float(x) for x in w[4:9])  # c2 c3 cab lco vel
t = uproot.open(fn)["T"]
br = []
for pl in PL:
    for x in ["GoodTopTdcTimeUnCorr", "GoodBtmTdcTimeUnCorr", "GoodTopAdcPulseAmp", "GoodBtmAdcPulseAmp", "GoodTopAdcPulseInt", "GoodBtmAdcPulseInt"]:
        br.append(f"P.ladhod.{pl}.{x}")
for x in ["TopTdcCounter", "BtmTdcCounter", "TopAdcCounter", "BtmAdcCounter"]:
    br.append(f"P.ladhod.101.{x}")
a = t.arrays(br, library="np")
n = len(a[br[0]])
def arr(pl, x):
    v = a[f"P.ladhod.{pl}.{x}"]
    return np.stack([np.pad(np.asarray(e, float), (0, 11 - len(e)), constant_values=1e30)[:11] for e in v])
T = {pl: {k: arr(pl, k) for k in ["GoodTopTdcTimeUnCorr", "GoodBtmTdcTimeUnCorr", "GoodTopAdcPulseAmp", "GoodBtmAdcPulseAmp", "GoodTopAdcPulseInt", "GoodBtmAdcPulseInt"]} for pl in PL}
def hit(pl, s):
    tt, tb = T[pl]["GoodTopTdcTimeUnCorr"], T[pl]["GoodBtmTdcTimeUnCorr"]
    at, ab = T[pl]["GoodTopAdcPulseAmp"], T[pl]["GoodBtmAdcPulseAmp"]
    ok = (np.abs(tt) < 1e8) & (np.abs(tb) < 1e8) & (at > 0) & (ab > 0) & (at < 1e8) & (ab < 1e8)
    th, y = np.full(tt.shape, np.nan), np.full(tt.shape, np.nan)
    for b in range(11):
        c2t, c3t, cab, lco, vel = par[(s, pl, "Top", b + 1)]
        c2b, c3b = par[(s, pl, "Btm", b + 1)][:2]
        tw = lambda A, c2, c3: c3 * ((A / thr[s]) ** -c2 - (200 / thr[s]) ** -c2)
        with np.errstate(all="ignore"):
            ct = tt[:, b] - tw(at[:, b], c2t, c3t) + lco
            cb = tb[:, b] - tw(ab[:, b], c2b, c3b) - 2 * cab + lco
        th[:, b] = np.where(ok[:, b], 0.5 * (ct + cb), np.nan)
        y[:, b] = np.where(ok[:, b], 0.5 * (cb - ct) * vel, np.nan)
    return th, y, ok
res = {}
ybins, tbins, ebins = np.arange(-300, 301, 2.0), np.arange(-10, 10.01, 0.1), np.arange(0, 250.1, 1)
DTB, IB, AB = np.linspace(-4, 12, 161), np.linspace(0, 250, 251), np.linspace(0, 1000, 251)
for s, tag in [(0, "old"), (1, "new")]:
    H = {pl: hit(pl, s) for pl in PL}
    for f_, b_ in [("000", "001"), ("100", "101")]:
        tf, yf, okf = H[f_]; tb_, yb, okb = H[b_]
        for b in range(11):
            m = okf[:, b] & okb[:, b] & (np.abs(tb_[:, b] - tf[:, b]) < 10)
            dy = yb[m, b] - yf[m, b]; dt = tb_[m, b] - tf[m, b]
            res[f"{tag}_dy_{b_}_{b+1}"] = np.histogram(dy, ybins)[0]
            res[f"{tag}_dt_{b_}_{b+1}"] = np.histogram(dt, tbins)[0]
            m2 = m.copy(); m2[m] = np.abs(dy) < 20
            ib = np.sqrt(T[b_]["GoodTopAdcPulseInt"][m2, b] * T[b_]["GoodBtmAdcPulseInt"][m2, b])
            res[f"{tag}_dte_{b_}_{b+1}"] = np.histogram2d((tb_ - tf)[m2, b], ib, [tbins, ebins])[0]
            if tag == "new":  # lad_edep_histos.C binning: dt -4..12 ns (160), integral 0-250 pC (250), amplitude 0-1000 mV (250)
                d2 = (tb_ - tf)[m2, b]
                for pl_, nm in [(b_, "dtb"), (f_, "dtf")]:
                    qi = np.sqrt(T[pl_]["GoodTopAdcPulseInt"][m2, b] * T[pl_]["GoodBtmAdcPulseInt"][m2, b])
                    qa = np.sqrt(T[pl_]["GoodTopAdcPulseAmp"][m2, b] * T[pl_]["GoodBtmAdcPulseAmp"][m2, b])
                    res[f"h_{nm}_int_{b_}_{b+1}"] = np.histogram2d(d2, qi, [DTB, IB])[0]
                    res[f"h_{nm}_amp_{b_}_{b+1}"] = np.histogram2d(d2, qa, [DTB, AB])[0]
            res[f"{tag}_yb_{b_}_{b+1}"] = np.histogram(yb[okb[:, b], b], ybins)[0]
for x in ["TopTdcCounter", "BtmTdcCounter", "TopAdcCounter", "BtmAdcCounter"]:
    v = np.concatenate([np.asarray(e) for e in a[f"P.ladhod.101.{x}"]])
    res[f"occ_101_{x}"] = np.histogram(v, np.arange(0.5, 12.5, 1))[0]
res["nev"] = np.array([n])
np.savez(out, ybins=ybins, tbins=tbins, ebins=ebins, DTB=DTB, IB=IB, AB=AB, **res)
print("wrote", out, n)
