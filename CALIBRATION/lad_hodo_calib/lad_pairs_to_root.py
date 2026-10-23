#!/usr/bin/env python3
"""Sum lad_pair_histos.py outputs and write their h_dtb/h_dtf_{int,amp} histograms into a ROOT file, in a period
directory, in the layout lad_edep_back.py and lad_edep_transfer.py read.

usage: lad_pairs_to_root.py out.root period_dir in1.npz [in2.npz ...]
"""
import sys

import numpy as np

import ROOT


def main():
    out, per, ins = sys.argv[1], sys.argv[2], sys.argv[3:]
    z0 = np.load(ins[0])
    keys = [k for k in z0.files if k.startswith("h_dt")]
    S = {k: 0 for k in keys}
    for fn in ins:
        z = np.load(fn)
        for k in keys:
            S[k] = S[k] + z[k]
    f = ROOT.TFile(out, "RECREATE")
    d = f.mkdir(per)
    d.cd()
    for k in keys:
        ye = z0["IB"] if "_int_" in k else z0["AB"]
        h = ROOT.TH2F(k, k, len(z0["DTB"]) - 1, z0["DTB"], len(ye) - 1, ye)
        for i in range(S[k].shape[0]):
            for j in range(S[k].shape[1]):
                if S[k][i, j]:
                    h.SetBinContent(i + 1, j + 1, S[k][i, j])
        h.SetEntries(S[k].sum())
        h.Write()
    f.Close()
    print(f"wrote {out}: {len(keys)} histograms from {len(ins)} files")


if __name__ == "__main__":
    main()
