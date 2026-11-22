// lad_tw_histos.C
//
// Time-walk histograms for the LAD hodoscope, one per PMT:
//   (TDC time - FADC pulse time) vs FADC pulse amplitude
// from the per-plane Good* arrays (X.ladhod.<plane>.Good{Top,Btm}AdcTdcDiffTime and ...AdcPulseAmp).
// The FADC pulse time is a constant-fraction time, so the amplitude dependence of the difference
// is the TDC (leading-edge discriminator) walk. Fit the output with lad_tw_fit.py.
//
// The TDC-ADC window is fixed (DT_LO..DT_HI, the ADC-TDC offsets put every PMT near 0) so that the
// outputs of different chunks can be hadd-ed.
//
// Usage:
//   root -l -b -q 'lad_tw_histos.C+("files.dat","out.root",4,"P")'
//   files.dat: one ROOT file per line ('#' comments allowed). spec: "P" or "H".

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TROOT.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace ladtw {
const int NPLANES                = 6;
const char *PLANES[NPLANES]      = {"000", "001", "100", "101", "200", "REFBAR"};
const int NPAD[NPLANES]          = {11, 11, 11, 11, 11, 1};
const char *SIDES[2]             = {"Top", "Btm"};
const double DT_LO = -15., DT_HI = 35.; // TDC - FADC time window (ns); walk makes the TDC late at low amplitude
const int NBIN_DT  = 1000;
const int NBIN_AMP               = 250;
const double AMP_MAX             = 1000.; // mV

typedef ROOT::RVec<double> RV;

int pmt_index(int ipl, int is, int ipad) { return (ipl * 2 + is) * 11 + ipad; }
} // namespace ladtw

void lad_tw_histos(const char *listfile, const char *outfile, int nthreads = 4, const char *spec = "P") {
  using namespace ladtw;
  TChain ch("T");
  std::ifstream in(listfile);
  std::string line;
  int nfiles = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    ch.Add(line.c_str());
    nfiles++;
  }
  std::cout << "files: " << nfiles << std::endl;
  if (nthreads > 1)
    ROOT::EnableImplicitMT(nthreads);
  ROOT::RDataFrame df(ch);
  const unsigned nslots = df.GetNSlots();

  // Column names, in pmt_index order: for each plane, side: (diff, amp)
  std::vector<std::string> cols;
  for (int ipl = 0; ipl < NPLANES; ipl++)
    for (int is = 0; is < 2; is++) {
      cols.push_back(Form("%s.ladhod.%s.Good%sAdcTdcDiffTime", spec, PLANES[ipl], SIDES[is]));
      cols.push_back(Form("%s.ladhod.%s.Good%sAdcPulseAmp", spec, PLANES[ipl], SIDES[is]));
    }

  // One function handles all 24 column pairs: pack them with a lambda taking a fixed argument list.
  auto loop = [&](auto &&fill) {
    df.ForeachSlot(
        [&](unsigned slot, const RV &d0, const RV &a0, const RV &d1, const RV &a1, const RV &d2, const RV &a2,
            const RV &d3, const RV &a3, const RV &d4, const RV &a4, const RV &d5, const RV &a5, const RV &d6,
            const RV &a6, const RV &d7, const RV &a7, const RV &d8, const RV &a8, const RV &d9, const RV &a9,
            const RV &d10, const RV &a10, const RV &d11, const RV &a11) {
          const RV *D[12] = {&d0, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11};
          const RV *A[12] = {&a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7, &a8, &a9, &a10, &a11};
          for (int k = 0; k < 12; k++) {
            int ipl = k / 2, is = k % 2;
            const RV &d = *D[k];
            const RV &a = *A[k];
            for (size_t ip = 0; ip < d.size() && ip < a.size() && (int)ip < NPAD[ipl]; ip++) {
              if (!(d[ip] > -1e5 && d[ip] < 1e5) || !(a[ip] > 0))
                continue;
              fill(slot, pmt_index(ipl, is, ip), d[ip], a[ip]);
            }
          }
        },
        cols);
  };

  const int NPMT = NPLANES * 2 * 11;

  // ---- 2D histograms, per slot then merged
  std::vector<std::vector<TH2F *>> h(nslots, std::vector<TH2F *>(NPMT, nullptr));
  for (unsigned s = 0; s < nslots; s++)
    for (int ipl = 0; ipl < NPLANES; ipl++)
      for (int is = 0; is < 2; is++)
        for (int ip = 0; ip < NPAD[ipl]; ip++) {
          int i  = pmt_index(ipl, is, ip);
          h[s][i] = new TH2F(Form("h_tw_%s_%s_%d_s%u", PLANES[ipl], SIDES[is], ip + 1, s),
                             Form("Plane %s %s paddle %d;FADC pulse amplitude (mV);TDC - FADC time (ns)",
                                  PLANES[ipl], SIDES[is], ip + 1),
                             NBIN_AMP, 0, AMP_MAX, NBIN_DT, DT_LO, DT_HI);
          h[s][i]->SetDirectory(nullptr);
        }
  loop([&](unsigned slot, int i, double d, double a) {
    if (h[slot][i])
      h[slot][i]->Fill(a, d);
  });

  TFile fout(outfile, "RECREATE");
  for (int ipl = 0; ipl < NPLANES; ipl++) {
    TDirectory *dir = fout.mkdir(PLANES[ipl]);
    dir->cd();
    for (int is = 0; is < 2; is++)
      for (int ip = 0; ip < NPAD[ipl]; ip++) {
        int i     = pmt_index(ipl, is, ip);
        TH2F *sum = (TH2F *)h[0][i]->Clone(Form("h_tw_%s_%s_%d", PLANES[ipl], SIDES[is], ip + 1));
        for (unsigned s = 1; s < nslots; s++)
          sum->Add(h[s][i]);
        sum->Write();
      }
  }
  fout.cd();
  TNamed("spec", spec).Write();
  TNamed("nfiles", Form("%d", nfiles)).Write();
  fout.Close();
  std::cout << "wrote " << outfile << std::endl;
}
