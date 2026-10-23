// lad_laser_histos.C
//
// LAD laser runs: per-PMT histograms of (TDC time - photodiode TDC time) vs FADC pulse amplitude.
// Laser runs have no spectrometer trigger, so the LAD TDCs have no reference time and hcana's
// Good* arrays are empty; this macro works from the raw per-plane arrays and uses the photodiode
// (same V1190 as the LAD TDCs) as the laser time reference.
//
// For each PMT the TDC hit closest to the laser peak and the largest FADC pulse are used.
// Pass 1 finds the per-PMT laser peak on the first NPEAK events, pass 2 fills the histograms.
// The output feeds lad_tw_fit.py (time walk: combine the attenuation scan) and lad_laser_offsets.py
// (bar-to-bar and top-bottom offsets).
//
// Usage:
//   root -l -b -q 'lad_laser_histos.C+("files.dat","out.root",1,"P")'   (nthreads is ignored)

#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <memory>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace ladlas {
const int NPLANES           = 6;
const char *PLANES[NPLANES] = {"000", "001", "100", "101", "200", "REFBAR"};
const int NPAD[NPLANES]     = {11, 11, 11, 11, 11, 1};
const char *SIDES[2]        = {"Top", "Btm"};
const double TDC2NS         = 0.09766; // V1190, ns per channel
const Long64_t NPEAK        = 20000;
const double WIN_LO = -8., WIN_HI = 22.;
const int NBIN_DT  = 600;
const int NBIN_AMP = 250;
const double AMP_MAX = 1000.;

typedef TTreeReaderArray<double> RV;
int pmt_index(int ipl, int is, int ipad) { return (ipl * 2 + is) * 11 + ipad; }

// For one plane/side: best (TDC - PD, amp, adc time) per paddle. TDC hit closest to target[paddle].
struct Best {
  double dt[11], amp[11], adct[11];
};
inline void pick(const ROOT::RVec<double> &tdc, const ROOT::RVec<double> &tcnt, const ROOT::RVec<double> &amp,
                 const ROOT::RVec<double> &acnt, const ROOT::RVec<double> &atime, double pd,
                 const double *target, int npad, Best &b) {
  for (int i = 0; i < 11; i++) {
    b.dt[i]   = NAN;
    b.amp[i]  = -1;
    b.adct[i] = NAN;
  }
  for (size_t k = 0; k < tdc.size() && k < tcnt.size(); k++) {
    int ip = (int)tcnt[k] - 1;
    if (ip < 0 || ip >= npad)
      continue;
    double dt = tdc[k] * TDC2NS - pd;
    if (std::isnan(b.dt[ip]) || std::fabs(dt - target[ip]) < std::fabs(b.dt[ip] - target[ip]))
      b.dt[ip] = dt;
  }
  for (size_t k = 0; k < amp.size() && k < acnt.size(); k++) {
    int ip = (int)acnt[k] - 1;
    if (ip < 0 || ip >= npad)
      continue;
    if (amp[k] > b.amp[ip]) {
      b.amp[ip]  = amp[k];
      b.adct[ip] = k < atime.size() ? atime[k] : NAN;
    }
  }
}
} // namespace ladlas

void lad_laser_histos(const char *listfile, const char *outfile, int /*nthreads*/ = 1, const char *spec = "P") {
  using namespace ladlas;
  TChain ch("T");
  {
    std::ifstream in(listfile);
    std::string line;
    int n = 0;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#') {
        ch.Add(line.c_str());
        n++;
      }
    std::cout << "files: " << n << ", entries: " << ch.GetEntries() << std::endl;
  }
  const int NPMT = NPLANES * 2 * 11;
  const std::string pdcol = Form("T.%s.photodiodeLAD_tdcTimeRaw", std::string(spec) == "P" ? "shms" : "hms");

  TTreeReader rd(&ch);
  TTreeReaderValue<double> pdv(rd, pdcol.c_str());
  // per plane/side: TdcTimeRaw, TdcCounter, AdcPulseAmp, AdcCounter, AdcPulseTime
  std::vector<std::unique_ptr<RV>> arr;
  for (int ipl = 0; ipl < NPLANES; ipl++)
    for (int is = 0; is < 2; is++)
      for (const char *v : {"TdcTimeRaw", "TdcCounter", "AdcPulseAmp", "AdcCounter", "AdcPulseTime"})
        arr.emplace_back(new RV(rd, Form("%s.ladhod.%s.%s%s", spec, PLANES[ipl], SIDES[is], v)));

  // ---- pass 1: peak of (TDC - PD) per PMT
  std::vector<double> peak(NPMT, 0.);
  {
    std::vector<TH1D *> h1(NPMT);
    for (int i = 0; i < NPMT; i++)
      h1[i] = new TH1D(Form("pk_%d", i), "", 20000, -1000, 1000);
    rd.SetEntriesRange(0, std::min<Long64_t>(NPEAK, ch.GetEntries()));
    while (rd.Next()) {
      if (!(*pdv > 0))
        continue;
      double pdt = *pdv * TDC2NS;
      for (int ipl = 0; ipl < NPLANES; ipl++)
        for (int is = 0; is < 2; is++) {
          int k = (ipl * 2 + is) * 5;
          RV &tdc = *arr[k], &tcnt = *arr[k + 1];
          for (size_t j = 0; j < tdc.GetSize() && j < tcnt.GetSize(); j++) {
            int ip = (int)tcnt[j] - 1;
            if (ip >= 0 && ip < NPAD[ipl])
              h1[pmt_index(ipl, is, ip)]->Fill(tdc[j] * TDC2NS - pdt);
          }
        }
    }
    for (int i = 0; i < NPMT; i++) {
      peak[i] = h1[i]->GetEntries() > 20 ? h1[i]->GetBinCenter(h1[i]->GetMaximumBin()) : 0.;
      delete h1[i];
    }
  }

  // ---- pass 2
  std::vector<TH2F *> h(NPMT, nullptr), hadc(NPMT, nullptr);
  std::vector<TH1D *> htdc(NPMT, nullptr); // TDC - PD without an FADC requirement (runs before 22590 have no laser FADC pulses)
  TH1D *hpd = new TH1D("h_pd", "photodiode raw TDC;ns", 4000, 0, 4000);
  hpd->SetDirectory(nullptr);
  for (int ipl = 0; ipl < NPLANES; ipl++)
    for (int is = 0; is < 2; is++)
      for (int ip = 0; ip < NPAD[ipl]; ip++) {
        int i = pmt_index(ipl, is, ip);
        h[i]  = new TH2F(Form("h_las_%s_%s_%d", PLANES[ipl], SIDES[is], ip + 1),
                        Form("Plane %s %s paddle %d;FADC pulse amplitude (mV);TDC - photodiode (ns)", PLANES[ipl],
                             SIDES[is], ip + 1),
                        NBIN_AMP, 0, AMP_MAX, NBIN_DT, peak[i] + WIN_LO, peak[i] + WIN_HI);
        hadc[i] = new TH2F(Form("h_adct_%s_%s_%d", PLANES[ipl], SIDES[is], ip + 1),
                           Form("Plane %s %s paddle %d;FADC pulse amplitude (mV);TDC - FADC time (ns)", PLANES[ipl],
                                SIDES[is], ip + 1),
                           NBIN_AMP, 0, AMP_MAX, 2000, -1000, 1000);
        h[i]->SetDirectory(nullptr);
        hadc[i]->SetDirectory(nullptr);
        htdc[i] = new TH1D(Form("h_tdc_%s_%s_%d", PLANES[ipl], SIDES[is], ip + 1),
                           Form("Plane %s %s paddle %d;TDC - photodiode (ns)", PLANES[ipl], SIDES[is], ip + 1), NBIN_DT,
                           peak[i] + WIN_LO, peak[i] + WIN_HI);
        htdc[i]->SetDirectory(nullptr);
      }
  rd.SetEntriesRange(0, ch.GetEntries());
  rd.Restart();
  Best b;
  std::vector<double> tdc, tcnt, amp, acnt, atime;
  while (rd.Next()) {
    if (!(*pdv > 0))
      continue;
    double pdt = *pdv * TDC2NS;
    hpd->Fill(pdt);
    for (int ipl = 0; ipl < NPLANES; ipl++)
      for (int is = 0; is < 2; is++) {
        int k = (ipl * 2 + is) * 5;
        ROOT::RVec<double> v[5];
        for (int m = 0; m < 5; m++)
          v[m] = ROOT::RVec<double>(arr[k + m]->begin(), arr[k + m]->end());
        pick(v[0], v[1], v[2], v[3], v[4], pdt, &peak[pmt_index(ipl, is, 0)], NPAD[ipl], b);
        for (int ip = 0; ip < NPAD[ipl]; ip++) {
          if (!std::isnan(b.dt[ip]))
            htdc[pmt_index(ipl, is, ip)]->Fill(b.dt[ip]);
          if (std::isnan(b.dt[ip]) || !(b.amp[ip] > 0))
            continue;
          int i = pmt_index(ipl, is, ip);
          h[i]->Fill(b.amp[ip], b.dt[ip]);
          if (!std::isnan(b.adct[ip]))
            hadc[i]->Fill(b.amp[ip], b.dt[ip] + pdt - b.adct[ip]);
        }
      }
  }

  TFile fout(outfile, "RECREATE");
  for (int ipl = 0; ipl < NPLANES; ipl++) {
    fout.mkdir(PLANES[ipl])->cd();
    for (int is = 0; is < 2; is++)
      for (int ip = 0; ip < NPAD[ipl]; ip++) {
        int i = pmt_index(ipl, is, ip);
        h[i]->Write();
        hadc[i]->Write();
        htdc[i]->Write();
      }
  }
  fout.cd();
  hpd->Write();
  TNamed("spec", spec).Write();
  fout.Close();
  std::cout << "wrote " << outfile << std::endl;
}
