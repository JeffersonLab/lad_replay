// lad_photon_sigma.C
//
// Width of the summed LAD photon flash for several timing calibrations, recomputed offline from the per-PMT
// arrays (X.ladhod.<plane>.Good{Top,Btm}TdcTimeUnCorr and ...AdcPulseAmp), as in lad_bar_timing_check.C:
//   t = 0.5 * [(t_top - tw_top) + (t_btm - tw_btm - 2 cableFit)] + LCoeff,  tw = c3 [(A/thr)^-c2 - (200/thr)^-c2]
//   tof - L/c = t - t_vertex - L/c (+ const), or t - t_vertex_RFcorr - z/c - L/c with the RF correction (as
//   THcLADKine::CalculateTOFRFcorr: t_vertex_RFcorr is the bunch time at the target centre), L from the vertex
//   (0, 0, z) to the hit (y from top-bottom).
// One histogram per parameter set x {no RF, RF} x hit selection, summed over all bars. Only bars that LADlib kept as a
// good hit (X.ladhod.goodhit_plane/paddle_0/1) in events with a vertex (react.ok, |z| < 20 cm) are used:
//   all     every such hit
//   back    back-plane hits (001, 101)
//   bveto   back-plane hits with no front PMT hit within +-1 paddle and 10 ns (photon-enriched)
//   bvetoE  bveto with E > 10 MeV (X.ladhod.<plane>.HodoHitEdep_MeV of that paddle)
// The veto uses the last set's times, so every set sees the same hits.
//
// params file: "<set> <plane> <side> <paddle> <c2> <c3> <cableFit> <LCoeff> <velFit>" and "<set> thr <value>"
// lines (lad_photon_sigma_params.py), sets 0..NSET-1.
//
// Usage: root -l -b -q 'lad_photon_sigma.C+("files.dat","out.root","params.txt","P",closure_set)'

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace lps {
const int NPL = 5, NSET = 4, NSEL = 4;
const char *PLANES[NPL]       = {"000", "001", "100", "101", "200"};
const double plane_zpos[NPL]  = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double plane_theta[NPL] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
const char *SELN[NSEL]        = {"all", "back", "bveto", "bvetoE"};
inline double paddle_centre(int pd) { return 110. - 22. * pd; } // pd from 0
const double C_LIGHT = 29.9792458, SENT = 1e9, EMIN = 10., ZMAX = 20.;

struct Par {
  double c2[2][NPL][11] = {}, c3[2][NPL][11] = {}, cab[NPL][11] = {}, lco[NPL][11] = {}, vel[NPL][11] = {};
  double thr = 120.;
};
} // namespace lps

void lad_photon_sigma(const char *listfile, const char *outfile, const char *parfile, const char *spec = "P",
                      int closure_set = lps::NSET - 1) {
  using namespace lps;
  Par par[NSET];
  for (int s = 0; s < NSET; s++)
    for (int p = 0; p < NPL; p++)
      for (int b = 0; b < 11; b++) {
        par[s].c3[0][p][b] = par[s].c3[1][p][b] = 1.;
        par[s].vel[p][b]                        = 16.2635;
      }
  {
    std::ifstream in(parfile);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream ss(line);
      int set;
      std::string pl;
      ss >> set >> pl;
      if (set < 0 || set >= NSET)
        continue;
      if (pl == "thr") {
        ss >> par[set].thr;
        continue;
      }
      int ip = -1;
      for (int i = 0; i < NPL; i++)
        if (pl == PLANES[i])
          ip = i;
      if (ip < 0)
        continue;
      std::string side;
      int pad;
      double c2, c3, cab, lco, vel;
      ss >> side >> pad >> c2 >> c3 >> cab >> lco >> vel;
      const int is = side == "Top" ? 0 : 1;
      par[set].c2[is][ip][pad - 1] = c2;
      par[set].c3[is][ip][pad - 1] = c3;
      par[set].cab[ip][pad - 1]    = cab;
      par[set].lco[ip][pad - 1]    = lco;
      par[set].vel[ip][pad - 1]    = vel;
    }
  }
  TChain ch("T");
  {
    std::ifstream in(listfile);
    std::string line;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#')
        ch.Add(line.c_str());
  }
  std::cout << "entries " << ch.GetEntries() << std::endl;
  TTreeReader rd(&ch);
  TTreeReaderValue<double> tv(rd, Form("%s.ladkin.t_vertex", spec));
  TTreeReaderValue<double> tvrf(rd, Form("%s.ladkin.t_vertex_RFcorr", spec));
  TTreeReaderValue<double> vok(rd, Form("%s.react.ok", spec)), vz(rd, Form("%s.react.z", spec));
  TTreeReaderArray<double> gpl0(rd, Form("%s.ladhod.goodhit_plane_0", spec)), gpl1(rd, Form("%s.ladhod.goodhit_plane_1", spec));
  TTreeReaderArray<double> gpd0(rd, Form("%s.ladhod.goodhit_paddle_0", spec)), gpd1(rd, Form("%s.ladhod.goodhit_paddle_1", spec));
  std::vector<std::unique_ptr<TTreeReaderArray<double>>> tt, tb, at, ab, avg, hpad, hE;
  for (int p = 0; p < NPL; p++) {
    tt.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodTopTdcTimeUnCorr", spec, PLANES[p])));
    tb.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodBtmTdcTimeUnCorr", spec, PLANES[p])));
    at.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodTopAdcPulseAmp", spec, PLANES[p])));
    ab.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodBtmAdcPulseAmp", spec, PLANES[p])));
    avg.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodHitTimeAvg", spec, PLANES[p])));
    hpad.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.HodoHitPaddleNum", spec, PLANES[p])));
    hE.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.HodoHitEdep_MeV", spec, PLANES[p])));
  }

  TH1::AddDirectory(kFALSE);
  TH1D *h[NSET][2][NSEL];
  for (int s = 0; s < NSET; s++)
    for (int r = 0; r < 2; r++)
      for (int c = 0; c < NSEL; c++)
        h[s][r][c] = new TH1D(Form("h_s%d_%s_%s", s, r ? "rf" : "norf", SELN[c]),
                              Form("set %d %s %s;t_{hit} - t_{vertex} - L/c (ns)", s, r ? "RF" : "no RF", SELN[c]),
                              8000, 1600, 2000);
  TH1D *hclos = new TH1D("h_closure", "recomputed - GoodHitTimeAvg;ns", 2000, -10, 10);

  auto twf = [](double A, double c2, double c3, double thr) {
    return c3 * (1. / std::pow(A / thr, c2) - 1. / std::pow(200. / thr, c2));
  };
  double t[NSET][NPL][11], E[NPL][11];
  bool ok[NPL][11], good[NPL][11];
  while (rd.Next()) {
    const bool okv = std::fabs(*tv) < SENT, okrf = std::fabs(*tvrf) < SENT;
    if ((!okv && !okrf) || *vok == 0 || !(std::fabs(*vz) < ZMAX))
      continue;
    for (int p = 0; p < NPL; p++)
      for (int b = 0; b < 11; b++)
        good[p][b] = false;
    auto mark = [&](const TTreeReaderArray<double> &pl, const TTreeReaderArray<double> &pd) {
      for (size_t i = 0; i < pl.GetSize() && i < pd.GetSize(); i++) {
        const int p = (int)std::lround(pl[i]), b = (int)std::lround(pd[i]); // plane index, 0-based paddle
        if (p >= 0 && p < NPL && b >= 0 && b < 11)
          good[p][b] = true;
      }
    };
    mark(gpl0, gpd0);
    mark(gpl1, gpd1);
    for (int p = 0; p < NPL; p++) {
      for (int b = 0; b < 11; b++) {
        E[p][b]  = -1.;
        ok[p][b] = false;
      }
      for (size_t i = 0; i < hpad[p]->GetSize() && i < hE[p]->GetSize(); i++) {
        const int b = (int)std::lround((*hpad[p])[i]) - 1; // HodoHitPaddleNum is 1-based
        if (b >= 0 && b < 11)
          E[p][b] = std::max(E[p][b], (*hE[p])[i]);
      }
      for (int b = 0; b < 11; b++) {
        if ((size_t)b >= tt[p]->GetSize())
          continue;
        const double t1 = (*tt[p])[b], t2 = (*tb[p])[b], a1 = (*at[p])[b], a2 = (*ab[p])[b];
        if (!(std::fabs(t1) < SENT && std::fabs(t2) < SENT && a1 > 0 && a2 > 0))
          continue;
        ok[p][b] = true;
        const double x0 = paddle_centre(b), z0 = plane_zpos[p];
        const double cth = std::cos(plane_theta[p]), sth = std::sin(plane_theta[p]);
        const double X = x0 * cth + z0 * sth, Z = -x0 * sth + z0 * cth;
        for (int s = 0; s < NSET; s++) {
          const Par &q    = par[s];
          const double ct = t1 - twf(a1, q.c2[0][p][b], q.c3[0][p][b], q.thr) + q.lco[p][b];
          const double cb = t2 - twf(a2, q.c2[1][p][b], q.c3[1][p][b], q.thr) - 2 * q.cab[p][b] + q.lco[p][b];
          const double y  = std::max(-220., std::min(220., 0.5 * (cb - ct) * q.vel[p][b]));
          const double dz = Z - *vz; // path from the vertex
          t[s][p][b]      = 0.5 * (ct + cb) - std::sqrt(X * X + y * y + dz * dz) / C_LIGHT;
          if (s == closure_set && (size_t)b < avg[p]->GetSize() && std::fabs((*avg[p])[b]) < SENT)
            hclos->Fill(0.5 * (ct + cb) - (*avg[p])[b]);
        }
      }
    }
    for (int p = 0; p < NPL; p++)
      for (int b = 0; b < 11; b++) {
        if (!ok[p][b] || !good[p][b])
          continue;
        const bool back = (p == 1 || p == 3);
        bool vetoed     = false;
        if (back)
          for (int bf = std::max(0, b - 1); bf <= std::min(10, b + 1) && !vetoed; bf++)
            if (ok[p - 1][bf] && std::fabs(t[NSET - 1][p - 1][bf] - t[NSET - 1][p][b]) < 10.)
              vetoed = true;
        const bool sel[NSEL] = {true, back, back && !vetoed, back && !vetoed && E[p][b] > EMIN};
        for (int s = 0; s < NSET; s++)
          for (int r = 0; r < 2; r++) {
            if (!(r ? okrf : okv))
              continue;
            const double v = t[s][p][b] - (r ? *tvrf + *vz / C_LIGHT : *tv);
            for (int c = 0; c < NSEL; c++)
              if (sel[c])
                h[s][r][c]->Fill(v);
          }
      }
  }
  TFile fout(outfile, "RECREATE");
  hclos->Write();
  for (int s = 0; s < NSET; s++)
    for (int r = 0; r < 2; r++)
      for (int c = 0; c < NSEL; c++)
        h[s][r][c]->Write();
  fout.Close();
  std::cout << "wrote " << outfile << std::endl;
}
