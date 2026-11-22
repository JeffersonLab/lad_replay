// lad_bar_timing_check.C
//
// Per-bar photon-flash check of a LAD timing calibration, recomputed offline from the per-PMT arrays
// (X.ladhod.<plane>.Good{Top,Btm}TdcTimeUnCorr and ...AdcPulseAmp) with two parameter sets:
//   set 0 ("old"): should reproduce the replay (closure: compared to X.ladhod.<plane>.GoodHitTimeAvg)
//   set 1 ("new"): the calibration under test
// For each set the hit time is
//   t = 0.5 * [(t_top - tw_top) + (t_btm - tw_btm - 2 cableFit)] + LCoeff,  tw = c3 [(A/thr)^-c2 - (200/thr)^-c2]
// and tof - L/c = t - t_vertex_RFcorr - L/c (+ const), with L from the target centre to the hit
// (y from the top-bottom difference and velFit). Histograms per plane/paddle, for all hits and for back-plane
// hits with no front hit within +-1 paddle and 10 ns (photon-enriched).
//
// params file: one line per PMT: "<set> <plane> <side> <paddle> <c2> <c3> <cableFit> <LCoeff> <velFit>",
// plus "<set> thr <value>" lines (made by lad_bar_timing_params.py).
//
// Usage: root -l -b -q 'lad_bar_timing_check.C+("files.dat","out.root",1,"params.txt","P")'

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace lbt {
const int NPL                = 5;
const char *PLANES[NPL]      = {"000", "001", "100", "101", "200"};
const double plane_zpos[NPL] = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double plane_theta[NPL] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
inline double paddle_centre(int pd) { return 110. - 22. * pd; } // pd from 0
const double C_LIGHT = 29.9792458;
const double SENT    = 1e9;

struct Par {
  double c2[2][NPL][11] = {}, c3[2][NPL][11] = {}, cab[NPL][11] = {}, lco[NPL][11] = {}, vel[NPL][11] = {};
  double thr = 120.;
};
} // namespace lbt

void lad_bar_timing_check(const char *listfile, const char *outfile, int /*nthreads*/ = 1, const char *parfile = "",
                          const char *spec = "P") {
  using namespace lbt;
  Par par[2];
  for (int s = 0; s < 2; s++)
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
      if (set < 0 || set > 1)
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
      int is = side == "Top" ? 0 : 1;
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
  std::vector<std::unique_ptr<TTreeReaderArray<double>>> tt, tb, at, ab, avg;
  for (int p = 0; p < NPL; p++) {
    tt.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodTopTdcTimeUnCorr", spec, PLANES[p])));
    tb.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodBtmTdcTimeUnCorr", spec, PLANES[p])));
    at.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodTopAdcPulseAmp", spec, PLANES[p])));
    ab.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodBtmAdcPulseAmp", spec, PLANES[p])));
    avg.emplace_back(new TTreeReaderArray<double>(rd, Form("%s.ladhod.%s.GoodHitTimeAvg", spec, PLANES[p])));
  }

  TH1::AddDirectory(kFALSE);
  const char *SETN[2] = {"old", "new"};
  const char *CATN[2] = {"all", "bveto"};
  TH1D *h[2][2][NPL][11];
  TH1D *hclos = new TH1D("h_closure", "old recomputed - GoodHitTimeAvg;ns", 2000, -10, 10);
  TH2D *hy[2][NPL];
  // photon time vs amplitude (geometric mean of top and bottom) and vs run, per plane; back planes use the
  // front-vetoed hits, the other planes all hits
  TH2D *hamp[2][NPL], *hrun[2][NPL];
  for (int s = 0; s < 2; s++)
    for (int p = 0; p < NPL; p++) {
      hy[s][p] = new TH2D(Form("h_y_%s_%s", SETN[s], PLANES[p]), ";paddle;y (cm)", 11, 0.5, 11.5, 120, -300, 300);
      hamp[s][p] = new TH2D(Form("h_amp_%s_%s", SETN[s], PLANES[p]), ";#sqrt{A_{top}A_{btm}} (mV);t (ns)", 100, 0, 500,
                            600, 1710, 1800);
      hrun[s][p] = new TH2D(Form("h_run_%s_%s", SETN[s], PLANES[p]), ";run;t (ns)", 1450, 22450, 23900, 180, 1710, 1800);
      for (int c = 0; c < 2; c++)
        for (int b = 0; b < 11; b++)
          h[s][c][p][b] = new TH1D(Form("h_%s_%s_%s_%d", SETN[s], CATN[c], PLANES[p], b + 1),
                                   Form("%s %s plane %s paddle %d;t_{hit} - t_{vertex,RF} - L/c (ns)", SETN[s], CATN[c],
                                        PLANES[p], b + 1),
                                   4000, 1600, 2000);
    }

  auto twf = [](double A, double c2, double c3, double thr) {
    return c3 * (1. / std::pow(A / thr, c2) - 1. / std::pow(200. / thr, c2));
  };

  double tnew[2][NPL][11], ageo[NPL][11];
  bool ok[NPL][11];
  TTreeReaderValue<double> runv(rd, "g.runnum");
  while (rd.Next()) {
    if (!(std::fabs(*tvrf) < SENT))
      continue;
    for (int p = 0; p < NPL; p++)
      for (int b = 0; b < 11; b++) {
        ok[p][b] = false;
        if ((size_t)b >= tt[p]->GetSize())
          continue;
        double t1 = (*tt[p])[b], t2 = (*tb[p])[b], a1 = (*at[p])[b], a2 = (*ab[p])[b];
        if (!(std::fabs(t1) < SENT && std::fabs(t2) < SENT && a1 > 0 && a2 > 0))
          continue;
        ok[p][b]   = true;
        ageo[p][b] = std::sqrt(a1 * a2);
        for (int s = 0; s < 2; s++) {
          const Par &q = par[s];
          double ct = t1 - twf(a1, q.c2[0][p][b], q.c3[0][p][b], q.thr) + q.lco[p][b];
          double cb = t2 - twf(a2, q.c2[1][p][b], q.c3[1][p][b], q.thr) - 2 * q.cab[p][b] + q.lco[p][b];
          double y  = 0.5 * (cb - ct) * q.vel[p][b];
          tnew[s][p][b] = 0.5 * (ct + cb);
          // path from the target centre
          const double x0 = paddle_centre(b), z0 = plane_zpos[p];
          const double c = std::cos(plane_theta[p]), sn = std::sin(plane_theta[p]);
          const double X = x0 * c + z0 * sn, Z = -x0 * sn + z0 * c;
          const double yy = std::max(-220., std::min(220., y));
          tnew[s][p][b] -= std::sqrt(X * X + yy * yy + Z * Z) / C_LIGHT;
          hy[s][p]->Fill(b + 1, y);
        }
        if ((size_t)b < avg[p]->GetSize() && std::fabs((*avg[p])[b]) < SENT) {
          const Par &q = par[0];
          double ct = t1 - twf(a1, q.c2[0][p][b], q.c3[0][p][b], q.thr) + q.lco[p][b];
          double cb = t2 - twf(a2, q.c2[1][p][b], q.c3[1][p][b], q.thr) - 2 * q.cab[p][b] + q.lco[p][b];
          hclos->Fill(0.5 * (ct + cb) - (*avg[p])[b]);
        }
      }
    for (int p = 0; p < NPL; p++)
      for (int b = 0; b < 11; b++) {
        if (!ok[p][b])
          continue;
        bool back = (p == 1 || p == 3), vetoed = false;
        if (back)
          for (int bf = std::max(0, b - 1); bf <= std::min(10, b + 1) && !vetoed; bf++)
            if (ok[p - 1][bf] && std::fabs(tnew[1][p - 1][bf] - tnew[1][p][b]) < 10.)
              vetoed = true;
        for (int s = 0; s < 2; s++) {
          double v = tnew[s][p][b] - *tvrf;
          h[s][0][p][b]->Fill(v);
          if (back && !vetoed)
            h[s][1][p][b]->Fill(v);
          if (!back || !vetoed) {
            hamp[s][p]->Fill(ageo[p][b], v);
            hrun[s][p]->Fill(*runv, v);
          }
        }
      }
  }
  TFile fout(outfile, "RECREATE");
  hclos->Write();
  for (int s = 0; s < 2; s++)
    for (int p = 0; p < NPL; p++) {
      hy[s][p]->Write();
      hamp[s][p]->Write();
      hrun[s][p]->Write();
      for (int c = 0; c < 2; c++)
        for (int b = 0; b < 11; b++)
          h[s][c][p][b]->Write();
    }
  fout.Close();
  std::cout << "wrote " << outfile << std::endl;
}
