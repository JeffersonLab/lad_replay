// lad_edep_histos.C
//
// Energy-deposition calibration histograms for the LAD hodoscope, per PMT-HV period, plane and paddle:
//   1/beta (from the RF-corrected ToF and the hit's path length) vs geometric-mean pulse integral
//   1/beta                                                     vs geometric-mean pulse amplitude
// for hits close to normal incidence (path through the bar < MAX_PATH_RATIO * thickness), for all hits,
// for hits that are part of a front-back pair, and for good hits matched to a GEM track. 1/beta runs to negative values (ToF before the photon flash)
// so that the accidental ADC spectrum can be measured per bar and subtracted.
// Stopping protons deposit their full kinetic energy (edep rises with T); above the punch-through
// energy edep falls again. lad_edep_fit.py fits that ridge with a range-energy model for every bar,
// including plane 200 (no back plane). MIPs (pions) sit at 1/beta ~ 1 in the same histograms.
//
// Also Delta-E - E histograms for same-paddle front-back pairs (front ADC vs back ADC per back paddle), which
// calibrate the back planes against the front ones without ToF.
//
// Hits are taken from the SHMS reference (P.ladhod) when the event has an SHMS vertex time, else from
// the HMS reference, so coincidence events are not counted twice.
//
// Optional exclusion list: lines "<plane> <paddle> <run>" (e.g. "100 1 23038"); hits of that bar in that run are
// skipped, and so are front-back pairs containing it (bars with a hardware fault in some runs).
//
// Optional path window: keep hits with path_lo <= path/thickness < path_hi (default: the MAX_PATH_RATIO cut), e.g.
// 1.03, 1.5 for the oblique hits the default cut drops.
//
// Usage: root -l -b -q 'lad_edep_histos.C+("files.dat","out.root",4,"exclude.txt")'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TROOT.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../lad_tof/lad_tof_offset.h"

namespace led {
const int N_PLANES                   = 5;
const int N_PADDLES                  = 11;
const char *PLANES[N_PLANES]         = {"000", "001", "100", "101", "200"};
const double plane_zpos[N_PLANES]    = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double plane_theta[N_PLANES]   = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
inline double paddle_centre(int pd) { return 110. - 22. * pd; }
const double C_LIGHT  = 29.9792458;
const double SENT     = 1e9;
const double VTX_ZMAX = 20.;
const double MAX_PATH_RATIO = 1.03; // keep hits with path/thickness < MAX_PATH_RATIO (incidence < ~14 deg)

// PMT HV periods (user, 2026-10-08)
const int N_PER                 = 6;
const int PER_START[N_PER]      = {0, 22375, 22533, 22536, 22695, 22698};
const char *PER_NAME[N_PER]     = {"p0_lt22375", "p1_22375", "p2_22533", "p3_22536", "p4_22695", "p5_22698"};
std::set<long> EXCL; // plane * 100 + paddle (0-based) + run * 1000
inline bool excluded(double run, int p, int b) { return !EXCL.empty() && EXCL.count((long)run * 1000 + p * 100 + b); }
inline int period(double run) {
  int p = 0;
  for (int i = 0; i < N_PER; i++)
    if (run >= PER_START[i])
      p = i;
  return p;
}

const int NB_IB = 280;
const double IB_LO = -2.0, IB_HI = 5.0; // 1/beta; 1/beta < ~0.9 is before the photon flash: accidentals only
const int NB_DT = 160;
const double DT_LO = -4.0, DT_HI = 12.0; // t_back - t_front (ns)
const int NB_INT = 250;
const double INT_HI = 250.; // pulse integral (pC)
const int NB_AMP = 250;
const double AMP_HI = 1000.; // mV

// categories: 0 all hits, 1 hits of a front-back pair (both goodhit planes valid), 2 good hits matched to a
// GEM track (goodhit_trackid >= 0)
const int N_CAT             = 3;
const char *CAT_NAME[N_CAT] = {"all", "pair", "trk"};

struct Hists {
  std::array<std::array<std::array<std::array<TH2F *, N_PADDLES>, N_PLANES>, N_PER>, N_CAT> hint{}, hamp{};
  std::array<std::array<TH1D *, N_PLANES>, N_PER> hpath{};
  // Delta-E - E: same-paddle front-back pairs, front ADC vs back ADC, per back plane (0: 001, 1: 101) and paddle
  std::array<std::array<std::array<TH2F *, N_PADDLES>, 2>, N_PER> hdee_int{}, hdee_amp{};
  // front-back Delta t (t_back - t_front) vs back ADC and vs front ADC, same-paddle pairs
  std::array<std::array<std::array<TH2F *, N_PADDLES>, 2>, N_PER> hdtb_int{}, hdtb_amp{}, hdtf_int{}, hdtf_amp{};
  void book(unsigned slot) {
    for (int ip = 0; ip < N_PER; ip++)
      for (int p = 0; p < N_PLANES; p++) {
        hpath[ip][p] = new TH1D(Form("h_path_%s_%s_s%u", PER_NAME[ip], PLANES[p], slot), ";path / thickness", 100, 1, 1.5);
        hpath[ip][p]->SetDirectory(nullptr);
        if (p == 1 || p == 3)
          for (int b = 0; b < N_PADDLES; b++) {
            const int k = p == 1 ? 0 : 1;
            hdee_int[ip][k][b] = new TH2F(Form("h_dee_int_%s_%s_%d_s%u", PER_NAME[ip], PLANES[p], b + 1, slot),
                                          Form("%s back plane %s paddle %d;front int (pC);back int (pC)", PER_NAME[ip],
                                               PLANES[p], b + 1),
                                          NB_INT, 0, INT_HI, NB_INT, 0, INT_HI);
            hdee_amp[ip][k][b] = new TH2F(Form("h_dee_amp_%s_%s_%d_s%u", PER_NAME[ip], PLANES[p], b + 1, slot),
                                          Form("%s back plane %s paddle %d;front amp (mV);back amp (mV)", PER_NAME[ip],
                                               PLANES[p], b + 1),
                                          NB_AMP, 0, AMP_HI, NB_AMP, 0, AMP_HI);
            hdee_int[ip][k][b]->SetDirectory(nullptr);
            hdee_amp[ip][k][b]->SetDirectory(nullptr);
            auto mk = [&](const char *nm, const char *yt, int ny, double yhi) {
              TH2F *hh = new TH2F(Form("h_%s_%s_%s_%d_s%u", nm, PER_NAME[ip], PLANES[p], b + 1, slot),
                                  Form("%s back plane %s paddle %d;t_{back} - t_{front} (ns);%s", PER_NAME[ip], PLANES[p],
                                       b + 1, yt),
                                  NB_DT, DT_LO, DT_HI, ny, 0, yhi);
              hh->SetDirectory(nullptr);
              return hh;
            };
            hdtb_int[ip][k][b] = mk("dtb_int", "back int (pC)", NB_INT, INT_HI);
            hdtb_amp[ip][k][b] = mk("dtb_amp", "back amp (mV)", NB_AMP, AMP_HI);
            hdtf_int[ip][k][b] = mk("dtf_int", "front int (pC)", NB_INT, INT_HI);
            hdtf_amp[ip][k][b] = mk("dtf_amp", "front amp (mV)", NB_AMP, AMP_HI);
          }
        for (int b = 0; b < N_PADDLES; b++)
          for (int c = 0; c < N_CAT; c++) {
            hint[c][ip][p][b] = new TH2F(Form("h_int_%s_%s_%s_%d_s%u", CAT_NAME[c], PER_NAME[ip], PLANES[p], b + 1, slot),
                                         Form("%s %s plane %s paddle %d;1/#beta;#sqrt{int_{top} int_{btm}} (pC)",
                                              CAT_NAME[c], PER_NAME[ip], PLANES[p], b + 1),
                                         NB_IB, IB_LO, IB_HI, NB_INT, 0, INT_HI);
            hamp[c][ip][p][b] = new TH2F(Form("h_amp_%s_%s_%s_%d_s%u", CAT_NAME[c], PER_NAME[ip], PLANES[p], b + 1, slot),
                                         Form("%s %s plane %s paddle %d;1/#beta;#sqrt{amp_{top} amp_{btm}} (mV)",
                                              CAT_NAME[c], PER_NAME[ip], PLANES[p], b + 1),
                                         NB_IB, IB_LO, IB_HI, NB_AMP, 0, AMP_HI);
            hint[c][ip][p][b]->SetDirectory(nullptr);
            hamp[c][ip][p][b]->SetDirectory(nullptr);
          }
      }
  }
};
} // namespace led

void lad_edep_histos(const char *listfile, const char *outfile, int nthreads = 4, const char *exclfile = "",
                     double path_lo = 0., double path_hi = led::MAX_PATH_RATIO) {
  using namespace led;
  using RVd = ROOT::VecOps::RVec<double>;
  TH1::AddDirectory(kFALSE);
  if (exclfile && exclfile[0]) {
    std::ifstream in(exclfile);
    std::string line, pl;
    int pd;
    long run;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream ss(line);
      if (!(ss >> pl >> pd >> run))
        continue;
      for (int p = 0; p < N_PLANES; p++)
        if (pl == PLANES[p])
          EXCL.insert(run * 1000 + p * 100 + (pd - 1));
    }
    std::cout << "excluding " << EXCL.size() << " (bar, run) combinations from " << exclfile << std::endl;
  }
  TChain ch("T");
  {
    std::ifstream in(listfile);
    std::string line;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#')
        ch.Add(line.c_str());
  }
  if (nthreads > 1)
    ROOT::EnableImplicitMT(nthreads);
  ROOT::RDataFrame df(ch);
  const unsigned nslots = df.GetNSlots();
  std::vector<std::unique_ptr<Hists>> H;
  for (unsigned s = 0; s < nslots; s++) {
    H.emplace_back(new Hists);
    H.back()->book(s);
  }

  for (char spec : {'P', 'H'}) {
    const std::string sp(1, spec);
    const std::string g = sp + ".ladhod.goodhit_";
    ROOT::RDF::RNode d = df;
    if (spec == 'H')
      d = d.Filter("std::fabs(P.ladkin.t_vertex) >= 1e9", "no SHMS vertex time");
    d = ladtof::define_tof(d, sp + "_tofrf0", spec, "0", true);
    d = ladtof::define_tof(d, sp + "_tofrf1", spec, "1", true);
    std::vector<std::string> cols = {g + "plane_0",    g + "plane_1",      g + "paddle_0",     g + "paddle_1",
                                     sp + "_tofrf0",   sp + "_tofrf1",     g + "hit_ypos_0",   g + "hit_ypos_1",
                                     g + "hitedep_0",  g + "hitedep_1",    g + "hitedep_amp_0", g + "hitedep_amp_1",
                                     g + "trackid",    g + "hittime_0",    g + "hittime_1",
                                     sp + ".react.ok", sp + ".react.z",    "g.runnum"};
    d.ForeachSlot(
        [&H, spec, path_lo, path_hi](unsigned slot, const RVd &pl0, const RVd &pl1, const RVd &pd0, const RVd &pd1, const RVd &t0,
             const RVd &t1, const RVd &y0, const RVd &y1, const RVd &e0, const RVd &e1, const RVd &a0, const RVd &a1,
             const RVd &trk, const RVd &ht0, const RVd &ht1, double vok, double vz, double run) {
          Hists &h       = *H[slot];
          const int per  = period(run);
          // Delta-E - E pairs (no ToF needed; SHMS-reference pass only, so nothing is filled twice)
          for (size_t i = 0; spec == 'P' && i < pl0.size() && i < pl1.size(); i++) {
            const int pf = (int)pl0[i], pb = (int)pl1[i];
            if (!((pf == 0 && pb == 1) || (pf == 2 && pb == 3)) || i >= pd0.size() || i >= pd1.size() ||
                pd0[i] != pd1[i])
              continue;
            const int b = (int)pd1[i];
            if (b < 0 || b >= N_PADDLES || excluded(run, pf, b) || excluded(run, pb, b))
              continue;
            const int k = pb == 1 ? 0 : 1;
            if (i < e0.size() && i < e1.size() && std::fabs(e0[i]) < SENT && std::fabs(e1[i]) < SENT)
              h.hdee_int[per][k][b]->Fill(e0[i], e1[i]);
            if (i < a0.size() && i < a1.size() && std::fabs(a0[i]) < SENT && std::fabs(a1[i]) < SENT)
              h.hdee_amp[per][k][b]->Fill(a0[i], a1[i]);
            if (i < ht0.size() && i < ht1.size() && std::fabs(ht0[i]) < SENT && std::fabs(ht1[i]) < SENT) {
              const double dt = ht1[i] - ht0[i];
              if (i < e0.size() && i < e1.size() && std::fabs(e0[i]) < SENT && std::fabs(e1[i]) < SENT) {
                h.hdtb_int[per][k][b]->Fill(dt, e1[i]);
                h.hdtf_int[per][k][b]->Fill(dt, e0[i]);
              }
              if (i < a0.size() && i < a1.size() && std::fabs(a0[i]) < SENT && std::fabs(a1[i]) < SENT) {
                h.hdtb_amp[per][k][b]->Fill(dt, a1[i]);
                h.hdtf_amp[per][k][b]->Fill(dt, a0[i]);
              }
            }
          }
          const double z = (vok != 0 && std::fabs(vz) < VTX_ZMAX) ? vz : 0.;
          for (int side = 0; side < 2; side++) {
            const RVd &PL = side ? pl1 : pl0, &PD = side ? pd1 : pd0, &T = side ? t1 : t0, &Y = side ? y1 : y0,
                      &E = side ? e1 : e0, &A = side ? a1 : a0;
            for (size_t i = 0; i < PL.size() && i < T.size(); i++) {
              const int p = (int)PL[i], b = (int)PD[i];
              if (p < 0 || p >= N_PLANES || b < 0 || b >= N_PADDLES || !(std::fabs(T[i]) < SENT) || excluded(run, p, b))
                continue;
              // pair: the other plane of this good hit is also there
              const RVd &PLo = side ? pl0 : pl1;
              const bool pair = i < PLo.size() && PLo[i] >= 0 && PLo[i] < N_PLANES;
              const double y  = std::fabs(Y[i]) < SENT ? Y[i] : 0.;
              const double x0 = paddle_centre(b), z0 = plane_zpos[p];
              const double c = std::cos(plane_theta[p]), s = std::sin(plane_theta[p]);
              const double X = x0 * c + z0 * s, Zz = -x0 * s + z0 * c - z;
              const double L = std::sqrt(X * X + y * y + Zz * Zz);
              // incidence: angle between the hit direction and the plane normal (sin th, 0, cos th)
              const double cosa = (X * s + Zz * c) / L;
              const double path = 1. / std::fabs(cosa);
              h.hpath[per][p]->Fill(path);
              if (path < path_lo || path >= path_hi)
                continue;
              const double ib = T[i] * C_LIGHT / L;
              const bool intrk = i < trk.size() && trk[i] >= 0;
              for (int cat = 0; cat < N_CAT; cat++) {
                if ((cat == 1 && !pair) || (cat == 2 && !intrk))
                  continue;
                if (i < E.size() && std::fabs(E[i]) < SENT)
                  h.hint[cat][per][p][b]->Fill(ib, E[i]);
                if (i < A.size() && std::fabs(A[i]) < SENT)
                  h.hamp[cat][per][p][b]->Fill(ib, A[i]);
              }
            }
          }
        },
        cols);
  }

  TFile fout(outfile, "RECREATE");
  for (int ip = 0; ip < N_PER; ip++) {
    bool any = false;
    for (unsigned s = 0; s < nslots && !any; s++)
      for (int p = 0; p < N_PLANES && !any; p++)
        any = H[s]->hpath[ip][p]->GetEntries() > 0;
    if (!any)
      continue;
    fout.mkdir(PER_NAME[ip])->cd();
    for (int p = 0; p < N_PLANES; p++) {
      TH1D *hp = (TH1D *)H[0]->hpath[ip][p]->Clone(Form("h_path_%s", PLANES[p]));
      for (unsigned s = 1; s < nslots; s++)
        hp->Add(H[s]->hpath[ip][p]);
      hp->Write();
      if (p == 1 || p == 3)
        for (int b = 0; b < N_PADDLES; b++) {
          const int k = p == 1 ? 0 : 1;
          TH2F *di = (TH2F *)H[0]->hdee_int[ip][k][b]->Clone(Form("h_dee_int_%s_%d", PLANES[p], b + 1));
          TH2F *da = (TH2F *)H[0]->hdee_amp[ip][k][b]->Clone(Form("h_dee_amp_%s_%d", PLANES[p], b + 1));
          for (unsigned s = 1; s < nslots; s++) {
            di->Add(H[s]->hdee_int[ip][k][b]);
            da->Add(H[s]->hdee_amp[ip][k][b]);
          }
          di->Write();
          da->Write();
          delete di;
          delete da;
          for (auto *arr : {&H[0]->hdtb_int, &H[0]->hdtb_amp, &H[0]->hdtf_int, &H[0]->hdtf_amp}) {
            const std::string base = std::string((*arr)[ip][k][b]->GetName());
            TH2F *hs = (TH2F *)(*arr)[ip][k][b]->Clone(base.substr(0, base.find("_" + std::string(PER_NAME[ip]))).c_str());
            hs->SetName(Form("%s_%s_%d", hs->GetName(), PLANES[p], b + 1));
            const size_t off = (char *)arr - (char *)H[0].get();
            for (unsigned s = 1; s < nslots; s++)
              hs->Add((*(decltype(arr))((char *)H[s].get() + off))[ip][k][b]);
            hs->Write();
            delete hs;
          }
        }
      for (int b = 0; b < N_PADDLES; b++)
        for (int c = 0; c < N_CAT; c++) {
          TH2F *hi = (TH2F *)H[0]->hint[c][ip][p][b]->Clone(Form("h_int_%s_%s_%d", CAT_NAME[c], PLANES[p], b + 1));
          TH2F *ha = (TH2F *)H[0]->hamp[c][ip][p][b]->Clone(Form("h_amp_%s_%s_%d", CAT_NAME[c], PLANES[p], b + 1));
          for (unsigned s = 1; s < nslots; s++) {
            hi->Add(H[s]->hint[c][ip][p][b]);
            ha->Add(H[s]->hamp[c][ip][p][b]);
          }
          hi->Write();
          ha->Write();
          delete hi;
          delete ha;
        }
    }
  }
  fout.cd();
  TNamed("tof", ladtof::signature().c_str()).Write();
  TNamed("path", Form("%g,%g", path_lo, path_hi)).Write();
  fout.Close();
  std::cout << "wrote " << outfile << std::endl;
}
