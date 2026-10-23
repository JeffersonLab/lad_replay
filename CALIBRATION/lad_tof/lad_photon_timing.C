// lad_photon_timing.C
//
// LAD timing calibration input: photon (tof - L/c) peak per bar, per plane and
// in total, for both spectrometer references (P = SHMS vertex time, H = HMS).
// Used to determine lglobal_time_offset (PARAM/LAD/LADKINE/lladkine.param).
//
// No GEM / LAD tracking requirement. Every LAD good hit is used (front-back
// paired and single-plane), with the hodoscope time = top/bottom average after
// all timing corrections (THcLADHodoHit::GetScinCorrectedTime), i.e.
// X.ladhod.goodhit_hit_tof_{0,1} = hittime - t_vertex + lglobal_time_offset.
// The ToF is re-referenced hit by hit to the calibrated per-spectrometer offsets in
// lad_tof_offset.h (any replay gives the same result); the offsets used are stored
// in the output as TNamed "tof_offsets" for lad_photon_timing_fit.py.
//
// Hit categories ("cat"):
//   all     every hit, any plane (000,001,100,101,200)
//   back    back-plane hits only (001,101), no veto
//   back_fv back-plane hits with a front veto: no hit in the matching front
//           plane (000 for 001, 100 for 101) within +-FV_DPAD paddles and
//           |dt| < FV_DT ns (FV_DT<=0: any time)
// Variables:
//   tof     goodhit_hit_tof (calibrated)         (ns)
//   tofc    tof - L/c                            (ns), photon peak -> 0 when calibrated
//   tofcrf  goodhit_hit_tof_rfcorr (calib.) - L/c (ns), RF-corrected vertex time
// L = |hit_lab - vertex|: hit_lab from lhodo_geom.param (paddle centre, ypos,
// plane zpos rotated by plane theta); vertex = (0,0,X.react.z) when
// X.react.ok && |z| < VTX_ZMAX, else the target centre.
//
// Output (histograms only, hadd-able):
//   <X>/<cat>/<var>_p<plane>_b<paddle>, <var>_p<plane>_sum, <var>_total
//   <X>/<cat>/<var>_vs_edep_p<plane>_b<paddle>, ..._sum, ..._total   (cats back, back_fv; tofc and tofcrf)
//   <X>/run/<var>_vs_run_<cat>                                        (tofc, tofcrf; cats all, back_fv)
//
// Usage:
//   root -l -b -q 'lad_photon_timing.C+("list.dat","out.root",4)'
// list.dat: one ROOT file per line (# comments allowed).

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDirectory.h>
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
#include <string>
#include <vector>

#include "lad_tof_offset.h" // calibrated LAD ToF convention (photon peak at tof-L/c = 0)

namespace lpt {

const int N_PLANES = 5, N_PADDLES = 11, N_CATS = 3, N_VARS = 3;
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const char *const cat_names[N_CATS]     = {"all", "back", "back_fv"};
const char *const var_names[N_VARS]     = {"tof", "tofc", "tofcrf"};

// lhodo_geom.param
const double plane_zpos[N_PLANES]  = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double plane_theta[N_PLANES] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
inline double paddle_centre(int pd) { return 110. - 22. * pd; }
const double C_LIGHT = 29.9792458; // cm/ns

const double VTX_ZMAX = 20.; // cm
const int FV_DPAD     = 1;   // front veto: +- paddles
const double FV_DT    = 10.; // front veto: |dt| window (ns); <=0 = any time

// Binning
const int NB_T = 3000;
const double T_LO[N_VARS] = {-50., -100., -100.};
const double T_HI[N_VARS] = {250., 200., 200.};
const int NB_T2   = 600;  const double T2_LO = -60., T2_HI = 60.; // 2D tofc axis (0.2 ns), photon peak at 0
const int NB_E    = 100;  const double E_LO = 0., E_HI = 100.;   // MeV
const int RUN_LO  = 22500, RUN_HI = 23800;

const double SENT = 1e9; // |value| above this = unset (1e30 / 1e38 sentinels)

struct Hists {
  // [cat][var][plane][paddle], [cat][var][plane], [cat][var]
  std::array<std::array<std::array<std::array<TH1D *, N_PADDLES>, N_PLANES>, N_VARS>, N_CATS> bar{};
  std::array<std::array<std::array<TH1D *, N_PLANES>, N_VARS>, N_CATS> pl{};
  std::array<std::array<TH1D *, N_VARS>, N_CATS> tot{};
  // 2D vs edep: cats back(1), back_fv(2); vars tofc(1), tofcrf(2); indexed [cat-1][var-1]
  std::array<std::array<std::array<std::array<TH2F *, N_PADDLES>, N_PLANES>, 2>, 2> e_bar{};
  std::array<std::array<std::array<TH2F *, N_PLANES>, 2>, 2> e_pl{};
  std::array<std::array<TH2F *, 2>, 2> e_tot{};
  // run: [cat all=0 / back_fv=1][var tofc=0 / tofcrf=1]
  std::array<std::array<TH2F *, 2>, 2> run{};
  std::vector<TH1 *> owned;

  void book(const std::string &sp, int slot) {
    auto sfx = "_s" + std::to_string(slot);
    for (int c = 0; c < N_CATS; ++c)
      for (int v = 0; v < N_VARS; ++v) {
        const std::string xt = std::string(";") + (v == 0 ? "tof (ns)" : "tof - L/c (ns)") + ";hits";
        for (int p = 0; p < N_PLANES; ++p) {
          for (int b = 0; b < N_PADDLES; ++b) {
            std::string n = std::string(var_names[v]) + "_p" + std::to_string(p) + "_b" + std::to_string(b);
            std::string t = sp + " " + cat_names[c] + " " + var_names[v] + " plane " + plane_names[p] + " bar " +
                            std::to_string(b) + xt;
            bar[c][v][p][b] = new TH1D((n + sfx).c_str(), t.c_str(), NB_T, T_LO[v], T_HI[v]);
            owned.push_back(bar[c][v][p][b]);
          }
          std::string n = std::string(var_names[v]) + "_p" + std::to_string(p) + "_sum";
          std::string t = sp + " " + cat_names[c] + " " + var_names[v] + " plane " + plane_names[p] + xt;
          pl[c][v][p]   = new TH1D((n + sfx).c_str(), t.c_str(), NB_T, T_LO[v], T_HI[v]);
          owned.push_back(pl[c][v][p]);
        }
        std::string n = std::string(var_names[v]) + "_total";
        std::string t = sp + " " + cat_names[c] + " " + var_names[v] + " all planes" + xt;
        tot[c][v]     = new TH1D((n + sfx).c_str(), t.c_str(), NB_T, T_LO[v], T_HI[v]);
        owned.push_back(tot[c][v]);
      }
    for (int c = 0; c < 2; ++c)
      for (int v = 0; v < 2; ++v) {
        const std::string vn = var_names[v + 1];
        const std::string cn = cat_names[c + 1];
        const std::string xt = ";tof - L/c (ns);edep (MeV)";
        for (int p = 0; p < N_PLANES; ++p) {
          if (p != 1 && p != 3)
            continue;
          for (int b = 0; b < N_PADDLES; ++b) {
            std::string n = vn + "_vs_edep_p" + std::to_string(p) + "_b" + std::to_string(b);
            std::string t = sp + " " + cn + " " + vn + " vs edep plane " + plane_names[p] + " bar " + std::to_string(b) + xt;
            e_bar[c][v][p][b] = new TH2F((n + sfx).c_str(), t.c_str(), NB_T2, T2_LO, T2_HI, NB_E, E_LO, E_HI);
            owned.push_back(e_bar[c][v][p][b]);
          }
          std::string n = vn + "_vs_edep_p" + std::to_string(p) + "_sum";
          std::string t = sp + " " + cn + " " + vn + " vs edep plane " + plane_names[p] + xt;
          e_pl[c][v][p] = new TH2F((n + sfx).c_str(), t.c_str(), NB_T2, T2_LO, T2_HI, NB_E, E_LO, E_HI);
          owned.push_back(e_pl[c][v][p]);
        }
        std::string n = vn + "_vs_edep_total";
        e_tot[c][v]   = new TH2F((n + sfx).c_str(), (sp + " " + cn + " " + vn + " vs edep" + xt).c_str(), NB_T2, T2_LO,
                                 T2_HI, NB_E, E_LO, E_HI);
        owned.push_back(e_tot[c][v]);
      }
    for (int c = 0; c < 2; ++c)
      for (int v = 0; v < 2; ++v) {
        const std::string cn = (c == 0) ? "all" : "back_fv";
        const std::string vn = var_names[v + 1];
        std::string n        = vn + "_vs_run_" + cn;
        run[c][v] = new TH2F((n + sfx).c_str(), (sp + " " + cn + " " + vn + " vs run;run;tof - L/c (ns)").c_str(),
                             RUN_HI - RUN_LO, RUN_LO, RUN_HI, NB_T2, T2_LO, T2_HI);
        owned.push_back(run[c][v]);
      }
  }
};

} // namespace lpt

// require_good_vertex: drop events whose reaction point is missing or outside |z| < VTX_ZMAX.
void lad_photon_timing(const char *dat_file, const char *out_file, int nthreads = 4, bool require_good_vertex = false) {
  using namespace lpt;
  using RVd = ROOT::VecOps::RVec<double>;
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0)
    ROOT::EnableImplicitMT(nthreads);

  TChain chain("T");
  {
    std::ifstream fin(dat_file);
    if (!fin.is_open()) {
      std::cerr << "cannot open " << dat_file << "\n";
      return;
    }
    std::string ln;
    while (std::getline(fin, ln)) {
      size_t a = ln.find_first_not_of(" \t\r\n");
      if (a == std::string::npos)
        continue;
      std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
      if (p.empty() || p[0] == '#')
        continue;
      chain.Add(p.c_str());
    }
  }
  std::cout << "[lad_photon_timing] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries())
    return;

  ROOT::RDataFrame df(chain);
  const unsigned nslots = df.GetNSlots();
  const char specs[2]   = {'P', 'H'};
  std::array<std::vector<std::unique_ptr<Hists>>, 2> H;
  for (int is = 0; is < 2; ++is)
    for (unsigned s = 0; s < nslots; ++s) {
      H[is].emplace_back(new Hists);
      H[is].back()->book(std::string(1, specs[is]), s);
    }

  // One event loop per spectrometer reference (different branches).
  for (int is = 0; is < 2; ++is) {
    const std::string sp(1, specs[is]);
    const std::string g = sp + ".ladhod.goodhit_";
    ROOT::RDF::RNode dfs = df;
    dfs = ladtof::define_tof(dfs, sp + "_tof0", specs[is], "0");
    dfs = ladtof::define_tof(dfs, sp + "_tof1", specs[is], "1");
    dfs = ladtof::define_tof(dfs, sp + "_tofrf0", specs[is], "0", true);
    dfs = ladtof::define_tof(dfs, sp + "_tofrf1", specs[is], "1", true);
    std::vector<std::string> cols = {g + "plane_0",        g + "plane_1",         g + "paddle_0",
                                     g + "paddle_1",       sp + "_tof0",          sp + "_tof1",
                                     sp + "_tofrf0",       sp + "_tofrf1",        g + "hit_ypos_0",
                                     g + "hit_ypos_1",     g + "hitedep_MeV_1",   sp + ".react.ok",
                                     sp + ".react.z",      "g.runnum"};
    auto &HS = H[is];
    auto fill = [&HS, require_good_vertex](unsigned slot, const RVd &pl0, const RVd &pl1, const RVd &pd0, const RVd &pd1, const RVd &t0,
                      const RVd &t1, const RVd &r0, const RVd &r1, const RVd &y0, const RVd &y1, const RVd &e1,
                      double vok, double vz, double runnum) {
      if (require_good_vertex && !(vok != 0 && std::fabs(vz) < VTX_ZMAX))
        return;
      Hists &h       = *HS[slot];
      const double z = (vok != 0 && std::fabs(vz) < VTX_ZMAX) ? vz : 0.;
      auto pathlen   = [z](int p, int b, double y) {
        const double x0 = paddle_centre(b), z0 = plane_zpos[p];
        const double c = std::cos(plane_theta[p]), s = std::sin(plane_theta[p]);
        const double X = x0 * c + z0 * s, Z = -x0 * s + z0 * c - z;
        return std::sqrt(X * X + y * y + Z * Z);
      };
      auto valid = [](double p, double pd, double t) {
        return p >= 0 && p < N_PLANES && pd >= 0 && pd < N_PADDLES && std::fabs(t) < SENT;
      };
      const size_t n = pl0.size();
      for (int slotside = 0; slotside < 2; ++slotside) {
        const RVd &PL = slotside ? pl1 : pl0, &PD = slotside ? pd1 : pd0, &T = slotside ? t1 : t0,
                  &R = slotside ? r1 : r0, &Y = slotside ? y1 : y0;
        for (size_t i = 0; i < n; ++i) {
          if (!valid(PL[i], PD[i], T[i]))
            continue;
          const int p = (int)PL[i], b = (int)PD[i];
          const double y   = std::fabs(Y[i]) < SENT ? Y[i] : 0.;
          const double Lc  = pathlen(p, b, y) / C_LIGHT;
          const double val[N_VARS] = {T[i], T[i] - Lc, std::fabs(R[i]) < SENT ? R[i] - Lc : -1e9};

          bool isback = (slotside == 1 && (p == 1 || p == 3));
          bool vetoed = false;
          if (isback) {
            const int pf = p - 1;
            for (size_t k = 0; k < n && !vetoed; ++k) {
              if (pl0[k] != pf || std::fabs(t0[k]) >= SENT)
                continue;
              if (std::abs((int)pd0[k] - b) > FV_DPAD)
                continue;
              if (FV_DT > 0 && std::fabs(t0[k] - T[i]) >= FV_DT)
                continue;
              vetoed = true;
            }
          }
          const bool incat[N_CATS] = {true, isback, isback && !vetoed};
          for (int c = 0; c < N_CATS; ++c) {
            if (!incat[c])
              continue;
            for (int v = 0; v < N_VARS; ++v) {
              h.bar[c][v][p][b]->Fill(val[v]);
              if (!((p == 2 || p == 3) && (b == 1 || b == 9))) { // as lad_tof_fast: 100/101 bars 1,9 out of sums
                h.pl[c][v][p]->Fill(val[v]);
                h.tot[c][v]->Fill(val[v]);
              }
            }
          }
          const bool insum = !((p == 2 || p == 3) && (b == 1 || b == 9));
          if (isback) {
            const double e = std::fabs(e1[i]) < SENT ? e1[i] : -1.;
            for (int c = 0; c < 2; ++c) {
              if (c == 1 && vetoed)
                continue;
              for (int v = 0; v < 2; ++v) {
                h.e_bar[c][v][p][b]->Fill(val[v + 1], e);
                if (insum) {
                  h.e_pl[c][v][p]->Fill(val[v + 1], e);
                  h.e_tot[c][v]->Fill(val[v + 1], e);
                }
              }
            }
          }
          if (insum) {
            for (int v = 0; v < 2; ++v) {
              h.run[0][v]->Fill(runnum, val[v + 1]);
              if (isback && !vetoed)
                h.run[1][v]->Fill(runnum, val[v + 1]);
            }
          }
        }
      }
    };
    std::cout << "[lad_photon_timing] event loop " << sp << "\n";
    dfs.ForeachSlot(fill, cols);
  }

  // Merge slots and write
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "cannot open " << out_file << "\n";
    return;
  }
  TNamed("tof_offsets", Form("P=%.4f,H=%.4f", ladtof::target_offset('P'), ladtof::target_offset('H'))).Write();
  for (int is = 0; is < 2; ++is) {
    auto &HS = H[is];
    for (size_t k = 0; k < HS[0]->owned.size(); ++k)
      for (unsigned s = 1; s < nslots; ++s)
        HS[0]->owned[k]->Add(HS[s]->owned[k]);
    TDirectory *sd = fout.mkdir(std::string(1, specs[is]).c_str());
    auto wr        = [](TDirectory *d, TH1 *h) {
      d->cd();
      std::string n = h->GetName();
      n             = n.substr(0, n.rfind("_s"));
      h->SetName(n.c_str());
      h->Write();
    };
    Hists &h0 = *HS[0];
    for (int c = 0; c < N_CATS; ++c) {
      TDirectory *cd = sd->mkdir(cat_names[c]);
      for (int v = 0; v < N_VARS; ++v) {
        for (int p = 0; p < N_PLANES; ++p) {
          for (int b = 0; b < N_PADDLES; ++b)
            wr(cd, h0.bar[c][v][p][b]);
          wr(cd, h0.pl[c][v][p]);
        }
        wr(cd, h0.tot[c][v]);
      }
      if (c >= 1)
        for (int v = 0; v < 2; ++v) {
          for (int p : {1, 3}) {
            for (int b = 0; b < N_PADDLES; ++b)
              wr(cd, h0.e_bar[c - 1][v][p][b]);
            wr(cd, h0.e_pl[c - 1][v][p]);
          }
          wr(cd, h0.e_tot[c - 1][v]);
        }
    }
    TDirectory *rd = sd->mkdir("run");
    for (int c = 0; c < 2; ++c)
      for (int v = 0; v < 2; ++v)
        wr(rd, h0.run[c][v]);
    for (auto &hs : HS)
      for (auto *o : hs->owned)
        delete o;
  }
  fout.Close();
  std::cout << "[lad_photon_timing] wrote " << out_file << "\n";
}
