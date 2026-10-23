// lad_rf_offset_check.C
//
// Per-run inputs for checking the RF offset table (l_rf_offset in lladkine.param)
// against the LAD photon flash.
//
// THcLADKine snaps the bunch time at the target centre, t_b = t_vertex' - z/c
// (t_vertex' includes the z*cos(theta_e)/c electron-path term), to the RF:
//   t_vertex_RFcorr = t_b - remainder(t_b - RF + rf_offset, RF_PERIOD)
// so for any table the RF-corrected ToF of a hit is
//   tof_RF - L/c = (tof - L/c) + remainder(u + rf_offset, RF_PERIOD),   u = fmod(t_b - RF, RF_PERIOD).
// Storing (tof - L/c, u) per hit therefore lets any RF offset table be applied
// afterwards (lad_rf_offset_check.py) without another pass over the data.
//
// Replays with X.ladkin.z_tof (current LADlib) already have the z*cos(theta_e) term in
// t_vertex and z_tof is the z THcLADKine used (KBIG = vertex outside lvertex_zmax, no ToF).
// For older replays both are emulated from X.react.ok/z and X.kin.scat_ang_rad.
//
// Per run and spectrometer X (P: T.shms.pRF_tdcTime, H: T.hms.hRF_tdcTime), only
// events with a valid X.ladkin.t_vertex (in practice that arm's own trigger):
//   run<N>/<X>_u_old   fmod(t_vertex - RF, T)       as get_RF_offset.C (t_vertex as replayed)
//   run<N>/<X>_u_new   fmod(t_b - RF, T)            as the current THcLADKine
//   run<N>/<X>_bfv     TH2F (tof - L/c, u_new), back-plane hits with a front veto
//   run<N>/<X>_all     TH2F (tof - L/c, u_new), all hits (bars 1,9 of 100/101 excluded)
// tof uses the calibrated per-spectrometer offsets and vertex-z convention of
// lad_tof_offset.h. Output is hadd-able (segments of a run add up).
//
// Usage: root -l -b -q 'lad_rf_offset_check.C+("list.dat","out.root",4)'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TROOT.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "lad_tof_offset.h"

namespace lrc {
const int N_PLANES = 5, N_PADDLES = 11;
const double plane_zpos[N_PLANES]  = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double plane_theta[N_PLANES] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
const double T_RF = ladtof::RF_PERIOD, C = ladtof::C_CM_NS, SENT = ladtof::SENTINEL;
const int FV_DPAD = 1;     // front veto: +- paddles
const double FV_DT = 10.;  // front veto: |dt| (ns)
const int NB_T = 800;  const double T_LO = -40., T_HI = 40.; // tof - L/c (0.1 ns)
const int NB_U = 200;                                         // u (0.02 ns)
const int NB_U1 = 400;                                        // 1D phase (0.01 ns)

struct RunH {
  TH1D *u_old, *u_new;
  TH2F *bfv, *all;
};
struct Slot {
  std::map<int, RunH> runs;
  RunH &get(int run, char sp, int slot) {
    auto it = runs.find(run);
    if (it != runs.end())
      return it->second;
    const std::string s = Form("_%c_%d_s%d", sp, run, slot);
    RunH h;
    h.u_old = new TH1D(("u_old" + s).c_str(), ";fmod(t_{vertex} - RF, T) (ns);events", NB_U1, 0., T_RF);
    h.u_new = new TH1D(("u_new" + s).c_str(), ";fmod(t_{b} - RF, T) (ns);events", NB_U1, 0., T_RF);
    h.bfv   = new TH2F(("bfv" + s).c_str(), ";tof - L/c (ns);u (ns)", NB_T, T_LO, T_HI, NB_U, 0., T_RF);
    h.all   = new TH2F(("all" + s).c_str(), ";tof - L/c (ns);u (ns)", NB_T, T_LO, T_HI, NB_U, 0., T_RF);
    return runs.emplace(run, h).first->second;
  }
};
inline double pathlen(int p, int b, double y, double z) {
  const double x0 = 110. - 22. * b, z0 = plane_zpos[p];
  const double c = std::cos(plane_theta[p]), s = std::sin(plane_theta[p]);
  const double X = x0 * c + z0 * s, Z = -x0 * s + z0 * c - z;
  return std::sqrt(X * X + y * y + Z * Z);
}
inline double fmodp(double x) { double r = std::fmod(x, T_RF); return r < 0 ? r + T_RF : r; }
} // namespace lrc

void lad_rf_offset_check(const char *dat_file, const char *out_file, int nthreads = 4) {
  using namespace lrc;
  using RVd = ROOT::VecOps::RVec<double>;
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0)
    ROOT::EnableImplicitMT(nthreads);

  TChain chain("T");
  {
    std::ifstream fin(dat_file);
    std::string ln;
    while (std::getline(fin, ln)) {
      size_t a = ln.find_first_not_of(" \t\r\n");
      if (a == std::string::npos || ln[a] == '#')
        continue;
      chain.Add(ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1).c_str());
    }
  }
  std::cout << "[lad_rf_offset_check] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries())
    return;
  ROOT::RDataFrame df(chain);
  const unsigned nslots = df.GetNSlots();

  const char specs[2] = {'P', 'H'};
  const char *rfcol[2] = {"T.shms.pRF_tdcTime", "T.hms.hRF_tdcTime"};
  std::vector<std::vector<Slot>> S(2, std::vector<Slot>(nslots));

  for (int is = 0; is < 2; ++is) {
    const char sp = specs[is];
    const std::string X(1, sp), g = X + ".ladhod.goodhit_";
    const double target = ladtof::target_offset(sp);
    auto &SS = S[is];
    // lrc_z: z used by THcLADKine (KBIG = rejected vertex); lrc_tvn: t_vertex with the z*cos(theta_e) term
    ROOT::RDF::RNode d = df;
    if (df.HasColumn(X + ".ladkin.z_tof")) {
      std::cout << "[lad_rf_offset_check] " << X << ": using " << X << ".ladkin.z_tof (t_vertex includes the z term)\n";
      d = d.Define("lrc_z_" + X, [](double z) { return z; }, {X + ".ladkin.z_tof"})
              .Define("lrc_tvn_" + X, [](double tv) { return tv; }, {X + ".ladkin.t_vertex"});
    } else {
      std::cout << "[lad_rf_offset_check] " << X << ": no z_tof branch, emulating the vertex-z correction\n";
      d = d.Define("lrc_z_" + X, [](double ok, double z) { return ladtof::z_tof(ok, z); }, {X + ".react.ok", X + ".react.z"})
              .Define("lrc_tvn_" + X,
                      [](double tv, double z, double th) {
                        const double zc = z < 0.1 * ladtof::KBIG ? z : 0.;
                        return tv + (std::isfinite(th) ? zc * std::cos(th) / C : 0.);
                      },
                      {X + ".ladkin.t_vertex", "lrc_z_" + X, X + ".kin.scat_ang_rad"});
    }
    auto fill = [&SS, sp, target](unsigned slot, double run, double rf, double tv, double tvn, double zt,
                                  const RVd &pl0, const RVd &pl1, const RVd &pd0, const RVd &pd1, const RVd &h0,
                                  const RVd &h1, const RVd &y0, const RVd &y1) {
      if (!(std::fabs(tv) < SENT) || !(std::fabs(rf) < SENT) || rf == 0.)
        return;
      RunH &h = SS[slot].get((int)run, sp, slot);
      const bool rejected = zt >= 0.1 * ladtof::KBIG; // vertex outside lvertex_zmax: no LAD ToF
      const double z = rejected ? 0. : zt;
      const double tb = tvn - z / C; // bunch time at the target centre (THcLADKine t_bunch0)
      h.u_old->Fill(fmodp(tv - rf));
      const double u = fmodp(tb - rf);
      h.u_new->Fill(u);
      if (rejected)
        return;
      const size_t n = pl0.size();
      for (int side = 0; side < 2; ++side) {
        const RVd &PL = side ? pl1 : pl0, &PD = side ? pd1 : pd0, &H = side ? h1 : h0, &Y = side ? y1 : y0;
        for (size_t i = 0; i < n; ++i) {
          if (!(PL[i] >= 0 && PL[i] < N_PLANES && PD[i] >= 0 && PD[i] < N_PADDLES && std::fabs(H[i]) < SENT))
            continue;
          const int p = (int)PL[i], b = (int)PD[i];
          const double y = std::fabs(Y[i]) < SENT ? Y[i] : 0.;
          const double tofc = H[i] - tvn + target - pathlen(p, b, y, z) / C;
          if (!((p == 2 || p == 3) && (b == 1 || b == 9)))
            h.all->Fill(tofc, u);
          if (side == 1 && (p == 1 || p == 3)) {
            bool veto = false;
            for (size_t k = 0; k < n && !veto; ++k)
              veto = pl0[k] == p - 1 && std::fabs(h0[k]) < SENT && std::abs((int)pd0[k] - b) <= FV_DPAD &&
                     std::fabs(h0[k] - H[i]) < FV_DT;
            if (!veto)
              h.bfv->Fill(tofc, u);
          }
        }
      }
    };
    std::cout << "[lad_rf_offset_check] event loop " << X << "\n";
    d.ForeachSlot(fill, {"g.runnum", rfcol[is], X + ".ladkin.t_vertex", "lrc_tvn_" + X, "lrc_z_" + X,
                         g + "plane_0", g + "plane_1", g + "paddle_0", g + "paddle_1",
                         g + "hittime_0", g + "hittime_1", g + "hit_ypos_0", g + "hit_ypos_1"});
  }

  TFile fout(out_file, "RECREATE");
  for (int is = 0; is < 2; ++is) {
    std::map<int, RunH> merged;
    for (auto &sl : S[is])
      for (auto &kv : sl.runs) {
        auto it = merged.find(kv.first);
        if (it == merged.end()) {
          merged.emplace(kv.first, kv.second);
          continue;
        }
        it->second.u_old->Add(kv.second.u_old);
        it->second.u_new->Add(kv.second.u_new);
        it->second.bfv->Add(kv.second.bfv);
        it->second.all->Add(kv.second.all);
      }
    for (auto &kv : merged) {
      const std::string dn = Form("run%d", kv.first);
      TDirectory *d = fout.GetDirectory(dn.c_str());
      if (!d)
        d = fout.mkdir(dn.c_str());
      d->cd();
      const std::string X(1, specs[is]);
      kv.second.u_old->Write((X + "_u_old").c_str());
      kv.second.u_new->Write((X + "_u_new").c_str());
      kv.second.bfv->Write((X + "_bfv").c_str());
      kv.second.all->Write((X + "_all").c_str());
    }
  }
  fout.Close();
  std::cout << "[lad_rf_offset_check] wrote " << out_file << "\n";
}
