// lad_tof_offset.h
//
// Shared LAD time-of-flight convention for the CALIBRATION macros.
//
// The replay computes X.ladhod.goodhit_hit_tof_<s> = hittime - X.ladkin.t_vertex + offset,
// with offset = lglobal_time_offset_shms / _hms (PARAM/LAD/LADKINE/lladkine.param).
// Older replays used one shared value (-1710 ns), which put the photon peak at
// tof - L/c = 18 ns (SHMS) and 42 ns (HMS). The calibrated per-spectrometer
// offsets below put photons at tof = L/c (tof - L/c = 0) for both spectrometers.
//
// define_tof() re-references the replay ToF to these offsets hit by hit
// (tof = hittime - t_vertex + offset), so a macro gives the same result on old
// and new replays, and on chains that mix them. The RF-corrected ToF is shifted
// by the same amount (target offset - offset used in the replay, recovered from
// hit_tof - (hittime - t_vertex)).
//
// Vertex z (LADlib with X.ladkin.z_tof): the electron path length uses the vertex
// z (t_vertex later by z*cos(theta_e)/c, also used for the z/c of the RF-corrected
// ToF). A vertex outside |z| < VERTEX_ZMAX gives no ToF (REJECT_BAD_Z), as in LADlib.
// For older replays (no z_tof branch) define_tof applies exactly the same
// correction from X.react.ok/z, X.kin.scat_ang_rad and X.ladkin.t_vertex(_RFcorr),
// including the RF bucket choice, so old and new replays still agree.
//
// Keep LAD_TOF_OFFSET_* in sync with lladkine.param. Every tof window in the
// macros is written in this convention (photon peak at tof - L/c = 0; the former
// SHMS-convention windows moved by -18 ns).

#ifndef LAD_TOF_OFFSET_H
#define LAD_TOF_OFFSET_H

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TF1.h>
#include <TH1.h>

#include <algorithm>

#include <cmath>
#include <iostream>
#include <string>

namespace ladtof {

// 2026-10-09: reset for the 2026 LAD recalibration (was -1728.0 / -1752.3), then SHMS -0.15 ns to match the
// HMS photon peak in the RF-corrected ToF (see lladkine.param)
constexpr double LAD_TOF_OFFSET_SHMS = -1728.01; // ns, lglobal_time_offset_shms
constexpr double LAD_TOF_OFFSET_HMS  = -1751.93; // ns, lglobal_time_offset_hms

inline double target_offset(char spec) { return (spec == 'H' || spec == 'h') ? LAD_TOF_OFFSET_HMS : LAD_TOF_OFFSET_SHMS; }

// Short tag for histogram-cache signatures, so caches made with another ToF
// convention are not reused.
constexpr double VERTEX_ZMAX = 20.0;     // cm, lvertex_zmax
constexpr bool REJECT_BAD_Z  = true;     // lvertex_bad_z_reject = 1: no ToF if the vertex is outside VERTEX_ZMAX
constexpr double KBIG        = 1e38;     // replay value of an uncomputed ToF
constexpr double SENTINEL_EDEP = 1e9;    // |edep| above this = unset
constexpr double RF_PERIOD   = 4.00801;  // ns, l_rf_period
constexpr double C_CM_NS     = 29.9792458;
// Proton selection of the CALIBRATION/lad_tof macros: the replay's back-plane flag isProton_1 (t_back - t_front and
// back-plane edep, THcLADKine::MakeProtonCut) AND a front-plane deposit above FRONT_EDEP_MIN (= lproton_front_edep_cut,
// the replay's isProton_0, computed here from hitedep_MeV_0 so it does not depend on the replay's parameter). On the
// 2026 calibration (MIP at 14 MeV) it keeps ~96% of the coincidence-peak protons and removes about half of the fast
// particles that leak through the back-plane cut at tof-L/c ~ 2-12 ns.
constexpr double FRONT_EDEP_MIN = 25.0; // MeV

inline std::string signature() {
  return "tofoff=" + std::to_string(LAD_TOF_OFFSET_SHMS) + "," + std::to_string(LAD_TOF_OFFSET_HMS) +
         ";tofz=" + std::to_string(VERTEX_ZMAX) + (REJECT_BAD_Z ? "rej" : "c0") +
         ";front=" + std::to_string(FRONT_EDEP_MIN);
}

// Per good hit: 1 if it is a proton for the macros (see FRONT_EDEP_MIN), else 0.
inline ROOT::VecOps::RVec<double> proton_flag(const ROOT::VecOps::RVec<double> &ip1,
                                              const ROOT::VecOps::RVec<double> &ef0) {
  ROOT::VecOps::RVec<double> r(ip1.size(), 0.);
  for (size_t i = 0; i < ip1.size() && i < ef0.size(); ++i)
    r[i] = (ip1[i] == 1. && ef0[i] > FRONT_EDEP_MIN && ef0[i] < SENTINEL_EDEP) ? 1. : 0.;
  return r;
}

// Vertex z used for the ToF corrections (as THcLADKine::fZTof): the reaction
// point, 0 without one, KBIG (= no ToF) outside VERTEX_ZMAX when REJECT_BAD_Z.
inline double z_tof(double react_ok, double z) {
  if (react_ok == 0 || !std::isfinite(z))
    return 0.;
  if (VERTEX_ZMAX <= 0 || std::fabs(z) < VERTEX_ZMAX)
    return z;
  return REJECT_BAD_Z ? KBIG : 0.;
}

constexpr double SENTINEL = 1e9; // |value| above this = unset (replay uses 1e30 / 1e38)

// Per hit: hittime - t_vertex + target. Unset slots / missing vertex keep the replay value.
inline ROOT::VecOps::RVec<double> retarget(const ROOT::VecOps::RVec<double> &tof,
                                           const ROOT::VecOps::RVec<double> &hittime, double tvertex,
                                           double target) {
  ROOT::VecOps::RVec<double> r(tof);
  if (std::fabs(tvertex) >= SENTINEL)
    return r;
  for (size_t i = 0; i < r.size() && i < hittime.size(); ++i)
    if (std::fabs(tof[i]) < SENTINEL && std::fabs(hittime[i]) < SENTINEL)
      r[i] = hittime[i] - tvertex + target;
  return r;
}

// Per hit: tof_rf + (target - offset used in the replay).
inline ROOT::VecOps::RVec<double> retarget_rf(const ROOT::VecOps::RVec<double> &tof_rf,
                                              const ROOT::VecOps::RVec<double> &tof,
                                              const ROOT::VecOps::RVec<double> &hittime, double tvertex,
                                              double target) {
  ROOT::VecOps::RVec<double> r(tof_rf);
  if (std::fabs(tvertex) >= SENTINEL)
    return r;
  for (size_t i = 0; i < r.size() && i < tof.size() && i < hittime.size(); ++i)
    if (std::fabs(tof_rf[i]) < SENTINEL && std::fabs(tof[i]) < SENTINEL && std::fabs(hittime[i]) < SENTINEL)
      r[i] = tof_rf[i] + target - (tof[i] - (hittime[i] - tvertex));
  return r;
}

// Define column `out` = spectrometer `spec`'s good-hit ToF for slot `side` ("0"/"1"),
// in the calibrated convention. rf=true gives the RF-corrected ToF instead.
// Falls back to the raw replay branch (with a warning) if hittime or t_vertex is absent.
inline ROOT::RDF::RNode define_tof(ROOT::RDF::RNode df, const std::string &out, char spec, const std::string &side,
                                   bool rf = false) {
  using RVd          = ROOT::VecOps::RVec<double>;
  const std::string sp(1, spec);
  const std::string g  = sp + ".ladhod.goodhit_";
  const std::string tc = g + "hit_tof_" + side, hc = g + "hittime_" + side, vc = sp + ".ladkin.t_vertex";
  const std::string rc = g + "hit_tof_rfcorr_" + side, vrc = sp + ".ladkin.t_vertex_RFcorr";
  const std::string okc = sp + ".react.ok", zc = sp + ".react.z", thc = sp + ".kin.scat_ang_rad";
  if (!df.HasColumn(hc) || !df.HasColumn(vc)) {
    std::cerr << "[lad_tof_offset] " << hc << " or " << vc << " absent: using " << (rf ? rc : tc)
              << " as replayed (ToF convention of that replay)\n";
    return df.Alias(out, rf ? rc : tc);
  }
  const double target = target_offset(spec);
  const bool replay_has_z = df.HasColumn(sp + ".ladkin.z_tof");
  const bool can_emulate  = df.HasColumn(okc) && df.HasColumn(zc) && df.HasColumn(thc) && (!rf || df.HasColumn(vrc));
  if (replay_has_z || !can_emulate) {
    if (!replay_has_z)
      std::cerr << "[lad_tof_offset] " << okc << "/" << zc << "/" << thc
                << " absent: vertex-z ToF correction NOT applied for " << out << "\n";
    if (!rf)
      return df.Define(out, [target](const RVd &t, const RVd &h, double v) { return retarget(t, h, v, target); },
                       {tc, hc, vc});
    return df.Define(out,
                     [target](const RVd &trf, const RVd &t, const RVd &h, double v) {
                       return retarget_rf(trf, t, h, v, target);
                     },
                     {rc, tc, hc, vc});
  }
  // Older replay: apply the vertex-z correction of the current THcLADKine.
  if (!rf)
    return df.Define(out,
                     [target](const RVd &t, const RVd &h, double v, double ok, double z, double th) {
                       RVd r(t);
                       if (std::fabs(v) >= SENTINEL)
                         return r;
                       const double zt = z_tof(ok, z);
                       const double dt = std::isfinite(th) ? zt * std::cos(th) / C_CM_NS : 0.;
                       for (size_t i = 0; i < r.size() && i < h.size(); ++i)
                         if (std::fabs(t[i]) < SENTINEL && std::fabs(h[i]) < SENTINEL)
                           r[i] = (zt >= SENTINEL) ? KBIG : h[i] - (v + dt) + target;
                       return r;
                     },
                     {tc, hc, vc, okc, zc, thc});
  return df.Define(out,
                   [target](const RVd &trf, const RVd &h, double v, double vrf, double ok, double z, double th) {
                     RVd r(trf);
                     if (std::fabs(v) >= SENTINEL || std::fabs(vrf) >= SENTINEL)
                       return r;
                     const double zt = z_tof(ok, z);
                     if (zt >= SENTINEL) { // vertex outside VERTEX_ZMAX: no ToF
                       for (size_t i = 0; i < r.size(); ++i)
                         if (std::fabs(trf[i]) < SENTINEL)
                           r[i] = KBIG;
                       return r;
                     }
                     const double dt = std::isfinite(th) ? zt * std::cos(th) / C_CM_NS : 0.;
                     const double tb = v + dt - zt / C_CM_NS; // bunch time at the target centre
                     const double rem_old = v - vrf;            // remainder(t_vertex - RF + rf_offset)
                     const double rem_new = (rem_old == 0.) ? 0. : std::remainder(rem_old + dt - zt / C_CM_NS, RF_PERIOD);
                     const double vrf_new = tb - rem_new;
                     for (size_t i = 0; i < r.size() && i < h.size(); ++i)
                       if (std::fabs(trf[i]) < SENTINEL && std::fabs(h[i]) < SENTINEL)
                         r[i] = h[i] - vrf_new + target - zt / C_CM_NS;
                     return r;
                   },
                   {rc, hc, vc, vrc, okc, zc, thc});
}

// tof - L/c regions (ns) used across the macros, calibrated convention.
// Former SHMS-convention values in brackets.
constexpr double TCORR_LO = -168., TCORR_HI = 157.;    // histogram range   [-150, 175]
constexpr double OOT_LO1 = -168., OOT_HI1 = -118.;     // out-of-time       [-150,-100]
constexpr double OOT_LO2 = 107., OOT_HI2 = 157.;       //                   [ 125, 175]
constexpr double IT_LO1 = -43., IT_HI1 = 12.;          // in-time sideband  [ -25,  30]
constexpr double IT_LO2 = 32., IT_HI2 = 107.;          //                   [  50, 125]
constexpr double PEAK_LO = 12., PEAK_HI = 32.;         // proton peak       [  30,  50]
constexpr double TRAP_C0 = -93., TRAP_C1 = -43.;       // trapezoid corners [ -75, -25,
constexpr double TRAP_C2 = 32., TRAP_C3 = 107.;        //                      50, 125]
constexpr double FIT_LO = -168., FIT_HI = 139.;        // fit range         [-150, 157]
constexpr double PEAK_MEAN0 = 23.;                     // gaus start mean   [41]
constexpr double ITB_LO = -43., ITB_HI = 152.;         // in-time-bg trapezoid integral [-25, 170]
constexpr double RATIO_LO = -68., RATIO_HI = 157.;     // track/total ratio x range     [-50, 175]
constexpr double SHIFT_FROM_OLD_SHMS = -18.;           // new - old SHMS convention (ns)

// In-time accidentals with a track (lad_tracking_eff, lad_hodo_eff, lad_hodo_dist): FLAT at the plateau level over
// [ACC_LO, ACC_HI] = the top of the former trapezoid, [TRAP_C1, TRAP_C2] = [-43, 32]. The model is only used inside
// that range: the tracked spectrum is fit with flat + gaussian there (FLATFIT_LO/HI), and the in-time sideband that
// scales the accidentals into the peak window is [ITF_LO, ITF_HI] = [-43, 0]: the flat plateau, stopping short of
// the pre-peak shoulder at ~0-12 ns (real coincidences, not accidentals) and of the falling edge [IT_LO2, IT_HI2],
// which is no longer used. Scale into the peak = (PEAK_HI - PEAK_LO) / (ITF_HI - ITF_LO) = 20/43.
// (IT_LO1/IT_HI1 = [-43, 12] are kept for the macros still on the trapezoid, e.g. lad_proton_tof_min.C.)
constexpr double ITF_LO = -43., ITF_HI = 0.;
constexpr double ACC_LO = TRAP_C1, ACC_HI = TRAP_C2;
constexpr double FLATFIT_LO = ACC_LO, FLATFIT_HI = ACC_HI; // fit range of the tracked spectrum
// Unit-height accidental shape integrated over [a, b) (0 outside [ACC_LO, ACC_HI]).
inline double acc_shape_int(double a, double b) { return std::max(0., std::min(b, ACC_HI) - std::max(a, ACC_LO)); }
// In-time accidentals in the peak window per in-time accidental in the IT window.
inline double acc_scale_peak() {
  const double it = acc_shape_int(ITF_LO, ITF_HI);
  return (it > 0.) ? acc_shape_int(PEAK_LO, PEAK_HI) / it : 0.;
}

// Start values for the proton-ToF fits, taken from the histogram itself, so one fit setup works for the SHMS- and
// HMS-referenced spectra at any statistics. (The former fixed starts 50/70/100 and the gaussian centre fixed at
// PEAK_MEAN0 only suited SHMS spectra; HMS fits failed.) The gaussian centre is free within the peak window.
// trapgaus parameters: flat, trap_top, gaus_h, gaus_mean, gaus_sigma; flatgaus: flat, gaus_h, gaus_mean, gaus_sigma.
inline double mean_in(TH1 *h, double lo, double hi) {
  if (!h) return 0.;
  const int b1 = h->FindBin(lo + 1e-6), b2 = h->FindBin(hi - 1e-6);
  return (b2 >= b1) ? h->Integral(b1, b2) / (b2 - b1 + 1) : 0.;
}
// Position and height of the highest 5-bin running mean inside the peak window.
inline double peak_pos(TH1 *h, double &height) {
  height = 0.;
  double x = PEAK_MEAN0;
  if (!h) return x;
  for (int b = h->FindBin(PEAK_LO + 1e-6); b <= h->FindBin(PEAK_HI - 1e-6); ++b) {
    double s = 0.;
    for (int k = -2; k <= 2; ++k) s += h->GetBinContent(b + k);
    if (s / 5. > height) { height = s / 5.; x = h->GetBinCenter(b); }
  }
  return x;
}
inline void init_trapgaus(TF1 *f, TH1 *h) {
  const double flat = 0.5 * (mean_in(h, OOT_LO1, OOT_HI1) + mean_in(h, OOT_LO2, OOT_HI2));
  const double top = std::max(flat, mean_in(h, IT_LO1, IT_HI1)); // in-time plateau left of the peak
  double hpk = 0.;
  const double x0 = peak_pos(h, hpk);
  const double gh = std::max(hpk - top, 0.1 * hpk + 1e-3);
  f->SetParameters(flat, top, gh, x0, 5.);
  f->SetParLimits(2, 0., 10. * hpk + 10.);
  f->SetParLimits(3, PEAK_LO, PEAK_HI);
  f->SetParLimits(4, 1.5, 12.);
}
// Tracked spectrum, flat + gaussian over [FLATFIT_LO, FLATFIT_HI]: flatgaus parameters (plateau, gaus_h, gaus_mean,
// gaus_sigma); plateau start = mean of the flat in-time sideband [ITF_LO, ITF_HI], centre free in the peak window.
inline void init_flatgaus_trk(TF1 *f, TH1 *h) {
  const double top = mean_in(h, ITF_LO, ITF_HI);
  double hpk = 0.;
  const double x0 = peak_pos(h, hpk);
  f->SetParameters(top, std::max(hpk - top, 0.1 * hpk + 1e-3), x0, 5.);
  f->SetParLimits(1, 0., 10. * hpk + 10.);
  f->SetParLimits(2, PEAK_LO, PEAK_HI);
  f->SetParLimits(3, 1.5, 12.);
}
// Total (no-track) spectrum: centre and width fixed to the proton+track fit.
inline void init_flatgaus(TF1 *f, TH1 *h, double mean, double sigma) {
  const double flat = 0.5 * (mean_in(h, OOT_LO1, OOT_HI1) + mean_in(h, OOT_LO2, OOT_HI2));
  double hpk = 0.;
  peak_pos(h, hpk);
  f->SetParameters(flat, std::max(hpk - flat, 1e-3), mean, sigma);
  f->FixParameter(2, mean);
  f->FixParameter(3, sigma);
}

// GEM acceptance of a LAD hit: does the straight line from the vertex (0, 0, vz) to the hit cross the GEM module(s)
// inside the active area? Same geometry as written_docs/background_subtraction/plot_gem_acceptance.py (validated
// against data there):
//  * GEMs: SHMS 13.5 deg alignment (PARAM/LAD/GEM/lgem_align_shms13p5.param) as lab origin + U/V strip axes; active
//    area |v| < 61.44 cm (3072 V strips x 0.4 mm), -30.72 < u < 24.96 cm (physical U strips 0-1391). The 17 deg
//    alignment (lgem_align_shms17.param) is NOT included: only valid for 13.5 deg run lists.
//  * LAD planes: lhodo_geom.param theta/zpos; a hit at paddle p (0-based, centre x = 110 - 22 p) and position y along
//    the paddle is at R_y(theta) (x, y, zpos) in the lab (THcLADHodoscope::GetHitPositionLab).
// margin (cm) shrinks the active area on the GEM plane, so hits whose real track can scatter (or whose vertex is off)
// across the edge are dropped too. Plane index as in the replay (0..4 = 000, 001, 100, 101, 200).
namespace gemacc {
constexpr double O[2][3] = {{68.2415, 1.3775, -50.5867}, {83.6552, 1.3082, -61.9537}};
constexpr double U[2][3] = {{0.015691, 0.999781, -0.013835}, {0.018632, 0.99971, -0.015249}};
constexpr double V[2][3] = {{0.609729, 0.001398, 0.792609}, {0.607471, 0.000795, 0.794341}};
constexpr double V_HALF = 61.44, U_LO = -30.72, U_HI = 24.96;
constexpr double THETA[5] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
constexpr double ZPOS[5]  = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};

inline bool crosses(int layer, double vz, int plane, double paddle, double ypos, double margin) {
  if (layer < 0 || layer > 1 || plane < 0 || plane > 4)
    return false;
  const double x = 110. - 22. * paddle, th = THETA[plane], z = ZPOS[plane];
  const double P[3] = {x * std::cos(th) + z * std::sin(th), ypos, -x * std::sin(th) + z * std::cos(th)};
  const double vtx[3] = {0., 0., vz};
  const double *u = U[layer], *v = V[layer], *o = O[layer];
  const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  double d[3], dn = 0., on = 0.;
  for (int k = 0; k < 3; ++k) {
    d[k] = P[k] - vtx[k];
    dn += d[k] * n[k];
    on += (o[k] - vtx[k]) * n[k];
  }
  if (std::fabs(dn) < 1e-12)
    return false;
  const double t = on / dn;
  if (!(t > 0.))
    return false;
  double uu = 0., vv = 0.;
  for (int k = 0; k < 3; ++k) {
    const double r = vtx[k] + t * d[k] - o[k];
    uu += r * u[k];
    vv += r * v[k];
  }
  return uu > U_LO + margin && uu < U_HI - margin && std::fabs(vv) < V_HALF - margin;
}

// mode 0: no requirement; 1: inside GEM0; 2: inside GEM0 and GEM1.
inline bool accept(int mode, double vz, int plane, double paddle, double ypos, double margin) {
  if (mode <= 0)
    return true;
  if (!crosses(0, vz, plane, paddle, ypos, margin))
    return false;
  return mode == 1 || crosses(1, vz, plane, paddle, ypos, margin);
}
} // namespace gemacc

} // namespace ladtof

#endif
