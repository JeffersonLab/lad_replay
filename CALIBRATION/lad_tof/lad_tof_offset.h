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

#include <cmath>
#include <iostream>
#include <string>

namespace ladtof {

constexpr double LAD_TOF_OFFSET_SHMS = -1728.0; // ns, lglobal_time_offset_shms
constexpr double LAD_TOF_OFFSET_HMS  = -1752.3; // ns, lglobal_time_offset_hms

inline double target_offset(char spec) { return (spec == 'H' || spec == 'h') ? LAD_TOF_OFFSET_HMS : LAD_TOF_OFFSET_SHMS; }

// Short tag for histogram-cache signatures, so caches made with another ToF
// convention are not reused.
constexpr double VERTEX_ZMAX = 20.0;     // cm, lvertex_zmax
constexpr bool REJECT_BAD_Z  = true;     // lvertex_bad_z_reject = 1: no ToF if the vertex is outside VERTEX_ZMAX
constexpr double KBIG        = 1e38;     // replay value of an uncomputed ToF
constexpr double RF_PERIOD   = 4.00801;  // ns, l_rf_period
constexpr double C_CM_NS     = 29.9792458;

inline std::string signature() {
  return "tofoff=" + std::to_string(LAD_TOF_OFFSET_SHMS) + "," + std::to_string(LAD_TOF_OFFSET_HMS) +
         ";tofz=" + std::to_string(VERTEX_ZMAX) + (REJECT_BAD_Z ? "rej" : "c0");
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

} // namespace ladtof

#endif
