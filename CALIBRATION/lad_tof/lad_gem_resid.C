// lad_gem_resid.C
// -------------------------------------------------------------------------
// GEM residual / alignment study and a window-based GEM layer efficiency, for
// the P spectrometer. Independent of the tracking fit: it predicts where a
// hodoscope proton should cross each GEM module and compares that prediction
// with EVERY GEM cluster in the event (not just the one the tracking picked,
// and not only the 3 highest-ADC clusters the 1D tracking considers).
//
// Prediction: straight line from the event vertex (P.react.x/y/z) to the
// back-plane hodoscope hit (plane 001 or 101, THcLADHodoscope::GetHitPositionLab
// geometry), crossed with each GEM module plane. Residual = cluster pos - the
// predicted coordinate along that module's U or V axis. The module frame
// (origin, U/V axes, normal) is recovered from the data: clust.lab = origin +
// clust.pos * axisHat, so a linear fit of lab vs pos over the first events
// gives both.
//
// Proton selection = lad_hodo_eff.C / lad_tracking_eff.C: good hits with
// isProton_1 == 1 on the back planes 001/101, event has a vertex in the target
// window (react.ok != 0, |react.z| < VTX_ZMAX).
//
// Histograms (per hodo plane 001/101, GEM layer 0/1, axis U/V):
//   resid   : residual (all clusters) vs corrected tof
//   rvp     : residual vs predicted coordinate, tof peak [12,32] only (calibrated ToF, lad_tof_offset.h)
//   pair    : layer-0 residual w.r.t. the line vertex -> layer-1 cluster
//             (layer-1 cluster within 8 cm of its hodo prediction; the unmeasured
//             layer-1 coordinate is taken from the hodo prediction). The lever arm
//             is short and the hodo drops out, so this tests layer 0 against
//             layer 1 + vertex with ~mm resolution.
//   ptof    : proton-cut corrected tof (one entry per proton hit)
//
// Plots / printout (draw step, single-threaded):
//   * residual in the tof peak with narrow gaus + broad gaus + line: narrow mean
//     = alignment offset along that axis, narrow sigma = prediction + GEM
//     resolution (layer-0 U also shows a broad ~8 cm correlated component).
//   * residual vs predicted coordinate profile + pol1: slope != 0 means a
//     rotation / scale error in that axis.
//   * pair residual with the same fit.
//   * window efficiency: for half-widths k*sigma (narrow residual width), count
//     clusters with |r - mean| < k*sigma, subtract the cluster density in local
//     residual sidebands 4-7 sigma (same events, so the proton's own random
//     clusters are removed), then extract the proton peak from the tof spectrum
//     with the lad_tracking_eff fits (flat + trapezoid + gaus for the numerator,
//     flat + gaus for the proton count). eff = gaus area ratio. Each proton has at
//     most one real cluster per (layer, axis), so this is a hit efficiency per
//     layer/axis that does not depend on the chi2 or on accidental matches.
//
// Usage (same signature as the other lad_tof macros, so it runs in the
// slurm/ split-merge workflow):
//   root -l -b -q 'lad_gem_resid.C+("input.dat","out.root",8)'
//   root -l -b -q 'lad_gem_resid.C+("input.dat","out.root",8,"cache.root")'
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>
#include <chrono>
#include <cmath>

#include "lad_tof_offset.h" // calibrated LAD ToF convention (photon peak at tof-L/c = 0)
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace gemres {
// tof axis (corrected tof-L/c), 2.5 ns bins for the per-cluster histograms,
// 0.5 ns for the proton count (as lad_tracking_eff).
const int TOF_NB = 130, PTOF_NB = 650;
const double TOF_LO = ladtof::TCORR_LO, TOF_HI = ladtof::TCORR_HI; // [-168, 157]
const int RES_NB = 400;
const double RES_LO = -20., RES_HI = 20.; // cm
const int PAIR_NB = 400;
const double PAIR_LO = -10., PAIR_HI = 10.; // cm
const int RVP_NB = 200, PRED_NB = 65;
const double RVP_LO = -10., RVP_HI = 10., PRED_LO = -65., PRED_HI = 65.;
const double PAIR_L1_WIN = 8.0; // cm, layer-1 cluster must be this close to its hodo prediction
const double PEAK_LO = ladtof::PEAK_LO, PEAK_HI = ladtof::PEAK_HI; // [12, 32]
// Event-vertex window: react.ok != 0 AND |react.z| < VTX_ZMAX. The target foils
// sit at z = -10, 0, +10 cm; reconstructed vertices far outside (some at |z| of
// metres) give meaningless vertex -> hodo lines for the GEM tracking.
const double VTX_ZMAX = 20.0; // cm
// Efficiency windows |r - mean| < k * sigma (sigma = narrow residual width), with
// the background under the window estimated from local sidebands
// SIDE_LO*sigma < |r - mean| < SIDE_HI*sigma (linear background -> the two-sided
// average is exact for a straight line and close for the smooth broad bump).
const double SIDE_LO = 4., SIDE_HI = 7.;
const std::vector<double> WINDOWS = {1., 1.5, 2., 2.5, 3.};

// Hodoscope geometry (PARAM/LAD/HODO/lhodo_geom.param), by plane index 0..3.
const double HODO_THETA[4] = {2.612452591, 2.615772334, 2.212325822, 2.213860777};
const double HODO_ZPOS[4] = {618.8125006, 658.5516251, 526.1277344, 568.079082};
const char *const PLANE_NAME[2] = {"001", "101"};
const int PLANE_IDX[2] = {1, 3};
const double HODO_R[5] = {615., 655.6, 523., 563.6, 615.}; // tof path radius, as lad_hodo_eff.C
const char *const AXN[2] = {"U", "V"};

struct Frame {
  TVector3 O, A[2], n; // origin, U/V axis, plane normal (pointing away from the target)
  bool ok = false;
};

// THcLADHodoscope::GetHitPositionLab: (centre, ypos, zpos) rotated about y by theta.
inline TVector3 hodo_lab(int pi, double paddle, double ypos) {
  TVector3 h(110. - 22. * paddle, ypos, HODO_ZPOS[pi]);
  h.RotateY(HODO_THETA[pi]);
  return h;
}
inline double corr_tof(int pi, double paddle, double ypos, double tof) {
  const double dx = 110. - 22. * paddle, p2d = std::sqrt(ypos * ypos + dx * dx);
  return tof - std::sqrt(p2d * p2d + HODO_R[pi] * HODO_R[pi]) / 100. / 0.3;
}
// Track through v with direction d crossing the plane of frame f.
inline bool cross(const Frame &f, const TVector3 &v, const TVector3 &d, TVector3 &P) {
  const double nd = f.n.Dot(d);
  if (std::fabs(nd) < 1e-9)
    return false;
  P = v + (f.n.Dot(f.O - v) / nd) * d;
  return true;
}

// Fit functions shared with lad_tracking_eff.C (trapezoid corners -93,-43,32,107).
double trapgaus(double *xx, double *p) {
  const double x = xx[0], amp = p[1] - p[0];
  double trap = 0.;
  if (x >= ladtof::TRAP_C0 && x < ladtof::TRAP_C1)
    trap = amp * (x - ladtof::TRAP_C0) / (ladtof::TRAP_C1 - ladtof::TRAP_C0);
  else if (x >= ladtof::TRAP_C1 && x < ladtof::TRAP_C2)
    trap = amp;
  else if (x >= ladtof::TRAP_C2 && x < ladtof::TRAP_C3)
    trap = amp * (ladtof::TRAP_C3 - x) / (ladtof::TRAP_C3 - ladtof::TRAP_C2);
  return p[0] + trap + p[2] * std::exp(-0.5 * std::pow((x - p[3]) / p[4], 2));
}
double flatgaus(double *xx, double *p) { return p[0] + p[1] * std::exp(-0.5 * std::pow((xx[0] - p[2]) / p[3], 2)); }
// narrow gaus (real hits) + broad gaus (correlated wide component) + line
double gauslin(double *xx, double *p) {
  const double x = xx[0];
  return p[0] * std::exp(-0.5 * std::pow((x - p[1]) / p[2], 2)) +
         p[5] * std::exp(-0.5 * std::pow((x - p[1]) / p[6], 2)) + p[3] + p[4] * x;
}

// Narrow + broad gaus + line fit of a residual projection; mean/sigma are the
// NARROW component's (the broad one soaks up the wide correlated bump seen on
// layer-0 U, so it can't pull the narrow fit).
void fit_resid(TH1D *h, double &mu, double &sig, double &emu, double &esig, TF1 *&f) {
  const double lo = h->GetXaxis()->GetXmin(), hi = h->GetXaxis()->GetXmax();
  // seed the mean from the maximum of a lightly smoothed copy
  TH1D *hs = (TH1D *)h->Clone("hs_tmp");
  hs->Smooth(2);
  const double x0 = hs->GetBinCenter(hs->GetMaximumBin());
  delete hs;
  const double edge = 0.5 * (h->GetBinContent(1) + h->GetBinContent(h->GetNbinsX()));
  f = new TF1(Form("%s_fit", h->GetName()), gauslin, lo, hi, 7);
  f->SetParNames("n_amp", "mean", "n_sigma", "c0", "c1", "b_amp", "b_sigma");
  const double pk = std::max(1., h->GetBinContent(h->FindBin(x0)) - edge);
  f->SetParameters(0.5 * pk, x0, 0.8, edge, 0., 0.5 * pk, 6.);
  f->SetParLimits(0, 0., 10 * pk + 10);
  f->SetParLimits(1, x0 - 2., x0 + 2.);
  f->SetParLimits(2, 0.05, 3.);
  f->SetParLimits(5, 0., 10 * pk + 10);
  f->SetParLimits(6, 3.5, 30.);
  f->SetNpx(800);
  h->Fit(f, "RQN0B");
  mu = f->GetParameter(1);
  sig = std::fabs(f->GetParameter(2));
  emu = f->GetParError(1);
  esig = f->GetParError(2);
}
} // namespace gemres

// Draw step: fits, canvases and the printed summary, from filled or cached histograms.
void lad_gem_resid_draw(TFile *fin, const char *out_file) {
  using namespace gemres;
  ROOT::DisableImplicitMT(); // fits independent of the thread count
  gStyle->SetOptStat(0);
  TFile fout(out_file, "RECREATE");
  auto get2 = [&](const std::string &n) { return dynamic_cast<TH2D *>(fin->Get(n.c_str())); };
  std::vector<TObject *> keep;

  printf("\n[lad_gem_resid] ===== P residuals: tof peak [%g,%g] ns, all clusters vs vertex->hodo line =====\n", PEAK_LO,
         PEAK_HI);
  printf("  plane layer axis   mean(cm)      sigma(cm)   slope vs pred (mrad)\n");
  double MU[2][2][2] = {{{0}}}, SG[2][2][2] = {{{0}}};
  for (int pl = 0; pl < 2; ++pl) {
    TCanvas *cr = new TCanvas(Form("c_resid_pl%s", PLANE_NAME[pl]), Form("residuals, hodo plane %s", PLANE_NAME[pl]),
                              1600, 1000);
    TCanvas *cp = new TCanvas(Form("c_resid_vs_pred_pl%s", PLANE_NAME[pl]),
                              Form("residual vs predicted position, hodo plane %s", PLANE_NAME[pl]), 1600, 1000);
    cr->Divide(2, 2);
    cp->Divide(2, 2);
    for (int L = 0; L < 2; ++L)
      for (int a = 0; a < 2; ++a) {
        const std::string tag = Form("pl%s_L%d_%s", PLANE_NAME[pl], L, AXN[a]);
        TH2D *h2 = get2("resid_" + tag);
        if (!h2)
          continue;
        const int t0 = h2->GetYaxis()->FindBin(PEAK_LO + 1e-6), t1 = h2->GetYaxis()->FindBin(PEAK_HI - 1e-6);
        TH1D *h = h2->ProjectionX(("rpk_" + tag).c_str(), t0, t1);
        h->SetTitle(Form("hodo %s, GEM layer %d, %s: cluster - prediction (tof peak);residual (cm);clusters",
                         PLANE_NAME[pl], L, AXN[a]));
        double mu, sg, emu, esg;
        TF1 *f;
        fit_resid(h, mu, sg, emu, esg, f);
        MU[pl][L][a] = mu;
        SG[pl][L][a] = sg;
        cr->cd(1 + 2 * L + a);
        h->Draw("hist");
        f->SetLineColor(kRed + 1);
        f->Draw("same");
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.045);
        tx.DrawLatex(0.14, 0.84, Form("mean %.2f #pm %.2f cm", mu, emu));
        tx.DrawLatex(0.14, 0.78, Form("#sigma %.2f #pm %.2f cm", sg, esg));
        keep.insert(keep.end(), {h, f});
        // residual vs predicted coordinate, within +-3 sigma of the peak
        TH2D *hp = get2("rvp_" + tag);
        double slope = NAN, eslope = NAN;
        if (hp) {
          cp->cd(1 + 2 * L + a);
          const int r0 = hp->GetXaxis()->FindBin(mu - 3 * sg), r1 = hp->GetXaxis()->FindBin(mu + 3 * sg);
          TProfile *pr = hp->ProfileY(("rvp_prof_" + tag).c_str(), r0, r1);
          pr->SetTitle(Form("hodo %s, GEM layer %d, %s: residual vs predicted (|r-mean|<3#sigma);predicted %s (cm);mean "
                            "residual (cm)",
                            PLANE_NAME[pl], L, AXN[a], AXN[a]));
          pr->SetMarkerStyle(20);
          pr->SetMarkerSize(0.6);
          pr->SetMinimum(mu - 2.5);
          pr->SetMaximum(mu + 2.5);
          // fit the interior only: near the module edges (U +-31 cm, V +-61 cm) the
          // residual window is truncated and the mean bends toward the inside.
          const double pmax = a == 0 ? 20. : 45.;
          TF1 *fl = new TF1(("rvp_fit_" + tag).c_str(), "pol1", -pmax, pmax);
          if (pr->GetEntries() > 20) {
            pr->Fit(fl, "QRN0");
            slope = fl->GetParameter(1) * 1e3;
            eslope = fl->GetParError(1) * 1e3;
          }
          pr->Draw();
          fl->SetLineColor(kRed + 1);
          fl->Draw("same");
          TLatex tx2;
          tx2.SetNDC();
          tx2.SetTextSize(0.045);
          tx2.DrawLatex(0.14, 0.84, Form("slope %.1f #pm %.1f mrad", slope, eslope));
          keep.insert(keep.end(), {pr, fl});
        }
        printf("  %s    %d     %s   %6.2f +- %4.2f   %5.2f +- %4.2f   %7.1f +- %5.1f\n", PLANE_NAME[pl], L, AXN[a], mu, emu,
               sg, esg, slope, eslope);
      }
    fout.cd();
    cr->Write();
    cp->Write();
  }

  printf("\n[lad_gem_resid] ===== P layer-0 residual w.r.t. vertex -> layer-1 cluster line (tof peak) =====\n");
  printf("  plane axis   mean(cm)      sigma(cm)\n");
  for (int pl = 0; pl < 2; ++pl) {
    TCanvas *c = new TCanvas(Form("c_pair_pl%s", PLANE_NAME[pl]), "layer 0 vs vertex->layer 1", 1600, 600);
    c->Divide(2, 1);
    for (int a = 0; a < 2; ++a) {
      const std::string tag = Form("pl%s_%s", PLANE_NAME[pl], AXN[a]);
      TH2D *h2 = get2("pair_" + tag);
      if (!h2)
        continue;
      const int t0 = h2->GetYaxis()->FindBin(PEAK_LO + 1e-6), t1 = h2->GetYaxis()->FindBin(PEAK_HI - 1e-6);
      TH1D *h = h2->ProjectionX(("ppk_" + tag).c_str(), t0, t1);
      h->SetTitle(Form("hodo %s, %s: layer-0 cluster - (vertex #rightarrow layer-1 cluster) (tof peak);residual "
                       "(cm);pairs",
                       PLANE_NAME[pl], AXN[a]));
      double mu, sg, emu, esg;
      TF1 *f;
      fit_resid(h, mu, sg, emu, esg, f);
      c->cd(a + 1);
      h->Draw("hist");
      f->SetLineColor(kRed + 1);
      f->Draw("same");
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.045);
      tx.DrawLatex(0.14, 0.84, Form("mean %.2f #pm %.2f cm", mu, emu));
      tx.DrawLatex(0.14, 0.78, Form("#sigma %.2f #pm %.2f cm", sg, esg));
      keep.insert(keep.end(), {h, f});
      printf("  %s   %s   %6.2f +- %4.2f   %5.2f +- %4.2f\n", PLANE_NAME[pl], AXN[a], mu, emu, sg, esg);
    }
    fout.cd();
    c->Write();
  }

  // ---- window efficiency ----
  printf("\n[lad_gem_resid] ===== P GEM hit efficiency per layer/axis (accidental- and tof-bg-subtracted) =====\n");
  printf("  window = |r - mean| < k*sigma; background from local sidebands %g-%g sigma\n", SIDE_LO, SIDE_HI);
  for (int pl = 0; pl < 2; ++pl) {
    TH1D *hpt = dynamic_cast<TH1D *>(fin->Get(Form("ptof_pl%s", PLANE_NAME[pl])));
    if (!hpt)
      continue;
    printf("  plane %s      ", PLANE_NAME[pl]);
    for (double k : WINDOWS)
      printf("  %-4.1fsig", k);
    printf("   (sigma)\n");
    TCanvas *ce = new TCanvas(Form("c_eff_pl%s", PLANE_NAME[pl]), "GEM hit efficiency vs window", 1000, 700);
    auto *mg = new TMultiGraph();
    auto *lg = new TLegend(0.6, 0.15, 0.89, 0.4);
    const int col[4] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
    for (int L = 0; L < 2; ++L)
      for (int a = 0; a < 2; ++a) {
        const std::string tag = Form("pl%s_L%d_%s", PLANE_NAME[pl], L, AXN[a]);
        TH2D *h2 = get2("resid_" + tag);
        if (!h2)
          continue;
        const double mu = MU[pl][L][a];
        auto *g = new TGraphErrors();
        printf("    L%d %s      ", L, AXN[a]);
        const double sg = SG[pl][L][a];
        for (double k : WINDOWS) {
          const double w = k * sg, slo = SIDE_LO * sg, shi = SIDE_HI * sg;
          // tof spectrum of (clusters in window) - (sideband density * window width)
          TH1D *hin = new TH1D(Form("win_%s_%g", tag.c_str(), k), "", TOF_NB, TOF_LO, TOF_HI);
          const double scale = w / (shi - slo);
          for (int bt = 1; bt <= h2->GetNbinsY(); ++bt) {
            double sin = 0., sside = 0., vin = 0., vside = 0.;
            for (int br = 1; br <= h2->GetNbinsX(); ++br) {
              const double r = std::fabs(h2->GetXaxis()->GetBinCenter(br) - mu);
              const double c = h2->GetBinContent(br, bt);
              if (r < w) {
                sin += c;
                vin += c;
              } else if (r > slo && r < shi) {
                sside += c;
                vside += c;
              }
            }
            hin->SetBinContent(bt, sin - scale * sside);
            hin->SetBinError(bt, std::sqrt(vin + scale * scale * vside));
          }
          // numerator: flat + trapezoid + gaus (mean fixed at 23 ns, as lad_tracking_eff)
          TF1 f2("f2", trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
          f2.SetParameters(10., 30., 50., ladtof::PEAK_MEAN0, 5.);
          f2.FixParameter(3, ladtof::PEAK_MEAN0);
          f2.SetParLimits(4, 1., 15.);
          hin->Fit(&f2, "RQN0");
          const double s2 = std::fabs(f2.GetParameter(4));
          const double num = f2.GetParameter(2) * s2 * std::sqrt(TMath::TwoPi()) / hin->GetBinWidth(1);
          const double enum_ = f2.GetParError(2) * s2 * std::sqrt(TMath::TwoPi()) / hin->GetBinWidth(1);
          // denominator: proton hits, flat + gaus with the same mean and width
          TF1 f1("f1", flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
          f1.SetParameters(50., 100., ladtof::PEAK_MEAN0, s2);
          f1.FixParameter(2, ladtof::PEAK_MEAN0);
          f1.FixParameter(3, s2);
          hpt->Fit(&f1, "RQN0");
          const double den = f1.GetParameter(1) * s2 * std::sqrt(TMath::TwoPi()) / hpt->GetBinWidth(1);
          const double eden = f1.GetParError(1) * s2 * std::sqrt(TMath::TwoPi()) / hpt->GetBinWidth(1);
          const double eff = den > 0 ? num / den : NAN;
          const double eeff = den > 0 ? eff * std::sqrt(std::pow(enum_ / std::max(num, 1e-9), 2) + std::pow(eden / den, 2)) : NAN;
          g->SetPoint(g->GetN(), k, eff);
          g->SetPointError(g->GetN() - 1, 0., eeff);
          printf("  %6.3f", eff);
          delete hin;
        }
        printf("   (%.2f cm)\n", SG[pl][L][a]);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(col[2 * L + a]);
        g->SetLineColor(col[2 * L + a]);
        mg->Add(g, "PL");
        lg->AddEntry(g, Form("layer %d %s (#sigma_{res} %.2f cm)", L, AXN[a], SG[pl][L][a]), "pl");
      }
    ce->cd();
    gPad->SetGridy();
    mg->SetTitle(Form("P hodo %s: GEM hit efficiency (accidental + tof-bg subtracted);window half-width (units of the "
                      "narrow residual #sigma);efficiency",
                      PLANE_NAME[pl]));
    mg->Draw("A");
    mg->SetMinimum(0.);
    mg->SetMaximum(1.2);
    lg->Draw();
    keep.insert(keep.end(), {mg, lg});
    fout.cd();
    ce->Write();
  }
  fout.Close();
  std::cout << "[lad_gem_resid] wrote " << out_file << "\n";
}

void lad_gem_resid(const char *dat_file = "../files/run-lists/all_C3_runlist_SHMS_13p5.dat",
                   const char *out_file = "files/gem_resid/gem_resid_C3_SHMS_13p5_P.root", int nthreads = 4,
                   const char *cache_file = "") {
  using namespace gemres;
  using RVd = ROOT::VecOps::RVec<double>;
  const auto t_start = std::chrono::steady_clock::now();
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  TChain chain("T");
  std::string datlist;
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
      datlist += p + "\n";
    }
  }

  // Cache (same scheme as the other lad_tof macros; fix_sig.C rewrites ";runlist=").
  std::string sig = "lad_gem_resid;gr_v3;" + ladtof::signature() + ";res=" + std::to_string(RES_NB) + ";pair=" + std::to_string(PAIR_NB) +
                    ";l1win=" + std::to_string(PAIR_L1_WIN) +
                    ";vtxz=" + std::to_string(VTX_ZMAX) +
                    ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));
  const bool cache_on = cache_file && cache_file[0];
  if (cache_on && !gSystem->AccessPathName(cache_file)) {
    std::unique_ptr<TFile> fc(TFile::Open(cache_file, "READ"));
    auto *s = fc ? dynamic_cast<TNamed *>(fc->Get("signature")) : nullptr;
    if (s && sig == s->GetTitle()) {
      std::cout << "[lad_gem_resid] cache HIT -> loading histograms, skipping event loop: " << cache_file << "\n";
      lad_gem_resid_draw(fc.get(), out_file);
      return;
    }
  }
  if (cache_on)
    std::cout << "[lad_gem_resid] cache MISS: " << cache_file << "\n";
  std::cout << "[lad_gem_resid] entries: " << chain.GetEntries() << "\n";

  // ---- module frames from the data: lab = O + pos * A, per (layer, axis) ----
  Frame fr[2];
  {
    TTreeReader r(&chain);
    TTreeReaderArray<double> lay(r, "P.gem.clust.layer"), ax(r, "P.gem.clust.axis"), pos(r, "P.gem.clust.pos"),
        lx(r, "P.gem.clust.labx"), ly(r, "P.gem.clust.laby"), lz(r, "P.gem.clust.labz");
    // per (layer, axis): sums for a least-squares line of each lab component vs pos
    double S[2][2][9] = {{{0}}}; // n, sp, spp, sx, sy, sz, spx, spy, spz
    long nread = 0;
    while (r.Next() && nread < 200000) {
      ++nread;
      for (size_t i = 0; i < pos.GetSize(); ++i) {
        const int L = (int)lay[i], a = (int)ax[i];
        if (L < 0 || L > 1 || a < 0 || a > 1)
          continue;
        double *s = S[L][a];
        const double p = pos[i];
        s[0] += 1;
        s[1] += p;
        s[2] += p * p;
        s[3] += lx[i];
        s[4] += ly[i];
        s[5] += lz[i];
        s[6] += p * lx[i];
        s[7] += p * ly[i];
        s[8] += p * lz[i];
      }
    }
    for (int L = 0; L < 2; ++L) {
      TVector3 O[2];
      bool ok = true;
      for (int a = 0; a < 2; ++a) {
        const double *s = S[L][a];
        const double D = s[0] * s[2] - s[1] * s[1];
        if (s[0] < 10 || std::fabs(D) < 1e-9) {
          ok = false;
          continue;
        }
        TVector3 A((s[0] * s[6] - s[1] * s[3]) / D, (s[0] * s[7] - s[1] * s[4]) / D, (s[0] * s[8] - s[1] * s[5]) / D);
        O[a] = TVector3((s[3] - A.X() * s[1]) / s[0], (s[4] - A.Y() * s[1]) / s[0], (s[5] - A.Z() * s[1]) / s[0]);
        fr[L].A[a] = A;
      }
      if (!ok)
        continue;
      fr[L].O = 0.5 * (O[0] + O[1]);
      fr[L].n = fr[L].A[0].Cross(fr[L].A[1]).Unit();
      if (fr[L].n.Dot(fr[L].O) < 0)
        fr[L].n = -fr[L].n;
      fr[L].ok = true;
      printf("[lad_gem_resid] layer %d frame (from %ld events): O=(%.2f,%.2f,%.2f) |O_U-O_V|=%.3f cm\n", L, nread,
             fr[L].O.X(), fr[L].O.Y(), fr[L].O.Z(), (O[0] - O[1]).Mag());
      for (int a = 0; a < 2; ++a)
        printf("    %s axis=(%.4f,%.4f,%.4f) |A|=%.4f\n", AXN[a], fr[L].A[a].X(), fr[L].A[a].Y(), fr[L].A[a].Z(),
               fr[L].A[a].Mag());
      printf("    normal=(%.4f,%.4f,%.4f)  U.V=%.4f\n", fr[L].n.X(), fr[L].n.Y(), fr[L].n.Z(),
             fr[L].A[0].Dot(fr[L].A[1]));
    }
    if (!fr[0].ok || !fr[1].ok) {
      std::cerr << "[lad_gem_resid] could not determine the GEM module frames\n";
      return;
    }
  }

  // ---- event loop: per-slot histograms, merged at the end ----
  if (nthreads > 0)
    ROOT::EnableImplicitMT(nthreads);
  ROOT::RDataFrame rdf(chain);
  const unsigned nslots = rdf.GetNSlots();
  struct Hists {
    std::unique_ptr<TH2D> res[2][2][2], rvp[2][2][2], pair[2][2];
    std::unique_ptr<TH1D> ptof[2];
  };
  std::vector<Hists> H(nslots);
  for (unsigned s = 0; s < nslots; ++s)
    for (int pl = 0; pl < 2; ++pl) {
      H[s].ptof[pl].reset(new TH1D(Form("ptof_pl%s", PLANE_NAME[pl]),
                                   Form("P proton tof, hodo %s;tof-L/c (ns);proton hits", PLANE_NAME[pl]), PTOF_NB,
                                   TOF_LO, TOF_HI));
      for (int a = 0; a < 2; ++a) {
        H[s].pair[pl][a].reset(new TH2D(Form("pair_pl%s_%s", PLANE_NAME[pl], AXN[a]),
                                        ";layer-0 residual w.r.t. vertex->layer-1 (cm);tof-L/c (ns)", PAIR_NB, PAIR_LO,
                                        PAIR_HI, TOF_NB, TOF_LO, TOF_HI));
        for (int L = 0; L < 2; ++L) {
          H[s].res[pl][L][a].reset(new TH2D(Form("resid_pl%s_L%d_%s", PLANE_NAME[pl], L, AXN[a]),
                                            ";cluster - prediction (cm);tof-L/c (ns)", RES_NB, RES_LO, RES_HI, TOF_NB,
                                            TOF_LO, TOF_HI));
          H[s].rvp[pl][L][a].reset(new TH2D(Form("rvp_pl%s_L%d_%s", PLANE_NAME[pl], L, AXN[a]),
                                            ";cluster - prediction (cm);predicted coordinate (cm)", RVP_NB, RVP_LO,
                                            RVP_HI, PRED_NB, PRED_LO, PRED_HI));
        }
      }
    }

  auto fill = [&](unsigned slot, double ok, double vx, double vy, double vz, const RVd &pl1, const RVd &pd1,
                  const RVd &yp1, const RVd &tf1, const RVd &ip1, const RVd &clay, const RVd &cax, const RVd &cpos) {
    if (ok == 0. || !(std::fabs(vz) < VTX_ZMAX))
      return;
    Hists &h = H[slot];
    const TVector3 v(vx, vy, vz);
    for (size_t i = 0; i < pl1.size(); ++i) {
      if (ip1[i] != 1.)
        continue;
      const int pi = (int)std::round(pl1[i]);
      const int pl = pi == 1 ? 0 : pi == 3 ? 1 : -1;
      if (pl < 0)
        continue;
      const double tof = corr_tof(pi, pd1[i], yp1[i], tf1[i]);
      h.ptof[pl]->Fill(tof);
      const TVector3 d = hodo_lab(pi, pd1[i], yp1[i]) - v;
      TVector3 P[2];
      double pred[2][2];
      bool okL[2];
      for (int L = 0; L < 2; ++L) {
        okL[L] = cross(fr[L], v, d, P[L]);
        for (int a = 0; a < 2; ++a)
          pred[L][a] = okL[L] ? fr[L].A[a].Dot(P[L] - fr[L].O) : 0.;
      }
      const bool peak = tof >= PEAK_LO && tof < PEAK_HI;
      for (size_t c = 0; c < cpos.size(); ++c) {
        const int L = (int)clay[c], a = (int)cax[c];
        if (L < 0 || L > 1 || a < 0 || a > 1 || !okL[L])
          continue;
        const double r = cpos[c] - pred[L][a];
        h.res[pl][L][a]->Fill(r, tof);
        if (peak)
          h.rvp[pl][L][a]->Fill(r, pred[L][a]);
      }
      // layer 0 vs the line vertex -> layer-1 cluster
      if (!okL[0] || !okL[1])
        continue;
      for (int a = 0; a < 2; ++a) {
        const int b = 1 - a;
        for (size_t c1 = 0; c1 < cpos.size(); ++c1) {
          if ((int)clay[c1] != 1 || (int)cax[c1] != a || std::fabs(cpos[c1] - pred[1][a]) > PAIR_L1_WIN)
            continue;
          const TVector3 P1 = fr[1].O + cpos[c1] * fr[1].A[a] + pred[1][b] * fr[1].A[b];
          TVector3 P0;
          if (!cross(fr[0], v, P1 - v, P0))
            continue;
          const double p0 = fr[0].A[a].Dot(P0 - fr[0].O);
          for (size_t c0 = 0; c0 < cpos.size(); ++c0)
            if ((int)clay[c0] == 0 && (int)cax[c0] == a)
              h.pair[pl][a]->Fill(cpos[c0] - p0, tof);
        }
      }
    }
  };
  std::cout << "[lad_gem_resid] running event loop (" << nslots << " slots)...\n";
  ROOT::RDF::RNode dfn = ladtof::define_tof(rdf, "P_tof_1", 'P', "1"); // calibrated ToF, any replay
  dfn.ForeachSlot(fill, {"P.react.ok", "P.react.x", "P.react.y", "P.react.z", "P.ladhod.goodhit_plane_1",
                         "P.ladhod.goodhit_paddle_1", "P.ladhod.goodhit_hit_ypos_1", "P_tof_1",
                         "P.ladhod.goodhit_isProton_1", "P.gem.clust.layer", "P.gem.clust.axis", "P.gem.clust.pos"});

  // merge the slots into slot 0 and write the cache
  Hists &M = H[0];
  for (unsigned s = 1; s < nslots; ++s)
    for (int pl = 0; pl < 2; ++pl) {
      M.ptof[pl]->Add(H[s].ptof[pl].get());
      for (int a = 0; a < 2; ++a) {
        M.pair[pl][a]->Add(H[s].pair[pl][a].get());
        for (int L = 0; L < 2; ++L) {
          M.res[pl][L][a]->Add(H[s].res[pl][L][a].get());
          M.rvp[pl][L][a]->Add(H[s].rvp[pl][L][a].get());
        }
      }
    }
  const std::string hist_file = cache_on ? cache_file : std::string(out_file) + ".hists.root";
  {
    TFile fc(hist_file.c_str(), "RECREATE");
    for (int pl = 0; pl < 2; ++pl) {
      M.ptof[pl]->Write();
      for (int a = 0; a < 2; ++a) {
        M.pair[pl][a]->Write();
        for (int L = 0; L < 2; ++L) {
          M.res[pl][L][a]->Write();
          M.rvp[pl][L][a]->Write();
        }
      }
    }
    TNamed("signature", sig.c_str()).Write();
  }
  std::cout << "[lad_gem_resid] wrote histogram " << (cache_on ? "cache: " : "file: ") << hist_file << "\n";
  {
    std::unique_ptr<TFile> fc(TFile::Open(hist_file.c_str(), "READ"));
    lad_gem_resid_draw(fc.get(), out_file);
  }
  const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  printf("[lad_gem_resid] Total time: %.1f s\n", dt);
}
