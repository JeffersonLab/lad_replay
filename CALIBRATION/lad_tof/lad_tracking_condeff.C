// lad_tracking_condeff.C
// -------------------------------------------------------------------------
// Conditional GEM tracking efficiencies for coincidence protons: given that a
// LAD proton hit has a track of 1D variant A, the probability that it also has a
// track of variant B, P(B | A) = eps_AB / eps_A, for every pair of the nine 1D
// projective variants (x / y / combined, GEM0 / GEM1 / both).
//
// Selection = lad_tracking_eff.C / proton_tof_plot.C (P spectrometer only):
// events with a vertex in the target window (P.react.ok != 0, |P.react.z| <
// VTX_ZMAX), good hits with isProton_1 == 1 on planes 001/101, calibrated
// tof - L/c (lad_tof_offset.h) with the paddle-centre path length.
//
// Filled (0.5 ns bins), for chi2 cuts CUTS (a track of variant v means
// 0 <= goodhit_chiSquare_<v> < cut):
//   P_ptof_tot                  all proton hits
//   P_ptof_chi<c>_<a>_<a>       hits with a track of variant a
//   P_ptof_chi<c>_<a>_<b>       hits with tracks of BOTH a and b (a < b)
//
// Extraction (draw step):
//   * every spectrum is fit with the lad_tracking_eff shapes (all hits: flat +
//     gaussian; tracked: flat + fixed-corner trapezoid + gaussian), gaussian mean
//     fixed, gaussian width fitted on all hits and fixed to it for the tracked
//     spectra, so the ratio of gaussian heights is the ratio of signal counts.
//     Binned likelihood fits (the AND spectra can be sparse).
//   * r_X = (gaussian height of X) / (gaussian height of all hits): fraction of
//     signal hits that pass X. f_X = out-of-time sideband counts of X / of all
//     hits: probability that GEM noise alone passes X.
//   * single variant (eq. 8):  eps_A = (r_A - f_A) / (1 - f_A)
//   * pair, no shared (GEM layer, strip axis) slot -- independent noise:
//       r_AB = eps_AB (1-f_A)(1-f_B) + eps_A f_B (1-f_A) + eps_B f_A (1-f_B) + f_AB
//     (true track on both; true A + noise B; noise A + true B; noise on both,
//     the last measured directly as f_AB), solved for eps_AB.
//   * pair sharing a slot (e.g. 1D_x_GEM0 and 1D_x_GEMboth, both using GEM0 V):
//     the noise is correlated and the formula above does not hold; the AND is
//     treated as a single variant, eps_AB ~ (r_AB - f_AB)/(1 - f_AB). This is
//     exact for A = B (so the diagonal is 1) but approximate otherwise; those
//     cells are hatched in the matrix.
//   * P(B | A) = eps_AB / eps_A. The quoted uncertainty propagates the fit
//     errors of the four gaussian heights as if independent (they are not: the
//     AND hits are a subset of A), so it is only indicative.
//
// Usage (same arguments as the other lad_tof macros; runs in slurm/ with
// -m lad_tracking_condeff):
//   root -l -b -q 'lad_tracking_condeff.C+("input.dat","out.root",4)'
//   root -l -b -q 'lad_tracking_condeff.C+("input.dat","out.root",4,"cache.root","png_dir")'
// Redraw from a cache:
//   root -l -b -q -e '.L lad_tracking_condeff.C+' -e 'lad_tracking_condeff_draw("cache.root","out.root")'
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "lad_tof_offset.h" // calibrated LAD ToF convention (photon peak at tof-L/c = 0)

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace cef {
const int NB = 650;
const double XLO = ladtof::TCORR_LO, XHI = ladtof::TCORR_HI;
const double VTX_ZMAX = 20.0; // cm, as lad_tracking_eff.C
const std::vector<int> CUTS = {50, 100};
const double hodo_radii[5] = {615., 655.6, 523., 563.6, 615.}; // cm, by plane index

// GEM (layer, strip axis) slots each variant's track uses: 0 = GEM0 U (y),
// 1 = GEM0 V (x), 2 = GEM1 U (y), 3 = GEM1 V (x).
struct Var {
  std::string dir, tsuf, label;
  std::vector<int> slots;
};
const std::vector<Var> VARS = {
    {"1D_x_GEM0", "_1D_xz_GEM0", "x GEM0", {1}},           {"1D_x_GEM1", "_1D_xz_GEM1", "x GEM1", {3}},
    {"1D_x_GEMboth", "_1D_xz_GEMboth", "x both", {1, 3}},  {"1D_y_GEM0", "_1D_y_GEM0", "y GEM0", {0}},
    {"1D_y_GEM1", "_1D_y_GEM1", "y GEM1", {2}},            {"1D_y_GEMboth", "_1D_y_GEMboth", "y both", {0, 2}},
    {"1D_GEM0", "_1D_GEM0", "x+y GEM0", {0, 1}},           {"1D_GEM1", "_1D_GEM1", "x+y GEM1", {2, 3}},
    {"1D_GEMboth", "_1D_GEMboth", "x+y both", {0, 1, 2, 3}},
};
const int NV = 9;

std::string hname_tot() { return "P_ptof_tot"; }
std::string hname(int cut, int a, int b) {
  if (a > b)
    std::swap(a, b);
  return "P_ptof_chi" + std::to_string(cut) + "_" + std::to_string(a) + "_" + std::to_string(b);
}
bool shares_slot(int a, int b) {
  for (int s : VARS[a].slots)
    for (int t : VARS[b].slots)
      if (s == t)
        return true;
  return false;
}

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

double range_mean(TH1D *h, double lo, double hi) {
  const int b1 = h->FindBin(lo + 1e-6), b2 = h->FindBin(hi - 1e-6);
  return (b2 >= b1) ? h->Integral(b1, b2) / (b2 - b1 + 1) : 0.;
}
double oot_sum(TH1D *h) {
  return h->Integral(h->FindBin(ladtof::OOT_LO1 + 1e-6), h->FindBin(ladtof::OOT_HI1 - 1e-6)) +
         h->Integral(h->FindBin(ladtof::OOT_LO2 + 1e-6), h->FindBin(ladtof::OOT_HI2 - 1e-6));
}

struct FitRes {
  double h = 0, eh = 0, sigma = 0; // gaussian height (counts/bin) and its error
  TF1 *fn = nullptr;
};
FitRes fit_total(TH1D *h) {
  FitRes r;
  r.fn = new TF1("fit_tot", flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
  const double flat = range_mean(h, ladtof::OOT_LO1, ladtof::OOT_HI1);
  r.fn->SetParameters(flat, std::max(1., h->GetMaximum() - flat), ladtof::PEAK_MEAN0, 4.4);
  r.fn->FixParameter(2, ladtof::PEAK_MEAN0);
  r.fn->SetParLimits(3, 1., 15.);
  h->Fit(r.fn, "RQN0L");
  r.h = r.fn->GetParameter(1);
  r.eh = r.fn->GetParError(1);
  r.sigma = std::fabs(r.fn->GetParameter(3));
  return r;
}
FitRes fit_tracked(TH1D *h, double sigma, const std::string &name) {
  FitRes r;
  r.fn = new TF1(name.c_str(), trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
  const double flat = range_mean(h, ladtof::OOT_LO1, ladtof::OOT_HI1);
  const double plat = range_mean(h, ladtof::TRAP_C1, ladtof::PEAK_LO - 5.);
  r.fn->SetParameters(flat, std::max(flat, plat), std::max(1., h->GetBinContent(h->FindBin(ladtof::PEAK_MEAN0)) - plat),
                      ladtof::PEAK_MEAN0, sigma);
  r.fn->FixParameter(3, ladtof::PEAK_MEAN0);
  r.fn->FixParameter(4, sigma);
  h->Fit(r.fn, "RQN0L");
  r.h = r.fn->GetParameter(2);
  r.eh = r.fn->GetParError(2);
  r.sigma = sigma;
  return r;
}

// eps_AB from the measured fractions (see header).
double eps_pair(bool independent, double rAB, double fA, double fB, double fAB, double eA, double eB) {
  if (!independent)
    return (fAB < 1.) ? (rAB - fAB) / (1. - fAB) : 0.;
  const double den = (1. - fA) * (1. - fB);
  return (den > 0.) ? (rAB - fAB - eA * fB * (1. - fA) - eB * fA * (1. - fB)) / den : 0.;
}

int draw_all(const std::function<TH1D *(const std::string &)> &get, const char *out_file, const char *png_dir) {
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");
  const bool png = png_dir && png_dir[0];
  if (png)
    gSystem->mkdir(png_dir, kTRUE);
  TFile fout(out_file, "RECREATE");
  TH1D *htot = get(hname_tot());
  if (!htot) {
    std::cerr << "[lad_tracking_condeff] " << hname_tot() << " missing\n";
    return 0;
  }
  htot = (TH1D *)htot->Clone("P_ptof_tot_fit");
  FitRes ft = fit_total(htot);
  const double oot_tot = oot_sum(htot);
  {
    TDirectory *d = fout.mkdir("spectra");
    d->cd();
    TCanvas c("P_c_tot", "P all proton hits", 800, 550);
    htot->Draw();
    ft.fn->Draw("same");
    c.Write();
  }
  printf("[lad_tracking_condeff] all hits: gaussian height %.1f +- %.1f, sigma %.2f ns, OOT counts %.0f\n", ft.h, ft.eh,
         ft.sigma, oot_tot);

  int ncanv = 0;
  for (int cut : CUTS) {
    // r, f and errors for every (a, b), a <= b
    double r[NV][NV] = {}, er[NV][NV] = {}, f[NV][NV] = {};
    TDirectory *dsp = fout.GetDirectory("spectra")->mkdir(Form("chi%d", cut));
    for (int a = 0; a < NV; ++a)
      for (int b = a; b < NV; ++b) {
        TH1D *h = get(hname(cut, a, b));
        if (!h) {
          std::cerr << "[lad_tracking_condeff] " << hname(cut, a, b) << " missing\n";
          return ncanv;
        }
        h = (TH1D *)h->Clone((hname(cut, a, b) + "_fit").c_str());
        FitRes fr = fit_tracked(h, ft.sigma, "fit_" + hname(cut, a, b));
        r[a][b] = r[b][a] = (ft.h > 0) ? fr.h / ft.h : 0.;
        er[a][b] = er[b][a] = (ft.h > 0) ? std::hypot(fr.eh / ft.h, fr.h * ft.eh / (ft.h * ft.h)) : 0.;
        f[a][b] = f[b][a] = (oot_tot > 0) ? oot_sum(h) / oot_tot : 0.;
        dsp->cd();
        TCanvas c(Form("c_%s_%s", VARS[a].dir.c_str(), VARS[b].dir.c_str()),
                  Form("P proton hits with %s AND %s, chi2 < %d", VARS[a].dir.c_str(), VARS[b].dir.c_str(), cut), 800,
                  550);
        h->SetTitle(Form("%s;tof - L/c (ns);Counts / 0.5 ns", c.GetTitle()));
        h->Draw();
        fr.fn->Draw("same");
        c.Write();
      }
    double eps[NV];
    for (int a = 0; a < NV; ++a)
      eps[a] = (f[a][a] < 1.) ? (r[a][a] - f[a][a]) / (1. - f[a][a]) : 0.;

    // P(B | A) and an indicative error from the gaussian-height errors
    TH2D hP(Form("P_condeff_chi%d", cut),
            Form("P(row track | column track), #chi^{2} < %d;given a track of;probability of also a track of", cut), NV,
            0, NV, NV, 0, NV);
    TH2D hE(Form("P_condeff_err_chi%d", cut), Form("uncertainty of P(row | column), #chi^{2} < %d", cut), NV, 0, NV,
            NV, 0, NV);
    TH1D hEps(Form("P_eff_chi%d", cut), Form("single-variant efficiency, #chi^{2} < %d", cut), NV, 0, NV);
    TH1D hF(Form("P_f_chi%d", cut), Form("noise fraction f, #chi^{2} < %d", cut), NV, 0, NV);
    printf("\n[lad_tracking_condeff] chi2 < %d\n%-14s %8s %8s\n", cut, "variant", "f", "eps");
    for (int a = 0; a < NV; ++a) {
      for (TH1 *hh : {(TH1 *)&hP, (TH1 *)&hE}) {
        hh->GetXaxis()->SetBinLabel(a + 1, VARS[a].label.c_str());
        hh->GetYaxis()->SetBinLabel(a + 1, VARS[a].label.c_str());
      }
      hEps.GetXaxis()->SetBinLabel(a + 1, VARS[a].dir.c_str());
      hF.GetXaxis()->SetBinLabel(a + 1, VARS[a].dir.c_str());
      hEps.SetBinContent(a + 1, eps[a]);
      hEps.SetBinError(a + 1, er[a][a] / std::max(1e-9, 1. - f[a][a]));
      hF.SetBinContent(a + 1, f[a][a]);
      printf("%-14s %8.3f %8.3f\n", VARS[a].dir.c_str(), f[a][a], eps[a]);
    }
    auto condP = [&](int a, int b, double rab, double ra, double rb) {
      const double fa = f[a][a], fb = f[b][b];
      const double ea = (ra - fa) / (1. - fa), eb = (rb - fb) / (1. - fb);
      const double eab = (a == b) ? ea : eps_pair(!shares_slot(a, b), rab, fa, fb, f[a][b], ea, eb);
      return (ea > 0.) ? eab / ea : 0.;
    };
    printf("P(row | column):\n%-10s", "");
    for (int a = 0; a < NV; ++a)
      printf(" %10s", VARS[a].label.c_str());
    printf("\n");
    for (int b = 0; b < NV; ++b) {     // row: the track asked about
      printf("%-10s", VARS[b].label.c_str());
      for (int a = 0; a < NV; ++a) { // column: the track given
        const double p = condP(a, b, r[a][b], r[a][a], r[b][b]);
        // error: vary each input height by its error, add in quadrature
        double e2 = 0.;
        const std::array<std::array<double, 3>, 3> d = {{{er[a][b], 0, 0}, {0, er[a][a], 0}, {0, 0, er[b][b]}}};
        for (const auto &dd : d) {
          const double pv = condP(a, b, r[a][b] + dd[0], r[a][a] + dd[1], r[b][b] + dd[2]);
          e2 += (pv - p) * (pv - p);
        }
        hP.SetBinContent(a + 1, b + 1, p);
        hE.SetBinContent(a + 1, b + 1, std::sqrt(e2));
        printf(" %4.2f+-%4.2f", p, std::sqrt(e2));
      }
      printf("\n");
    }

    fout.cd();
    hP.Write();
    hE.Write();
    hEps.Write();
    hF.Write();
    for (auto [hh, nm] : {std::pair{&hP, "P_c_condeff"}, std::pair{&hE, "P_c_condeff_err"}}) {
      TCanvas c(Form("%s_chi%d", nm, cut), hh->GetTitle(), 1000, 850);
      c.SetLeftMargin(0.15);
      c.SetBottomMargin(0.13);
      c.SetRightMargin(0.13);
      hh->SetMinimum(0.);
      hh->SetMaximum(hh == &hP ? 1.2 : 0.3);
      hh->SetMarkerSize(1.3);
      hh->GetXaxis()->SetTitleOffset(1.4);
      hh->GetYaxis()->SetTitleOffset(1.9);
      hh->GetXaxis()->SetLabelSize(0.028);
      hh->GetYaxis()->SetLabelSize(0.028);
      hh->Draw("COLZ");
      // cell values as explicit labels (the TEXT option's number format is a global
      // style setting, lost when the saved canvas is opened elsewhere)
      for (int a = 0; a < NV; ++a)
        for (int b = 0; b < NV; ++b) {
          TLatex *tv = new TLatex(a + 0.5, b + (hh == &hP ? 0.58 : 0.5), Form("%.2f", hh->GetBinContent(a + 1, b + 1)));
          tv->SetTextAlign(22);
          tv->SetTextSize(0.026);
          tv->Draw();
          if (hh == &hP && a != b) {
            TLatex *te = new TLatex(a + 0.5, b + 0.28, Form("#pm%.2f", hE.GetBinContent(a + 1, b + 1)));
            te->SetTextAlign(22);
            te->SetTextSize(0.017);
            te->Draw();
          }
        }
      // hatch the pairs that share a GEM slot (correlated noise, approximate)
      for (int a = 0; a < NV; ++a)
        for (int b = 0; b < NV; ++b)
          if (a != b && shares_slot(a, b)) {
            TBox *bx = new TBox(a, b, a + 1, b + 1);
            bx->SetFillStyle(3004);
            bx->SetFillColor(kGray + 2);
            bx->Draw();
          }
      TLatex *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.022);
      tx->DrawLatex(0.15, 0.015, "hatched: the two variants share a GEM plane/strip axis (correlated noise), approximate");
      c.Write();
      if (png)
        c.SaveAs(Form("%s/%s_chi%d.png", png_dir, nm, cut));
      ++ncanv;
    }
  }
  fout.Close();
  std::cout << "[lad_tracking_condeff] wrote " << out_file << "\n";
  return ncanv;
}
} // namespace cef

void lad_tracking_condeff_draw(const char *hist_file, const char *out_file = "condeff.root", const char *png_dir = "") {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  TFile fin(hist_file, "READ");
  if (fin.IsZombie()) {
    std::cerr << "[lad_tracking_condeff] cannot open " << hist_file << "\n";
    return;
  }
  std::vector<std::unique_ptr<TH1D>> keep;
  cef::draw_all(
      [&](const std::string &n) -> TH1D * {
        TH1D *h = dynamic_cast<TH1D *>(fin.Get(n.c_str()));
        if (h)
          keep.emplace_back(h);
        return h;
      },
      out_file, png_dir);
}

void lad_tracking_condeff(const char *dat_file, const char *out_file = "condeff.root", int nthreads = 4,
                          const char *cache_file = "", const char *png_dir = "") {
  using namespace cef;
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0)
    ROOT::EnableImplicitMT(nthreads);

  // Runlist -> TChain (same parsing as the other lad_tof macros; datlist feeds the
  // cache signature, so slurm/fix_sig.C works on this cache too).
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
  if (!chain.GetNtrees()) {
    std::cerr << "[lad_tracking_condeff] empty runlist\n";
    return;
  }
  chain.LoadTree(0);
  for (const auto &v : VARS)
    if (!chain.GetBranch(("P.ladhod.goodhit_chiSquare" + v.tsuf).c_str())) {
      std::cerr << "[lad_tracking_condeff] P.ladhod.goodhit_chiSquare" << v.tsuf
                << " absent: needs a replay with the LADlib 1D-cluster tracking\n";
      return;
    }

  // Cache signature; ";runlist=<hash>" must stay the LAST field (slurm/fix_sig.C rewrites it).
  std::string sig = "lad_tracking_condeff;v2;" + ladtof::signature() + ";tof=" + std::to_string(NB) + "," +
                    std::to_string(XLO) + "," + std::to_string(XHI) + ";cuts=";
  for (int c : CUTS)
    sig += std::to_string(c) + ",";
  sig += ";vtxz=" + std::to_string(VTX_ZMAX);
  sig += ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));

  const bool cache_on = cache_file && cache_file[0];
  if (cache_on && !gSystem->AccessPathName(cache_file)) {
    TFile fc(cache_file, "READ");
    auto *s = dynamic_cast<TNamed *>(fc.Get("signature"));
    const bool hit = s && sig == s->GetTitle();
    std::cout << "[lad_tracking_condeff] cache " << (hit ? "HIT -> loading histograms, skipping event loop" : "MISS")
              << ": " << cache_file << "\n";
    if (hit) {
      fc.Close();
      lad_tracking_condeff_draw(cache_file, out_file, png_dir);
      return;
    }
  } else if (cache_on) {
    std::cout << "[lad_tracking_condeff] cache MISS: " << cache_file << "\n";
  }

  using RVd = ROOT::VecOps::RVec<double>;
  using RVi = ROOT::VecOps::RVec<int>;
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf.Filter([](double ok, double z) { return ok != 0. && std::fabs(z) < VTX_ZMAX; },
                                   {"P.react.ok", "P.react.z"}, "has_vertex_P");
  df = ladtof::define_tof(df, "P_tof_1", 'P', "1"); // calibrated ToF, any replay
  const std::string pfx = "P.ladhod.goodhit_";
  // tof - L/c of each good hit; NaN unless it is a proton on plane 001 or 101
  df = df.Define("P_ptof",
                 [](const RVd &pl1, const RVd &pd1, const RVd &y1, const RVd &t1, const RVd &ip1) {
                   RVd r(pl1.size(), std::nan(""));
                   for (size_t i = 0; i < pl1.size(); ++i) {
                     const int pi = (int)std::round(pl1[i]);
                     if (ip1[i] != 1. || (pi != 1 && pi != 3) || i >= t1.size())
                       continue;
                     const double dx = 110. - 22. * pd1[i], p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                     r[i] = t1[i] - std::sqrt(p2d * p2d + hodo_radii[pi] * hodo_radii[pi]) / 100. / 0.3;
                   }
                   return r;
                 },
                 {pfx + "plane_1", pfx + "paddle_1", pfx + "hit_ypos_1", "P_tof_1", pfx + "isProton_1"});
  std::vector<std::string> chicols;
  for (const auto &v : VARS)
    chicols.push_back(pfx + "chiSquare" + v.tsuf);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> booked;
  // only the selected proton hits (P_ptof is NaN for the others)
  df = df.Define("P_ptof_sel", [](const RVd &t) { return t[ROOT::VecOps::Map(t, [](double x) { return std::isfinite(x); })]; },
                 {"P_ptof"});
  booked.push_back(df.Histo1D(
      {hname_tot().c_str(), "P all proton hits;tof - L/c (ns);Counts", NB, XLO, XHI}, "P_ptof_sel"));
  for (int cut : CUTS) {
    // bit v of the mask: the hit has a track of variant v (0 <= chi2 < cut)
    const std::string mcol = "P_mask_chi" + std::to_string(cut);
    const double c = cut;
    df = df.Define(mcol,
                   [c](const RVd &c0, const RVd &c1, const RVd &c2, const RVd &c3, const RVd &c4, const RVd &c5,
                       const RVd &c6, const RVd &c7, const RVd &c8) {
                     const RVd *cs[9] = {&c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8};
                     RVi m(c0.size(), 0);
                     for (int v = 0; v < 9; ++v)
                       for (size_t i = 0; i < m.size() && i < cs[v]->size(); ++i)
                         if ((*cs[v])[i] >= 0. && (*cs[v])[i] < c)
                           m[i] |= (1 << v);
                     return m;
                   },
                   chicols);
    for (int a = 0; a < NV; ++a)
      for (int b = a; b < NV; ++b) {
        const int want = (1 << a) | (1 << b);
        const std::string col = hname(cut, a, b);
        df = df.Define(col,
                       [want](const RVd &t, const RVi &m) {
                         RVd r;
                         for (size_t i = 0; i < t.size() && i < m.size(); ++i)
                           if (std::isfinite(t[i]) && (m[i] & want) == want)
                             r.push_back(t[i]);
                         return r;
                       },
                       {"P_ptof", mcol});
        booked.push_back(df.Histo1D({col.c_str(),
                                     (std::string("P proton hits with ") + VARS[a].dir + " AND " + VARS[b].dir +
                                      ";tof - L/c (ns);Counts")
                                         .c_str(),
                                     NB, XLO, XHI},
                                    col));
      }
  }
  std::cout << "[lad_tracking_condeff] " << chain.GetNtrees() << " files, " << booked.size()
            << " histograms; running event loop...\n";
  booked.front().GetValue();

  if (cache_on) {
    TFile fc(cache_file, "RECREATE");
    TNamed("signature", sig.c_str()).Write();
    for (auto &h : booked)
      h->Write();
    fc.Close();
    std::cout << "[lad_tracking_condeff] wrote histogram cache: " << cache_file << "\n";
  }
  ROOT::DisableImplicitMT(); // fits single-threaded: results independent of the thread count
  cef::draw_all(
      [&](const std::string &n) -> TH1D * {
        for (auto &h : booked)
          if (n == h->GetName())
            return h.GetPtr();
        return nullptr;
      },
      out_file, png_dir);
}
