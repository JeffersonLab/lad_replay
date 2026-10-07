// proton_tof_plot.C
// -------------------------------------------------------------------------
// Standalone version of lad_tracking_eff.C's <spec>_c_proton_tof canvas. Reads
// the raw replay ROOT files (a runlist, like the other lad_tof macros), fills
// only the histograms this canvas needs, and draws it. The x, noTrackVertex and
// noTrackVertex_x tracking variants are left out; the standard and 1D variants
// are drawn for every chi-square cut.
//
// Filled per spectrometer (P, H), with the same selection as lad_tracking_eff.C
// (event vertex react.ok != 0 with |react.z| < VTX_ZMAX; isProton_1 == 1 hits on planes 001/101; tof
// corrected for the paddle-centre path length):
//   <spec>_tof_corr_proton                          proton total
//   <spec>_tof_corr_proton_track<variant>_cut<ic>   proton + track in the variant's
//                                                   chiSquare window [lo, hi*scale)
//
// Pads (identical to lad_tracking_eff.C):
//   1: proton total tof, fit flat + gaussian (width fixed to pad 2's)
//   2: proton+track tof, fit flat + fixed-corner trapezoid + gaussian
//   3: (track - fit bg)/(total - fit bg), both rebinned x5, GEM-noise corrected (eq. 8)
//   4: raw track/total, no background subtraction
//
// Output layout matches lad_tracking_eff.C:
//   <spec>/proton_id/chi2cut_<val>/<variant>/<spec>_c_proton_tof
//
// Usage (same argument convention as lad_tracking_eff / lad_hodo_eff / lad_hodo_dist,
// so it also runs in the slurm/ split-merge workflow; see README_subMIT_slurm.md):
//   root -l -b -q 'proton_tof_plot.C+("input.dat","out.root")'
//   root -l -b -q 'proton_tof_plot.C+("input.dat","out.root",8)'                   // 8 threads
//   root -l -b -q 'proton_tof_plot.C+("input.dat","out.root",8,"cache.root")'      // + histogram cache
//   root -l -b -q 'proton_tof_plot.C+("input.dat","out.root",8,"","png_dir")'      // + one PNG per canvas
//
// Redraw only, from any file holding these histograms (this macro's cache or a
// lad_tracking_eff cache -- the histogram names are the same):
//   root -l -b -q -e '.L proton_tof_plot.C+' -e 'proton_tof_draw("cache.root","out.root","png_dir")'
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <array>
#include <cmath>

#include "lad_tof_offset.h" // calibrated LAD ToF convention (photon peak at tof-L/c = 0)
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace ptof {
// Same constants as lad_tracking_eff.C.
const int NBINS_TCORR = 650;
const double XMIN_TCORR = ladtof::TCORR_LO, XMAX_TCORR = ladtof::TCORR_HI; // [-168, 157], photon peak at 0
const double SB_LO1 = ladtof::OOT_LO1, SB_HI1 = ladtof::OOT_HI1, SB_LO2 = ladtof::OOT_LO2, SB_HI2 = ladtof::OOT_HI2;
const double CHI_CUT_2D = 100.0, CHI_CUT_1D = 100.0;
// Event-vertex window: react.ok != 0 AND |react.z| < VTX_ZMAX. The target foils
// sit at z = -10, 0, +10 cm; reconstructed vertices far outside (some at |z| of
// metres) give meaningless vertex -> hodo lines for the GEM tracking.
const double VTX_ZMAX = 20.0; // cm
const double CHI_CUT_BASE = 100.0;
const int N_CUTS = 5;
const std::array<double, N_CUTS> CHI_CUT_SCALES = {0.01, 0.1, 0.5, 1.0, 5.0}; // chi2 cuts 1, 10, 50, 100, 500
const std::array<char, 2> specs = {'P', 'H'};
const double hodo_radii[5] = {615., 655.6, 523., 563.6, 615.}; // cm, by plane index

struct Variant {
  std::string dir, tsuf;
  double chi_lo, chi_hi;
};
// lad_tracking_eff.C's variants minus x / noTrackVertex / noTrackVertex_x.
const std::vector<Variant> VARIANTS = {
    {"standard", "", -1e30, CHI_CUT_2D},
    {"1D_x_GEM0", "_1D_xz_GEM0", 0., CHI_CUT_1D},       {"1D_x_GEM1", "_1D_xz_GEM1", 0., CHI_CUT_1D},
    {"1D_x_GEMboth", "_1D_xz_GEMboth", 0., CHI_CUT_1D}, {"1D_y_GEM0", "_1D_y_GEM0", 0., CHI_CUT_1D},
    {"1D_y_GEM1", "_1D_y_GEM1", 0., CHI_CUT_1D},        {"1D_y_GEMboth", "_1D_y_GEMboth", 0., CHI_CUT_1D},
    {"1D_GEM0", "_1D_GEM0", 0., CHI_CUT_1D},            {"1D_GEM1", "_1D_GEM1", 0., CHI_CUT_1D},
    {"1D_GEMboth", "_1D_GEMboth", 0., CHI_CUT_1D},
};

std::string h_tot_name(const std::string &sp) { return sp + "_tof_corr_proton"; }
std::string h_trk_name(const std::string &sp, const Variant &v, int ic) {
  return sp + "_tof_corr_proton_track" + v.tsuf + "_cut" + std::to_string(ic);
}

double trapgaus(double *xx, double *p) {
  const double x = xx[0];
  const double amp = p[1] - p[0];
  double trap;
  if (x < ladtof::TRAP_C0)
    trap = 0.;
  else if (x < ladtof::TRAP_C1)
    trap = amp * (x - ladtof::TRAP_C0) / (ladtof::TRAP_C1 - ladtof::TRAP_C0);
  else if (x < ladtof::TRAP_C2)
    trap = amp;
  else if (x < ladtof::TRAP_C3)
    trap = amp * (ladtof::TRAP_C3 - x) / (ladtof::TRAP_C3 - ladtof::TRAP_C2);
  else
    trap = 0.;
  return p[0] + trap + p[2] * std::exp(-0.5 * std::pow((x - p[3]) / p[4], 2));
}

double flatgaus(double *xx, double *p) { return p[0] + p[1] * std::exp(-0.5 * std::pow((xx[0] - p[2]) / p[3], 2)); }

// Hatched region bands under the histogram (see lad_tracking_eff.C for why hatches).
void draw_shaded(TH1D *h, const std::vector<std::array<double, 4>> &bands) {
  h->DrawCopy();
  gPad->Update();
  const double y1 = gPad->GetUymin(), y2 = gPad->GetUymax();
  for (const auto &bnd : bands) {
    TBox *b = new TBox(bnd[0], y1, bnd[1], y2);
    b->SetFillColor((int)bnd[2]);
    b->SetFillStyle((int)bnd[3]);
    b->SetLineColor((int)bnd[2]);
    b->SetLineWidth(1);
    b->Draw();
  }
  h->DrawCopy("same");
  gPad->RedrawAxis();
}

void draw_gaus_integral(TH1D *h, double amp, double sigma) {
  const double binw = h->GetXaxis()->GetBinWidth(1);
  const double area = amp * std::fabs(sigma) * std::sqrt(TMath::TwoPi());
  TLatex *tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.035);
  tx->SetTextColor(kBlue + 2);
  tx->DrawLatex(0.14, 0.84, Form("Gaus integral: %.0f", binw > 0. ? area / binw : 0.));
}

void subtract_fit_bg(TH1D *h, TF1 *fbg, double orig_binw) {
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double xlo = h->GetXaxis()->GetBinLowEdge(b), xhi = h->GetXaxis()->GetBinUpEdge(b);
    h->SetBinContent(b, h->GetBinContent(b) - (orig_binw > 0. ? fbg->Integral(xlo, xhi) / orig_binw : 0.));
  }
}

// GEM-noise correction of the pad-3 ratio (written_docs/background_subtraction,
// eq. 8; same as lad_tracking_eff.C). f = B_noise/B is the with-track fraction of
// the OOT sidebands, where every tracked hit is GEM noise. The fit-bg-subtracted
// ratio is (C_t - B_GEM)/S = f + (1-f) eps; each bin is replaced by the true-track
// efficiency eps = (r - f)/(1 - f). Bins left empty by the divide are kept.
double noise_frac(TH1 *hAll, TH1 *hTrk) {
  auto sb = [](TH1 *h) {
    return h->Integral(h->FindBin(SB_LO1 + 1e-6), h->FindBin(SB_HI1 - 1e-6)) +
           h->Integral(h->FindBin(SB_LO2 + 1e-6), h->FindBin(SB_HI2 - 1e-6));
  };
  const double b = sb(hAll);
  return b > 0. ? sb(hTrk) / b : 0.;
}

void noise_correct_ratio(TH1 *h, double f) {
  if (!(f < 1.))
    return;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinContent(b) == 0. && h->GetBinError(b) == 0.)
      continue;
    h->SetBinContent(b, (h->GetBinContent(b) - f) / (1. - f));
    h->SetBinError(b, h->GetBinError(b) / (1. - f));
  }
}

TDirectory *mkdirs(TDirectory *base, const std::vector<std::string> &segs) {
  TDirectory *d = base;
  for (const auto &s : segs) {
    TDirectory *n = d->GetDirectory(s.c_str());
    d = n ? n : d->mkdir(s.c_str());
  }
  return d;
}

TCanvas *make_canvas(const std::string &sp, const Variant &v, int ic, TH1D *h_tot, TH1D *h_trk) {
  const std::string &tu = v.tsuf;
  const std::string cc = "_cut" + std::to_string(ic);
  const std::string cutstr = std::to_string((int)std::lround(CHI_CUT_SCALES[ic] * CHI_CUT_BASE));
  const int kSB = kGray + 2, kPk = kRed, fsSB = 3004, fsPk = 3005;
  const std::vector<std::array<double, 4>> sb_pk = {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                                                    {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                                                    {ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}};

  TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(),
                           (sp + " proton tof [" + v.dir + ", chi2<" + cutstr + "]").c_str(), 1400, 1000);
  c->Divide(2, 2);

  // Pad-2 fit first: its gaussian width is reused (fixed) in the pad-1 fit.
  TH1D *ht2 = (TH1D *)h_trk->Clone((sp + "_proton_track_tof_p2" + tu).c_str());
  TF1 *f2 = new TF1((sp + "_fit_trapgaus" + tu + cc).c_str(), trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
  f2->SetParNames("flat", "trap_top", "gaus_h", "gaus_mean", "gaus_sigma");
  f2->SetParameters(50., 70., 100., ladtof::PEAK_MEAN0, 5.);
  f2->FixParameter(3, ladtof::PEAK_MEAN0);
  f2->SetLineColor(kGreen + 2);
  f2->SetNpx(600);
  ht2->Fit(f2, "RQN0");
  const double sig2 = f2->GetParameter(4);

  c->cd(1);
  TH1D *hp1 = (TH1D *)h_tot->Clone((sp + "_proton_tof_p1" + tu).c_str());
  TF1 *f1 = new TF1((sp + "_fit_flatgaus" + tu + cc).c_str(), flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
  f1->SetParNames("flat", "gaus_h", "gaus_mean", "gaus_sigma");
  f1->SetParameters(50., 100., ladtof::PEAK_MEAN0, sig2);
  f1->FixParameter(2, ladtof::PEAK_MEAN0);
  f1->FixParameter(3, sig2);
  f1->SetLineColor(kGreen + 2);
  f1->SetNpx(600);
  hp1->Fit(f1, "RQN0");
  draw_shaded(hp1, sb_pk);
  f1->Draw("same");
  draw_gaus_integral(hp1, f1->GetParameter(1), f1->GetParameter(3));
  delete hp1;

  c->cd(2);
  draw_shaded(ht2, sb_pk);
  f2->Draw("same");
  draw_gaus_integral(ht2, f2->GetParameter(2), f2->GetParameter(4));
  delete ht2;

  // Pad 3: rebin x5, subtract the non-gaussian part of each fit, divide, then
  // correct for GEM noise (eq. 8).
  c->cd(3);
  const double orig_binw = h_trk->GetXaxis()->GetBinWidth(1);
  TF1 fbg_trk((sp + "_bg_trk" + tu + cc).c_str(), trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
  fbg_trk.SetParameters(f2->GetParameter(0), f2->GetParameter(1), 0., f2->GetParameter(3), f2->GetParameter(4));
  TF1 fbg_tot((sp + "_bg_tot" + tu + cc).c_str(), flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
  fbg_tot.SetParameters(f1->GetParameter(0), 0., f1->GetParameter(2), f1->GetParameter(3));
  TH1D *ht_sb2 = (TH1D *)h_trk->Clone((sp + "_proton_track_ratio_trk_fitbg" + tu).c_str());
  TH1D *hp_sb2 = (TH1D *)h_tot->Clone((sp + "_proton_track_ratio_tot_fitbg" + tu).c_str());
  ht_sb2->Rebin(5);
  hp_sb2->Rebin(5);
  subtract_fit_bg(ht_sb2, &fbg_trk, orig_binw);
  subtract_fit_bg(hp_sb2, &fbg_tot, orig_binw);
  TH1D *hratio = (TH1D *)ht_sb2->Clone((sp + "_proton_track_ratio" + tu).c_str());
  const double f_noise = noise_frac(h_tot, h_trk);
  hratio->SetTitle((sp + " proton (track-fitbg)/(total-fitbg), GEM-noise corrected (f=" + Form("%.3f", f_noise) +
                    ");tof-L/c(ns);ratio")
                       .c_str());
  hratio->Divide(hp_sb2);
  noise_correct_ratio(hratio, f_noise);
  hratio->GetXaxis()->SetRangeUser(ladtof::RATIO_LO, ladtof::RATIO_HI);
  hratio->GetYaxis()->SetRangeUser(0., 3);
  draw_shaded(hratio, {{ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}});
  delete hratio;
  delete ht_sb2;
  delete hp_sb2;

  // Pad 4: raw track/total, full binning, no background subtraction.
  c->cd(4);
  TH1D *hratio_raw = (TH1D *)h_trk->Clone((sp + "_proton_track_ratio_raw" + tu).c_str());
  hratio_raw->SetTitle((sp + " proton track/total (no bg-sub);tof-L/c(ns);ratio").c_str());
  hratio_raw->Divide(h_tot);
  hratio_raw->GetYaxis()->SetRangeUser(0., 1.);
  draw_shaded(hratio_raw, sb_pk);
  delete hratio_raw;
  return c; // f1/f2 are canvas primitives and are intentionally not deleted
}

// Draw every available canvas from histograms looked up by name via 'get'.
int draw_all(const std::function<TH1D *(const std::string &)> &get, const char *out_file, const char *png_dir) {
  // Fit single-threaded: with implicit MT on, ROOT parallelizes the fit sums, and
  // the changed summation order moves the (poorly constrained) fits, so the same
  // histograms would give different plots for different thread counts.
  ROOT::DisableImplicitMT();
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[proton_tof_plot] cannot open output " << out_file << "\n";
    return 0;
  }
  const bool png = png_dir && png_dir[0];
  if (png)
    gSystem->mkdir(png_dir, kTRUE);
  int ndone = 0;
  for (char spc : specs) {
    const std::string sp(1, spc);
    TH1D *h_tot = get(h_tot_name(sp));
    if (!h_tot) {
      std::cerr << "[proton_tof_plot] " << h_tot_name(sp) << " missing; skipping " << sp << "\n";
      continue;
    }
    for (int ic = 0; ic < N_CUTS; ++ic) {
      const std::string cutstr = std::to_string((int)std::lround(CHI_CUT_SCALES[ic] * CHI_CUT_BASE));
      for (const auto &v : VARIANTS) {
        TH1D *h_trk = get(h_trk_name(sp, v, ic));
        if (!h_trk)
          continue; // variant absent in the data
        TCanvas *c = make_canvas(sp, v, ic, h_tot, h_trk);
        mkdirs(&fout, {sp, "proton_id", "chi2cut_" + cutstr, v.dir})->cd();
        c->Write();
        if (png)
          c->SaveAs(Form("%s/%s_c_proton_tof_chi2cut_%s_%s.png", png_dir, sp.c_str(), cutstr.c_str(), v.dir.c_str()));
        delete c;
        ++ndone;
      }
    }
  }
  fout.Close();
  std::cout << "[proton_tof_plot] wrote " << ndone << " canvases to " << out_file << "\n";
  return ndone;
}
} // namespace ptof

// Redraw from any file holding the histograms (this macro's cache, or a
// lad_tracking_eff cache).
void proton_tof_draw(const char *hist_file, const char *out_file = "proton_tof.root", const char *png_dir = "") {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  TFile fin(hist_file, "READ");
  if (fin.IsZombie()) {
    std::cerr << "[proton_tof_plot] cannot open " << hist_file << "\n";
    return;
  }
  std::vector<std::unique_ptr<TH1D>> keep;
  ptof::draw_all(
      [&](const std::string &n) -> TH1D * {
        TH1D *h = dynamic_cast<TH1D *>(fin.Get(n.c_str()));
        if (h)
          keep.emplace_back(h);
        return h;
      },
      out_file, png_dir);
}

void proton_tof_plot(const char *dat_file, const char *out_file = "proton_tof.root", int nthreads = 4,
                     const char *cache_file = "", const char *png_dir = "") {
  using namespace ptof;
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0)
    ROOT::EnableImplicitMT(nthreads);
  else
    ROOT::EnableImplicitMT();

  // Runlist -> TChain (same parsing as the other lad_tof macros; datlist feeds
  // the cache signature, so slurm/fix_sig.C works on this cache too).
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
    std::cerr << "[proton_tof_plot] empty runlist\n";
    return;
  }
  chain.LoadTree(0);
  auto has_branch = [&chain](const std::string &n) { return chain.GetBranch(n.c_str()) != nullptr; };

  // Variants present in the data (chiSquare branch in both spectrometers).
  std::vector<Variant> vars;
  for (const auto &v : VARIANTS) {
    bool ok = true;
    for (char spc : specs)
      ok = ok && has_branch(std::string(1, spc) + ".ladhod.goodhit_chiSquare" + v.tsuf);
    if (ok)
      vars.push_back(v);
    else
      std::cout << "[proton_tof_plot] tracking variant '" << v.dir << "' absent in data; skipping\n";
  }

  // Cache signature; ";runlist=<hash>" must stay the LAST field (slurm/fix_sig.C rewrites it).
  std::string sig = "proton_tof_plot;v3;" + ladtof::signature() + ";tof=" + std::to_string(NBINS_TCORR) + "," + std::to_string(XMIN_TCORR) + "," +
                    std::to_string(XMAX_TCORR) + ";cuts=";
  for (double s : CHI_CUT_SCALES)
    sig += std::to_string(s) + ",";
  sig += ";vars=";
  for (const auto &v : vars)
    sig += v.dir + "|" + std::to_string(v.chi_lo) + "|" + std::to_string(v.chi_hi) + ",";
  sig += ";vtxz=" + std::to_string(VTX_ZMAX);
  sig += ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));

  const bool cache_on = cache_file && cache_file[0];
  if (cache_on && !gSystem->AccessPathName(cache_file)) {
    TFile fc(cache_file, "READ");
    auto *s = dynamic_cast<TNamed *>(fc.Get("signature"));
    const bool hit = s && sig == s->GetTitle();
    std::cout << "[proton_tof_plot] cache " << (hit ? "HIT -> loading histograms, skipping event loop" : "MISS")
              << ": " << cache_file << "\n";
    if (hit) {
      fc.Close();
      proton_tof_draw(cache_file, out_file, png_dir);
      return;
    }
  } else if (cache_on) {
    std::cout << "[proton_tof_plot] cache MISS: " << cache_file << "\n";
  }

  // Event loop: only the columns this canvas needs.
  using RVd = ROOT::VecOps::RVec<double>;
  ROOT::RDataFrame rdf(chain);
  std::vector<ROOT::RDF::RResultPtr<TH1D>> booked;
  for (char spc : specs) {
    const std::string sp(1, spc), pfx = sp + ".ladhod.goodhit_";
    ROOT::RDF::RNode df = rdf;
    // histograms use only events with a vertex in the target window, as in lad_tracking_eff.C
    if (has_branch(sp + ".react.ok") && has_branch(sp + ".react.z"))
      df = df.Filter([](double ok, double z) { return ok != 0. && std::fabs(z) < VTX_ZMAX; },
                     {sp + ".react.ok", sp + ".react.z"}, "has_vertex_" + sp);
    else
      std::cout << "[proton_tof_plot] " << sp << ".react.ok/.react.z absent; vertex requirement NOT applied\n";
    df = ladtof::define_tof(df, sp + "_tof_1", spc, "1"); // calibrated ToF, any replay

    // Corrected tof of proton-flagged hits on planes 001/101, optionally requiring
    // chiSquare in [clo, chi) -- identical to lad_tracking_eff.C's mk_proton.
    auto proton_tof = [](bool req_track, double clo, double chi) {
      return [req_track, clo, chi](const RVd &pl1, const RVd &pd1, const RVd &y1, const RVd &t1, const RVd &ip1,
                                   const RVd &chisq) {
        RVd r;
        for (size_t i = 0; i < pl1.size(); ++i) {
          if (ip1[i] != 1.)
            continue;
          if (req_track && !(chisq[i] >= clo && chisq[i] < chi))
            continue;
          const int pi = (int)std::round(pl1[i]);
          if (pi != 1 && pi != 3)
            continue;
          const double dx = 110. - 22. * pd1[i]; // paddle-centre offset, 0-based paddle index
          const double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
          r.push_back(t1[i] - std::sqrt(p2d * p2d + hodo_radii[pi] * hodo_radii[pi]) / 100. / 0.3);
        }
        return r;
      };
    };
    const std::vector<std::string> base = {pfx + "plane_1", pfx + "paddle_1", pfx + "hit_ypos_1", sp + "_tof_1",
                                           pfx + "isProton_1"};
    auto cols = [&](const std::string &chicol) {
      auto c = base;
      c.push_back(chicol);
      return c;
    };

    const std::string ctot = sp + "_ptof_tot";
    df = df.Define(ctot, proton_tof(false, 0., 0.), cols(pfx + "chiSquare"));
    booked.push_back(df.Histo1D({h_tot_name(sp).c_str(), (sp + " tof corr proton;tof-L/c(ns);Counts").c_str(),
                                 NBINS_TCORR, XMIN_TCORR, XMAX_TCORR},
                                ctot));
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (const auto &v : vars) {
        const std::string col = sp + "_ptof_trk" + v.tsuf + "_cut" + std::to_string(ic);
        df = df.Define(col, proton_tof(true, v.chi_lo, v.chi_hi * CHI_CUT_SCALES[ic]), cols(pfx + "chiSquare" + v.tsuf));
        booked.push_back(df.Histo1D({h_trk_name(sp, v, ic).c_str(),
                                     (sp + " tof corr proton+track [" + v.dir + "];tof-L/c(ns);Counts").c_str(),
                                     NBINS_TCORR, XMIN_TCORR, XMAX_TCORR},
                                    col));
      }
  }
  std::cout << "[proton_tof_plot] " << chain.GetNtrees() << " files, " << booked.size()
            << " histograms; running event loop...\n";
  booked.front().GetValue(); // one event loop fills every booked histogram

  if (cache_on) {
    TFile fc(cache_file, "RECREATE");
    TNamed("signature", sig.c_str()).Write();
    for (auto &h : booked)
      h->Write();
    fc.Close();
    std::cout << "[proton_tof_plot] wrote histogram cache: " << cache_file << "\n";
  }
  draw_all(
      [&](const std::string &n) -> TH1D * {
        for (auto &h : booked)
          if (n == h->GetName())
            return h.GetPtr();
        return nullptr;
      },
      out_file, png_dir);
}
