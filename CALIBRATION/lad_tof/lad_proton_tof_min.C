// lad_proton_tof_min.C
// -------------------------------------------------------------------------
// Minimal, fast reproduction of ONE canvas: P_c_proton_tof for the
// 1D_x_GEM0 tracking variant.
//
// Purpose: a small, fast stand-in for the full macros when working on this
// canvas. It carries the same three changes they do (event-vertex
// requirement, filter/compile speed work, OOT*f background-subtraction
// scaling), transcribed from lad_hodo_eff.C / lad_hodo_dist.C, so it tracks
// their behaviour rather than a frozen earlier state:
//
//   vertex   fill only for events with a reconstructed P vertex
//            (P.react.ok != 0), as one compiled RDataFrame Filter. Skipped
//            with a warning when the branch is absent.
//   speed    the ">=1 proton-cut hit" test is folded into that same Filter,
//            so the per-hit Defines never run on events that cannot
//            contribute (every histogram here is proton-gated). Identical
//            output, less work per event.
//   OOT*f    the track-cut OOT suppression factor oot_scale_f(). NOTE: this
//            canvas subtracts fitted backgrounds (subtract_fit_bg) and never
//            forms sflat_it/sflat_pk, so f changes nothing that is drawn
//            here -- it is computed from the two histograms this macro
//            already builds and PRINTED, so the number the parent macros
//            fold into their other plot families can be checked cheaply.
//
// The canvas is built by exactly the code path the full macros use:
//   pad 1  proton total tof, flat+gaus fit (gaus sigma fixed to the pad-2 fit)
//   pad 2  proton+track tof, flat+trapezoid+gaus fit
//   pad 3  (track - fitbg) / (total - fitbg), rebinned by 5
//   pad 4  raw track/total ratio, no bg subtraction
// Shaded bands: grey = OOT sidebands [-168,-118] u [107,157], red = peak [12,32]
// (calibrated ToF, photon peak at tof-L/c = 0; see lad_tof_offset.h).
//
// Deliberately omitted (everything the canvas does not need): the H
// spectrometer, the other 12 tracking variants, all GEM/cluster-ADC/hodo
// quantity plots, the efficiency funnel, and the histogram cache. That is
// what makes it fast -- 1 or 3 RVec Defines and 2-4 TH1D fills per event
// instead of ~80 Defines and several hundred histograms.
//
// Usage (from CALIBRATION/lad_tof/):
//   root -l -b -q 'lad_proton_tof_min.C+("runlist.dat","out.root")'
//   root -l -b -q 'lad_proton_tof_min.C+("some_replay.root","out.root")'
//   root -l -b -q 'lad_proton_tof_min.C+("runlist.dat","out.root",8,200000)'
// dat_file may be either a run list (one ROOT file path per line, '#'
// comments allowed) or a single .root file. max_entries > 0 caps the event
// loop (this forces single-thread, since RDataFrame's Range needs MT off).
// cut_index picks the chi-square cut: 0,1,2 -> chi2 < 50,100,200; -1 = all.
// The canvas is also written next to out_file as a .png for quick eyeballing.
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <array>
#include <chrono>
#include <cmath>

#include "lad_tof_offset.h" // calibrated LAD ToF convention (photon peak at tof-L/c = 0)
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ===== constants, copied verbatim from lad_hodo_eff.C =====
const int NBINS_TCORR = 650; // 0.5 ns bins over [-168, 157] (photon peak at 0, see lad_tof_offset.h)
const double XMIN_TCORR = ladtof::TCORR_LO, XMAX_TCORR = ladtof::TCORR_HI;

// Two-sided sidebands (tof-L/c ns) used by the proton_tof ratio pads.
const double SB_LO1 = ladtof::OOT_LO1, SB_HI1 = ladtof::OOT_HI1, SB_LO2 = ladtof::OOT_LO2, SB_HI2 = ladtof::OOT_HI2;

const double CHI_CUT_1D = 100.0;
// Event-vertex window: react.ok != 0 AND |react.z| < VTX_ZMAX. The target foils
// sit at z = -10, 0, +10 cm; reconstructed vertices far outside (some at |z| of
// metres) give meaningless vertex -> hodo lines for the GEM tracking.
const double VTX_ZMAX = 20.0; // cm
const double CHI_CUT_BASE = 100.0;
const int N_CUTS = 3;
const std::array<double, N_CUTS> CHI_CUT_SCALES = {0.5, 1.0, 2.0};

// tof-L/c regions: all / out-of-time sideband / in-time sideband / peak.
// Only used by oot_scale_f below; the canvas itself does not project regions.
const std::vector<std::vector<std::array<double, 2>>> GREG_INT = {
    {{-1e9, 1e9}},                        // all
    {{SB_LO1, SB_HI1}, {SB_LO2, SB_HI2}}, // oot: [-168,-118] u [107,157]
    {{ladtof::IT_LO1, ladtof::IT_HI1}, {ladtof::IT_LO2, ladtof::IT_HI2}}, // it: [-43,12] u [32,107]
    {{ladtof::PEAK_LO, ladtof::PEAK_HI}}}; // peak: [12,32]

const double hodo_radii[5] = {615., 655.6, 523., 563.6, 615.}; // cm, by plane index

// The one variant this macro plots.
const char *const VAR_DIR = "1D_x_GEM0";
const char *const VAR_TSUF = "_1D_xz_GEM0";

const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_SHMS_13p5.dat";
const char *DEFAULT_OUT_FILE = "files/proton_tof_min/proton_tof_min_P_1D_x_GEM0.root";

// ===== fit models, copied verbatim from lad_hodo_eff.C =====
// flat + fixed-corner trapezoid (corners -93,-43,32,107) + gaussian.
static double trapgaus(double *xx, double *p) {
  const double x = xx[0];
  const double flat = p[0];
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
  const double g = p[2] * std::exp(-0.5 * std::pow((x - p[3]) / p[4], 2));
  return flat + trap + g;
}
static double flatgaus(double *xx, double *p) {
  const double x = xx[0];
  return p[0] + p[1] * std::exp(-0.5 * std::pow((x - p[2]) / p[3], 2));
}

// Track-cut OOT (accidental/flat) suppression factor, from the proton corrected-tof
// spectrum. Per tof bin b, f(b) = (all(b) - h(b)) / (all(b) - oot), where all = the
// no-track (proton-cut) spectrum (hAll), h = the with-track spectrum (hTrk), and oot
// = the with-track OOT-window mean level. f is the fraction of the flat/accidental
// background that survives the track requirement, so wherever the OOT template is
// subtracted from tof window W the parent macros use OOT * f_W instead of OOT. f_W is
// formed as a ratio of window sums over region 'reg' (GREG_INT index: 2 = IT, 3 = peak):
// f_W = (sum_W all - sum_W h) / (sum_W all - oot * N_W). The result is clamped to
// [0,1] (the track cut can only suppress accidentals) and falls back to 1.0 if the
// denominator is degenerate or the inputs are missing (unmodified subtraction).
static double oot_scale_f(const TH1 *hAll, const TH1 *hTrk, int reg) {
  if (!hAll || !hTrk)
    return 1.0;
  const TAxis *ax = hAll->GetXaxis();
  double sH = 0.;
  int nO = 0;
  for (const auto &iv : GREG_INT[1]) { // OOT window -> with-track flat level
    int b1 = ax->FindBin(iv[0] + 1e-6), b2 = ax->FindBin(iv[1] - 1e-6);
    for (int b = b1; b <= b2; ++b) {
      sH += hTrk->GetBinContent(b);
      ++nO;
    }
  }
  const double oot = (nO > 0) ? sH / nO : 0.;
  // ratio of window sums: f_W = (sum all - sum h) / (sum all - oot * N_W)
  double sumA = 0., sumH = 0.;
  int nW = 0;
  for (const auto &iv : GREG_INT[reg]) {
    int b1 = ax->FindBin(iv[0] + 1e-6), b2 = ax->FindBin(iv[1] - 1e-6);
    for (int b = b1; b <= b2; ++b) {
      sumA += hAll->GetBinContent(b);
      sumH += hTrk->GetBinContent(b);
      ++nW;
    }
  }
  const double den = sumA - oot * nW;
  if (std::fabs(den) < 1e-9)
    return 1.0; // sum all -> oot * N_W over the window: undefined
  double f = (sumA - sumH) / den;
  if (f < 0.)
    f = 0.;
  if (f > 1.)
    f = 1.;
  return f;
}

void lad_proton_tof_min(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE,
                        int nthreads = 100, long long max_entries = 0, int cut_index = 1) {

  const auto t_start = std::chrono::steady_clock::now();
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  // Range() (the max_entries cap) is incompatible with implicit MT.
  const bool use_range = (max_entries > 0);
  if (use_range) {
    std::cout << "[proton_tof_min] max_entries=" << max_entries << " -> single thread\n";
  } else if (nthreads > 0) {
    ROOT::EnableImplicitMT(nthreads);
    std::cout << "[proton_tof_min] implicit MT: " << nthreads << " threads\n";
  } else {
    ROOT::EnableImplicitMT();
    std::cout << "[proton_tof_min] implicit MT: all cores\n";
  }

  // ---------------------------------------------------------------
  // 1. TChain -- accept either a run list or a single .root file.
  // ---------------------------------------------------------------
  TChain chain("T");
  {
    const std::string in(dat_file);
    if (in.size() > 5 && in.compare(in.size() - 5, 5, ".root") == 0) {
      chain.Add(in.c_str());
      std::cout << "[proton_tof_min] single file: " << in << "\n";
    } else {
      std::ifstream fin(dat_file);
      if (!fin.is_open()) {
        std::cerr << "cannot open " << dat_file << "\n";
        return;
      }
      std::string ln;
      int nf = 0;
      while (std::getline(fin, ln)) {
        size_t a = ln.find_first_not_of(" \t\r\n");
        if (a == std::string::npos)
          continue;
        std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
        if (p.empty() || p[0] == '#')
          continue;
        chain.Add(p.c_str());
        ++nf;
      }
      std::cout << "[proton_tof_min] run list: " << nf << " files\n";
    }
  }
  std::cout << "[proton_tof_min] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries()) {
    std::cerr << "empty chain\n";
    return;
  }

  // ---------------------------------------------------------------
  // 2. Branch check
  // ---------------------------------------------------------------
  chain.LoadTree(0);
  auto has_branch = [&chain](const std::string &n) { return chain.GetBranch(n.c_str()) != nullptr; };
  const std::string chib = std::string("P.ladhod.goodhit_chiSquare") + VAR_TSUF;
  if (!has_branch(chib)) {
    std::cerr << "[proton_tof_min] " << chib << " absent -- cannot plot " << VAR_DIR << "\n";
    return;
  }

  const int ic_lo = (cut_index < 0) ? 0 : cut_index;
  const int ic_hi = (cut_index < 0) ? N_CUTS - 1 : cut_index;
  if (ic_lo < 0 || ic_hi >= N_CUTS) {
    std::cerr << "[proton_tof_min] cut_index must be -1, 0, 1 or 2\n";
    return;
  }

  // ---------------------------------------------------------------
  // 3. Columns + booking (P spectrometer, 1D_x_GEM0 only)
  // ---------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;
  if (use_range)
    df = df.Range(0, (unsigned int)max_entries);

  using RVd = ROOT::VecOps::RVec<double>;

  const std::string sp = "P";
  const std::string pfx = sp + ".ladhod.goodhit_";
  df = df.Alias(sp + "_plane_1", pfx + "plane_1");
  df = df.Alias(sp + "_paddle_1", pfx + "paddle_1");
  df = ladtof::define_tof(df, sp + "_tof_1", sp[0], "1"); // calibrated ToF, any replay
  df = df.Alias(sp + "_ypos_1", pfx + "hit_ypos_1");
  df = df.Alias(sp + "_isProton_1", pfx + "isProton_1");
  df = df.Alias(sp + "_chiSquare" + VAR_TSUF, chib);

  // Require an event vertex AND >=1 proton-cut hit, exactly as lad_hodo_dist.C
  // does. Every histogram here is proton-gated (isProton_1==1 on plane 001/101),
  // so an event with no such hit contributes nothing -- requiring one upstream
  // skips the per-hit Defines on those events with identical output. One
  // compiled Filter node (explicit column list, so no JIT'd string expression)
  // is shared by everything downstream.
  const bool has_react = has_branch(sp + ".react.ok") && has_branch(sp + ".react.z");
  if (!has_react)
    std::cout << "[proton_tof_min] " << sp << ".react.ok/.react.z absent; vertex requirement NOT applied\n";
  else
    df = df.Filter(
        [](double ok, double z, const RVd &pl1, const RVd &ip1) {
          if (ok == 0. || !(std::fabs(z) < VTX_ZMAX))
            return false;
          for (size_t i = 0; i < pl1.size(); ++i) {
            int p = (int)std::round(pl1[i]);
            if ((p == 1 || p == 3) && ip1[i] == 1.)
              return true;
          }
          return false;
        },
        {sp + ".react.ok", sp + ".react.z", sp + "_plane_1", sp + "_isProton_1"}, "has_vertex_proton_" + sp);

  // ---- proton-tagged corrected-tof column (planes 001 & 101 combined) ----
  // Verbatim from lad_hodo_eff.C's mk_proton.
  auto mk_proton = [&](const std::string &col, bool req_track, const std::string &chicol, double clo, double chi_hi) {
    df = df.Define(col,
                   [req_track, clo, chi_hi](const RVd &pl1, const RVd &pd1, const RVd &y1, const RVd &t1,
                                            const RVd &ip1, const RVd &chi) {
                     RVd r;
                     for (size_t i = 0; i < pl1.size(); ++i) {
                       if (ip1[i] != 1.)
                         continue;
                       if (req_track && !(chi[i] >= clo && chi[i] < chi_hi))
                         continue;
                       int pi = (int)std::round(pl1[i]);
                       if (pi != 1 && pi != 3)
                         continue;
                       double R = hodo_radii[pi];
                       // Transverse offset of the paddle CENTRE (cm), 0-BASED paddle index.
                       double dx = (110. - 22. * pd1[i]);
                       double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                       r.push_back(t1[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3);
                     }
                     return r;
                   },
                   {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", chicol});
  };

  // The chiSquare column is only read when req_track is true; for the no-track
  // total any existing column of the right type will do, so reuse the variant's.
  mk_proton(sp + "_ptof", false, sp + "_chiSquare" + VAR_TSUF, 0., 0.);
  auto r_tot = df.Histo1D({(sp + "_ptof").c_str(), (sp + " proton tof;tof-L/c(ns);Counts").c_str(), NBINS_TCORR,
                           XMIN_TCORR, XMAX_TCORR},
                          sp + "_ptof");

  std::array<ROOT::RDF::RResultPtr<TH1D>, N_CUTS> r_trk;
  for (int ic = ic_lo; ic <= ic_hi; ++ic) {
    const std::string col = sp + "_pttof" + VAR_TSUF + "_cut" + std::to_string(ic);
    mk_proton(col, true, sp + "_chiSquare" + VAR_TSUF, 0., CHI_CUT_1D * CHI_CUT_SCALES[ic]);
    r_trk[ic] = df.Histo1D(
        {col.c_str(), (sp + " proton+track tof;tof-L/c(ns);Counts").c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, col);
  }

  std::cout << "[proton_tof_min] running event loop...\n";
  TH1D *h_proton_tof = (TH1D *)r_tot.GetPtr()->Clone("h_proton_tof");
  TH1D *h_proton_track_tof[N_CUTS] = {nullptr};
  for (int ic = ic_lo; ic <= ic_hi; ++ic)
    h_proton_track_tof[ic] = (TH1D *)r_trk[ic].GetPtr()->Clone(("h_proton_track_tof_cut" + std::to_string(ic)).c_str());
  {
    const auto t_loop = std::chrono::steady_clock::now();
    std::printf("[proton_tof_min] event loop done: %.1f s\n",
                std::chrono::duration<double>(t_loop - t_start).count());
  }

  // ---------------------------------------------------------------
  // 4. Plotting -- the _c_proton_tof block, verbatim
  // ---------------------------------------------------------------
  // Make sure the output directory exists (mkdir -p on out_file's dirname).
  {
    const std::string of(out_file);
    const size_t slash = of.rfind('/');
    if (slash != std::string::npos && slash > 0)
      gSystem->mkdir(of.substr(0, slash).c_str(), kTRUE);
  }
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "cannot open output " << out_file << "\n";
    return;
  }
  auto mkpath = [](TDirectory *base, const std::string &path) -> TDirectory * {
    TDirectory *dcur = base;
    std::string seg;
    std::istringstream ss(path);
    while (std::getline(ss, seg, '/')) {
      if (seg.empty())
        continue;
      TDirectory *nd = dcur->GetDirectory(seg.c_str());
      if (!nd)
        nd = dcur->mkdir(seg.c_str());
      dcur = nd;
    }
    return dcur ? dcur : base;
  };
  auto draw_shaded = [](TH1D *h, const std::vector<std::array<double, 4>> &bands) {
    h->DrawCopy();
    gPad->Update();
    double y1 = gPad->GetUymin(), y2 = gPad->GetUymax();
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
  };
  auto draw_gaus_integral = [](TH1D *h, double amp, double sigma) {
    const double binw = h->GetXaxis()->GetBinWidth(1);
    const double area = amp * std::fabs(sigma) * std::sqrt(TMath::TwoPi());
    const double nev = (binw > 0.) ? area / binw : 0.;
    TLatex *tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.035);
    tx->SetTextColor(kBlue + 2);
    tx->DrawLatex(0.14, 0.84, Form("Gaus integral: %.0f", nev));
  };
  auto subtract_fit_bg = [](TH1D *h, TF1 *fbg, double orig_binw) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double xlo = h->GetXaxis()->GetBinLowEdge(b), xhi = h->GetXaxis()->GetBinUpEdge(b);
      double bg = (orig_binw > 0.) ? fbg->Integral(xlo, xhi) / orig_binw : 0.;
      h->SetBinContent(b, h->GetBinContent(b) - bg);
    }
  };
  const int kSB = kGray + 2, kPk = kRed, fsSB = 3004, fsPk = 3005;

  TDirectory *sdir = fout.mkdir(sp.c_str());
  for (int ic = ic_lo; ic <= ic_hi; ++ic) {
    const double cutval = CHI_CUT_SCALES[ic] * CHI_CUT_BASE;
    const std::string cutstr = std::to_string((int)std::lround(cutval));
    TDirectory *td = mkpath(sdir, "chi2cut_" + cutstr + "/" + VAR_DIR);
    td->cd();
    const std::string tu = VAR_TSUF;
    const std::string cc = "_cut" + std::to_string(ic);

    // ---- fit the proton+track tof (pad-2 model) ----
    TH1D *hpt = h_proton_track_tof[ic];
    TF1 *f2 = new TF1((sp + "_fit_trapgaus" + tu + cc).c_str(), trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
    f2->SetParNames("flat", "trap_top", "gaus_h", "gaus_mean", "gaus_sigma");
    f2->SetParameters(50., 70., 100., ladtof::PEAK_MEAN0, 5.);
    f2->FixParameter(3, ladtof::PEAK_MEAN0);
    f2->SetLineColor(kGreen + 2);
    f2->SetNpx(600);
    if (hpt)
      hpt->Fit(f2, "RQN0");
    const double sig2 = f2->GetParameter(4);

    TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(),
                             (sp + " proton tof [" + VAR_DIR + ", chi2<" + cutstr + "]").c_str(), 1400, 1000);
    c->Divide(2, 2);
    // pad 1: proton total, gaussian width fixed to the pad-2 width.
    c->cd(1);
    TH1D *hp1 = (TH1D *)h_proton_tof->Clone((sp + "_ptof_p1" + tu + cc).c_str());
    TF1 *f1 = new TF1((sp + "_fit_flatgaus" + tu + cc).c_str(), flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
    f1->SetParNames("flat", "gaus_h", "gaus_mean", "gaus_sigma");
    f1->SetParameters(50., 100., ladtof::PEAK_MEAN0, sig2);
    f1->FixParameter(2, ladtof::PEAK_MEAN0);
    f1->FixParameter(3, sig2);
    f1->SetLineColor(kGreen + 2);
    f1->SetNpx(600);
    hp1->Fit(f1, "RQN0");
    draw_shaded(hp1, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                      {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                      {ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}});
    f1->Draw("same");
    draw_gaus_integral(hp1, f1->GetParameter(1), f1->GetParameter(3));
    delete hp1;
    // pad 2: proton+track.
    c->cd(2);
    TH1D *ht2 = (TH1D *)hpt->Clone((sp + "_pttof_p2" + tu + cc).c_str());
    draw_shaded(ht2, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                      {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                      {ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}});
    f2->Draw("same");
    draw_gaus_integral(ht2, f2->GetParameter(2), f2->GetParameter(4));
    delete ht2;
    // pad 3: (track - fitbg) / (total - fitbg), rebinned by 5.
    c->cd(3);
    const double orig_binw = hpt->GetXaxis()->GetBinWidth(1);
    TF1 *fbg_trk = new TF1((sp + "_bg_trk" + tu + cc).c_str(), trapgaus, ladtof::FIT_LO, ladtof::FIT_HI, 5);
    fbg_trk->SetParameters(f2->GetParameter(0), f2->GetParameter(1), 0., f2->GetParameter(3), f2->GetParameter(4));
    TF1 *fbg_tot = new TF1((sp + "_bg_tot" + tu + cc).c_str(), flatgaus, ladtof::FIT_LO, ladtof::FIT_HI, 4);
    fbg_tot->SetParameters(f1->GetParameter(0), 0., f1->GetParameter(2), f1->GetParameter(3));
    TH1D *ht_rb5 = (TH1D *)hpt->Clone((sp + "_ratio_trk_rb5" + tu + cc).c_str());
    ht_rb5->Rebin(5);
    TH1D *hp_rb5 = (TH1D *)h_proton_tof->Clone((sp + "_ratio_tot_rb5" + tu + cc).c_str());
    hp_rb5->Rebin(5);
    subtract_fit_bg(ht_rb5, fbg_trk, orig_binw);
    subtract_fit_bg(hp_rb5, fbg_tot, orig_binw);
    TH1D *hratio = (TH1D *)ht_rb5->Clone((sp + "_ratio" + tu + cc).c_str());
    hratio->SetTitle((sp + " proton (track-fitbg)/(total-fitbg);tof-L/c(ns);ratio").c_str());
    hratio->Divide(hp_rb5);
    hratio->GetXaxis()->SetRangeUser(ladtof::RATIO_LO, ladtof::RATIO_HI);
    hratio->GetYaxis()->SetRangeUser(0., 3.);
    draw_shaded(hratio, {{ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}});
    delete hratio;
    delete ht_rb5;
    delete hp_rb5;
    delete fbg_trk;
    delete fbg_tot;
    // pad 4: raw track/total ratio (no bg-sub, full binning).
    c->cd(4);
    TH1D *ht_raw = (TH1D *)hpt->Clone((sp + "_ratio_raw_trk" + tu + cc).c_str());
    TH1D *hp_raw = (TH1D *)h_proton_tof->Clone((sp + "_ratio_raw_tot" + tu + cc).c_str());
    TH1D *hratio_raw = (TH1D *)ht_raw->Clone((sp + "_ratio_raw" + tu + cc).c_str());
    hratio_raw->SetTitle((sp + " proton track/total (no bg-sub);tof-L/c(ns);ratio").c_str());
    hratio_raw->Divide(hp_raw);
    hratio_raw->GetYaxis()->SetRangeUser(0., 1.);
    draw_shaded(hratio_raw, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                             {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                             {ladtof::PEAK_LO, ladtof::PEAK_HI, (double)kPk, (double)fsPk}});
    delete hratio_raw;
    delete ht_raw;
    delete hp_raw;

    // Fit parameters, printed so successive runs can be compared numerically
    // without opening the canvas.
    std::printf("[proton_tof_min] chi2<%s  pad2 trapgaus: flat=%.2f trap_top=%.2f gaus_h=%.2f mean=%.2f sigma=%.3f\n",
                cutstr.c_str(), f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2), f2->GetParameter(3),
                f2->GetParameter(4));
    std::printf("[proton_tof_min] chi2<%s  pad1 flatgaus: flat=%.2f gaus_h=%.2f mean=%.2f sigma=%.3f\n", cutstr.c_str(),
                f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3));

    // Track-cut OOT suppression factors. Nothing on this canvas consumes them
    // (pads 1-4 subtract fitted backgrounds, not the OOT template), but these
    // are the numbers the parent macros fold into sflat_it / sflat_pk for their
    // other plot families, so print them for comparison.
    {
      const double fIT = oot_scale_f(h_proton_tof, hpt, 2);
      const double fPK = oot_scale_f(h_proton_tof, hpt, 3);
      auto regW = [](int r) { double w = 0.; for (const auto &iv : GREG_INT[r]) w += iv[1] - iv[0]; return w; };
      const double wOOT = regW(1), wIT = regW(2), wPK = regW(3);
      std::printf("[proton_tof_min] chi2<%s  OOT scale: f_IT=%.4f f_pk=%.4f  ->  sflat_it=%.4f sflat_pk=%.4f "
                  "(unscaled %.4f / %.4f)\n",
                  cutstr.c_str(), fIT, fPK, ((wOOT > 0.) ? wIT / wOOT : 0.) * fIT,
                  ((wOOT > 0.) ? wPK / wOOT : 0.) * fPK, (wOOT > 0.) ? wIT / wOOT : 0.,
                  (wOOT > 0.) ? wPK / wOOT : 0.);
    }

    // Write the canvas into the file, and drop a .png beside out_file.
    c->Write();
    {
      std::string png(out_file);
      const size_t dot = png.rfind(".root");
      if (dot != std::string::npos)
        png = png.substr(0, dot);
      png += "_chi2cut" + cutstr + ".png";
      c->SaveAs(png.c_str());
    }
    delete c;
    delete f1;
    delete f2;
  }

  // Keep the two input histograms in the file too -- handy for diffing runs.
  fout.cd();
  h_proton_tof->Write();
  for (int ic = ic_lo; ic <= ic_hi; ++ic)
    h_proton_track_tof[ic]->Write();
  fout.Close();

  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  std::cout << "[proton_tof_min] Done. Wrote " << out_file << "\n";
  std::printf("[proton_tof_min] Total time: %.1f s\n", elapsed);
}
