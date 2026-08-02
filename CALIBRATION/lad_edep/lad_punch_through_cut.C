// lad_punch_through_cut.C
//
// Multi-threaded LAD punch-through analysis using RDataFrame.
//
// Reads per-plane HodoHit branches (one set per plane: 000, 001, 100, 101, 200).
// For front-back pairs (000&001, 100&101), finds hits on the same paddle where
// |YPos_front - YPos_back| < YPOS_DIFF_CUT, then plots:
//   - PT time (back_time - front_time)  [1D, per paddle + sums]
//   - PT time vs front Edep             [2D, per paddle]
//   - PT time vs back  Edep             [2D, per paddle]
//   - PT time vs front EdepAmp          [2D, per paddle]
//   - PT time vs back  EdepAmp          [2D, per paddle]
//
// Sum canvases exclude planes 100/101 paddles 1 and 9.
//
// Usage:
//   root -l -b -q 'lad_punch_through_cut.C("input.dat","out.root")'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
#if __has_include(<ROOT/RDFHelpers.hxx>)
#include <ROOT/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#elif __has_include(<ROOT/RDF/RDFHelpers.hxx>)
#include <ROOT/RDF/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#endif
#endif

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// =====================================================================
// Constants
// =====================================================================
const int N_PADDLES        = 11;
const double YPOS_DIFF_CUT = 20.0; // cm

const int NBINS_PT   = 150;
const double XMIN_PT = -5., XMAX_PT = 10.;
const int NBINS_EDEP   = 200;
const double XMIN_EDEP = 0., XMAX_EDEP = 100.;
const int NBINS_EDEP_AMP   = 400;
const double XMIN_EDEP_AMP = 0., XMAX_EDEP_AMP = 400.;

// const char *DEFAULT_DAT_FILE = "../files/run-lists/test.dat";
const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_22745-23590.dat";
const char *DEFAULT_OUT_FILE = "files/punch_through_cut/punch_through_cut.root";

// =====================================================================
// Per-pair tuple helpers
// =====================================================================
using RVd = ROOT::VecOps::RVec<double>;

struct PTuple {
  RVd pt, edep_f, edep_b, edepamp_f, edepamp_b;
};

static PTuple compute_pairs(double paddle, const RVd &pad_f, const RVd &t_f, const RVd &y_f, const RVd &e_f,
                            const RVd &ea_f, const RVd &pad_b, const RVd &t_b, const RVd &y_b, const RVd &e_b,
                            const RVd &ea_b) {
  PTuple r;
  for (size_t i = 0; i < pad_f.size(); ++i) {
    if (pad_f[i] != paddle)
      continue;
    for (size_t j = 0; j < pad_b.size(); ++j) {
      if (pad_b[j] != paddle)
        continue;
      if (std::abs(y_f[i] - y_b[j]) >= YPOS_DIFF_CUT)
        continue;
      r.pt.push_back(t_b[j] - t_f[i]);
      r.edep_f.push_back(e_f[i]);
      r.edep_b.push_back(e_b[j]);
      r.edepamp_f.push_back(ea_f[i]);
      r.edepamp_b.push_back(ea_b[j]);
    }
  }
  return r;
}

static PTuple compute_sum(bool excl, const RVd &pad_f, const RVd &t_f, const RVd &y_f, const RVd &e_f,
                          const RVd &ea_f, const RVd &pad_b, const RVd &t_b, const RVd &y_b, const RVd &e_b,
                          const RVd &ea_b) {
  PTuple r;
  for (size_t i = 0; i < pad_f.size(); ++i) {
    if (excl && (pad_f[i] == 1. || pad_f[i] == 9.))
      continue;
    for (size_t j = 0; j < pad_b.size(); ++j) {
      if (pad_b[j] != pad_f[i])
        continue;
      if (std::abs(y_f[i] - y_b[j]) >= YPOS_DIFF_CUT)
        continue;
      r.pt.push_back(t_b[j] - t_f[i]);
      r.edep_f.push_back(e_f[i]);
      r.edep_b.push_back(e_b[j]);
      r.edepamp_f.push_back(ea_f[i]);
      r.edepamp_b.push_back(ea_b[j]);
    }
  }
  return r;
}

// =====================================================================
void lad_punch_through_cut(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE) {

  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  ROOT::EnableImplicitMT();

  // ------------------------------------------------------------------
  // 1. TChain
  // ------------------------------------------------------------------
  TChain chain("T");
  int nfiles = 0;
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
      ++nfiles;
    }
  }
  if (nfiles == 0) {
    std::cerr << "empty file list\n";
    return;
  }
  const Long64_t nentries = chain.GetEntries(); // single scan; result cached for progress bar
  std::cout << "[lad_punch_through_cut] files: " << nfiles << "  entries: " << nentries << "\n";
  if (nentries == 0) {
    std::cerr << "empty chain\n";
    return;
  }

  // ------------------------------------------------------------------
  // 2. RDataFrame + aliases  (dots in branch names -> safe column names)
  // ------------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

#ifdef LAD_HAS_RDF_PROGRESSBAR
  ROOT::RDF::Experimental::AddProgressBar(df);
#else
  {
    const ULong64_t total = static_cast<ULong64_t>(nentries);
    const ULong64_t step  = std::max<ULong64_t>(1ULL, total / 200ULL);
    auto t0               = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());
    auto cnt              = df.Count();
    cnt.OnPartialResult(step, [total, t0](ULong64_t n) {
      double f = total ? double(n) / double(total) : 0.;
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - *t0).count();
      std::fprintf(stderr, "\r[lad_punch_through_cut] %6.2f%% (%llu/%llu) %.1fs eta %.1fs   ", 100. * f,
                   (unsigned long long)n, (unsigned long long)total, s, f > 0 ? s * (1. / f - 1.) : 0.);
      std::fflush(stderr);
    });
  }
#endif

  const char *const planes[5] = {"000", "001", "100", "101", "200"};
  for (int ip = 0; ip < 5; ++ip) {
    const std::string pl(planes[ip]);
    const std::string pfx = "P.ladhod." + pl + ".";
    df                    = df.Alias("p" + pl + "_time", pfx + "HodoHitTime");
    df                    = df.Alias("p" + pl + "_paddle_raw", pfx + "HodoHitPaddleNum");
    df                    = df.Define("p" + pl + "_paddle", "p" + pl + "_paddle_raw - 1.");
    df                    = df.Alias("p" + pl + "_edep", pfx + "HodoHitEdep_MeV");
    df                    = df.Alias("p" + pl + "_edepamp", pfx + "HodoHitEdepAmp_MeV");
    df                    = df.Alias("p" + pl + "_ypos", pfx + "HodoHitYPos");
  }

  // ------------------------------------------------------------------
  // 3. Column definitions
  //    pair 01: front=000, back=001
  //    pair 23: front=100, back=101  (sum excludes paddles 1 and 9)
  // ------------------------------------------------------------------
  for (int pair = 0; pair < 2; ++pair) {
    const std::string fp = (pair == 0) ? "000" : "100";
    const std::string bp = (pair == 0) ? "001" : "101";
    const std::string pp = (pair == 0) ? "01" : "23";
    const bool excl      = (pair == 1);

    // Per-paddle columns
    for (int pa = 0; pa < N_PADDLES; ++pa) {
      const double pv                     = pa;
      const std::string ps                = std::to_string(pa);
      const std::vector<std::string> srcs = {
          "p" + fp + "_paddle", "p" + fp + "_time", "p" + fp + "_ypos", "p" + fp + "_edep", "p" + fp + "_edepamp",
          "p" + bp + "_paddle", "p" + bp + "_time", "p" + bp + "_ypos", "p" + bp + "_edep", "p" + bp + "_edepamp"};

      df = df.Define(
          "pt_" + pp + "_b" + ps,
          [pv](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e, const RVd &f, const RVd &g,
               const RVd &h, const RVd &ii,
               const RVd &j) { return compute_pairs(pv, a, b, c, d, e, f, g, h, ii, j).pt; },
          srcs);
      df = df.Define(
          "edep_f_" + pp + "_b" + ps,
          [pv](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e, const RVd &f, const RVd &g,
               const RVd &h, const RVd &ii,
               const RVd &j) { return compute_pairs(pv, a, b, c, d, e, f, g, h, ii, j).edep_f; },
          srcs);
      df = df.Define(
          "edep_b_" + pp + "_b" + ps,
          [pv](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e, const RVd &f, const RVd &g,
               const RVd &h, const RVd &ii,
               const RVd &j) { return compute_pairs(pv, a, b, c, d, e, f, g, h, ii, j).edep_b; },
          srcs);
      df = df.Define(
          "edepamp_f_" + pp + "_b" + ps,
          [pv](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e, const RVd &f, const RVd &g,
               const RVd &h, const RVd &ii,
               const RVd &j) { return compute_pairs(pv, a, b, c, d, e, f, g, h, ii, j).edepamp_f; },
          srcs);
      df = df.Define(
          "edepamp_b_" + pp + "_b" + ps,
          [pv](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e, const RVd &f, const RVd &g,
               const RVd &h, const RVd &ii,
               const RVd &j) { return compute_pairs(pv, a, b, c, d, e, f, g, h, ii, j).edepamp_b; },
          srcs);
    }

    // Sum columns (all paddles; pair 23 excludes paddles 1 and 9)
    const std::vector<std::string> sum_srcs = {
        "p" + fp + "_paddle", "p" + fp + "_time", "p" + fp + "_ypos", "p" + fp + "_edep", "p" + fp + "_edepamp",
        "p" + bp + "_paddle", "p" + bp + "_time", "p" + bp + "_ypos", "p" + bp + "_edep", "p" + bp + "_edepamp"};

    df = df.Define("_ptup_" + pp + "_sum",
        [excl](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e,
               const RVd &f, const RVd &g, const RVd &h, const RVd &ii, const RVd &j) {
          return compute_sum(excl, a, b, c, d, e, f, g, h, ii, j);
        }, sum_srcs);
    df = df.Define("pt_"        + pp + "_sum", [](const PTuple &t){ return t.pt;        }, {"_ptup_"+pp+"_sum"});
    df = df.Define("edep_f_"    + pp + "_sum", [](const PTuple &t){ return t.edep_f;    }, {"_ptup_"+pp+"_sum"});
    df = df.Define("edep_b_"    + pp + "_sum", [](const PTuple &t){ return t.edep_b;    }, {"_ptup_"+pp+"_sum"});
    df = df.Define("edepamp_f_" + pp + "_sum", [](const PTuple &t){ return t.edepamp_f; }, {"_ptup_"+pp+"_sum"});
    df = df.Define("edepamp_b_" + pp + "_sum", [](const PTuple &t){ return t.edepamp_b; }, {"_ptup_"+pp+"_sum"});
  }

  // Totals (concatenate 01 + 23 sums)
  auto cat2 = [](const RVd &a, const RVd &b){ return ROOT::VecOps::Concatenate(a, b); };
  df = df.Define("pt_total",        cat2, {"pt_01_sum",        "pt_23_sum"});
  df = df.Define("edep_f_total",    cat2, {"edep_f_01_sum",    "edep_f_23_sum"});
  df = df.Define("edep_b_total",    cat2, {"edep_b_01_sum",    "edep_b_23_sum"});
  df = df.Define("edepamp_f_total", cat2, {"edepamp_f_01_sum", "edepamp_f_23_sum"});
  df = df.Define("edepamp_b_total", cat2, {"edepamp_b_01_sum", "edepamp_b_23_sum"});

  // ------------------------------------------------------------------
  // 4. Histogram booking
  // ------------------------------------------------------------------
  using RH1 = ROOT::RDF::RResultPtr<TH1D>;
  using RH2 = ROOT::RDF::RResultPtr<TH2D>;

  auto bk1 = [&](const std::string &col, const std::string &ttl, int nb, double lo, double hi) -> RH1 {
    return df.Histo1D({col.c_str(), ttl.c_str(), nb, lo, hi}, col);
  };
  auto bk2 = [&](const std::string &xc, const std::string &yc, const std::string &name, const std::string &ttl, int nx,
                 double lx, double hx, int ny, double ly,
                 double hy) -> RH2 { return df.Histo2D({name.c_str(), ttl.c_str(), nx, lx, hx, ny, ly, hy}, xc, yc); };

  std::vector<RH1> h_pt_01(N_PADDLES), h_pt_23(N_PADDLES);
  std::vector<RH2> h_pt_edep_f_01(N_PADDLES), h_pt_edep_b_01(N_PADDLES);
  std::vector<RH2> h_pt_edepamp_f_01(N_PADDLES), h_pt_edepamp_b_01(N_PADDLES);
  std::vector<RH2> h_pt_edep_f_23(N_PADDLES), h_pt_edep_b_23(N_PADDLES);
  std::vector<RH2> h_pt_edepamp_f_23(N_PADDLES), h_pt_edepamp_b_23(N_PADDLES);

  for (int pa = 0; pa < N_PADDLES; ++pa) {
    const std::string ps = std::to_string(pa);

    h_pt_01[pa] = bk1("pt_01_b" + ps, "PT 000-001 paddle " + ps + ";#Deltat (ns);Counts", NBINS_PT, XMIN_PT, XMAX_PT);
    h_pt_23[pa] = bk1("pt_23_b" + ps, "PT 100-101 paddle " + ps + ";#Deltat (ns);Counts", NBINS_PT, XMIN_PT, XMAX_PT);

    h_pt_edep_f_01[pa]    = bk2("pt_01_b" + ps, "edep_f_01_b" + ps, "h_pt_edep_f_01_b" + ps,
                                "PT vs Edep front 000-001 pad" + ps + ";#Deltat (ns);Edep front (MeV)", NBINS_PT, XMIN_PT,
                                XMAX_PT, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP);
    h_pt_edep_b_01[pa]    = bk2("pt_01_b" + ps, "edep_b_01_b" + ps, "h_pt_edep_b_01_b" + ps,
                                "PT vs Edep back 000-001 pad" + ps + ";#Deltat (ns);Edep back (MeV)", NBINS_PT, XMIN_PT,
                                XMAX_PT, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP);
    h_pt_edepamp_f_01[pa] = bk2("pt_01_b" + ps, "edepamp_f_01_b" + ps, "h_pt_edepamp_f_01_b" + ps,
                                "PT vs EdepAmp front 000-001 pad" + ps + ";#Deltat (ns);EdepAmp front (MeV)", NBINS_PT,
                                XMIN_PT, XMAX_PT, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP);
    h_pt_edepamp_b_01[pa] = bk2("pt_01_b" + ps, "edepamp_b_01_b" + ps, "h_pt_edepamp_b_01_b" + ps,
                                "PT vs EdepAmp back 000-001 pad" + ps + ";#Deltat (ns);EdepAmp back (MeV)", NBINS_PT,
                                XMIN_PT, XMAX_PT, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP);

    h_pt_edep_f_23[pa]    = bk2("pt_23_b" + ps, "edep_f_23_b" + ps, "h_pt_edep_f_23_b" + ps,
                                "PT vs Edep front 100-101 pad" + ps + ";#Deltat (ns);Edep front (MeV)", NBINS_PT, XMIN_PT,
                                XMAX_PT, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP);
    h_pt_edep_b_23[pa]    = bk2("pt_23_b" + ps, "edep_b_23_b" + ps, "h_pt_edep_b_23_b" + ps,
                                "PT vs Edep back 100-101 pad" + ps + ";#Deltat (ns);Edep back (MeV)", NBINS_PT, XMIN_PT,
                                XMAX_PT, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP);
    h_pt_edepamp_f_23[pa] = bk2("pt_23_b" + ps, "edepamp_f_23_b" + ps, "h_pt_edepamp_f_23_b" + ps,
                                "PT vs EdepAmp front 100-101 pad" + ps + ";#Deltat (ns);EdepAmp front (MeV)", NBINS_PT,
                                XMIN_PT, XMAX_PT, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP);
    h_pt_edepamp_b_23[pa] = bk2("pt_23_b" + ps, "edepamp_b_23_b" + ps, "h_pt_edepamp_b_23_b" + ps,
                                "PT vs EdepAmp back 100-101 pad" + ps + ";#Deltat (ns);EdepAmp back (MeV)", NBINS_PT,
                                XMIN_PT, XMAX_PT, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP);
  }

  RH1 h_pt_01_sum = bk1("pt_01_sum", "PT 000-001 sum;#Deltat (ns);Counts", NBINS_PT, XMIN_PT, XMAX_PT);
  RH1 h_pt_23_sum = bk1("pt_23_sum", "PT 100-101 sum;#Deltat (ns);Counts", NBINS_PT, XMIN_PT, XMAX_PT);
  RH1 h_pt_tot    = bk1("pt_total",  "PT total;#Deltat (ns);Counts",        NBINS_PT, XMIN_PT, XMAX_PT);

  // 2D sum and total histograms
  RH2 h_pt_edep_f_01_sum    = bk2("pt_01_sum","edep_f_01_sum",       "h_pt_edep_f_01_sum",
    "PT vs Edep front 000-001 sum;#Deltat (ns);Edep front (MeV)",       NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);
  RH2 h_pt_edep_f_23_sum    = bk2("pt_23_sum","edep_f_23_sum",       "h_pt_edep_f_23_sum",
    "PT vs Edep front 100-101 sum;#Deltat (ns);Edep front (MeV)",       NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);
  RH2 h_pt_edep_f_tot       = bk2("pt_total", "edep_f_total",        "h_pt_edep_f_total",
    "PT vs Edep front total;#Deltat (ns);Edep front (MeV)",             NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);

  RH2 h_pt_edep_b_01_sum    = bk2("pt_01_sum","edep_b_01_sum",       "h_pt_edep_b_01_sum",
    "PT vs Edep back 000-001 sum;#Deltat (ns);Edep back (MeV)",         NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);
  RH2 h_pt_edep_b_23_sum    = bk2("pt_23_sum","edep_b_23_sum",       "h_pt_edep_b_23_sum",
    "PT vs Edep back 100-101 sum;#Deltat (ns);Edep back (MeV)",         NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);
  RH2 h_pt_edep_b_tot       = bk2("pt_total", "edep_b_total",        "h_pt_edep_b_total",
    "PT vs Edep back total;#Deltat (ns);Edep back (MeV)",               NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,XMIN_EDEP,XMAX_EDEP);

  RH2 h_pt_edepamp_f_01_sum = bk2("pt_01_sum","edepamp_f_01_sum",    "h_pt_edepamp_f_01_sum",
    "PT vs EdepAmp front 000-001 sum;#Deltat (ns);EdepAmp front (MeV)", NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
  RH2 h_pt_edepamp_f_23_sum = bk2("pt_23_sum","edepamp_f_23_sum",    "h_pt_edepamp_f_23_sum",
    "PT vs EdepAmp front 100-101 sum;#Deltat (ns);EdepAmp front (MeV)", NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
  RH2 h_pt_edepamp_f_tot    = bk2("pt_total", "edepamp_f_total",     "h_pt_edepamp_f_total",
    "PT vs EdepAmp front total;#Deltat (ns);EdepAmp front (MeV)",       NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);

  RH2 h_pt_edepamp_b_01_sum = bk2("pt_01_sum","edepamp_b_01_sum",    "h_pt_edepamp_b_01_sum",
    "PT vs EdepAmp back 000-001 sum;#Deltat (ns);EdepAmp back (MeV)",   NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
  RH2 h_pt_edepamp_b_23_sum = bk2("pt_23_sum","edepamp_b_23_sum",    "h_pt_edepamp_b_23_sum",
    "PT vs EdepAmp back 100-101 sum;#Deltat (ns);EdepAmp back (MeV)",   NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
  RH2 h_pt_edepamp_b_tot    = bk2("pt_total", "edepamp_b_total",     "h_pt_edepamp_b_total",
    "PT vs EdepAmp back total;#Deltat (ns);EdepAmp back (MeV)",         NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);

  // ------------------------------------------------------------------
  // 5. Trigger event loop
  // ------------------------------------------------------------------
  std::cout << "[lad_punch_through_cut] Running event loop...\n";
  (void)h_pt_tot->GetEntries();
#ifndef LAD_HAS_RDF_PROGRESSBAR
  std::fprintf(stderr, "\n");
#endif

  // ------------------------------------------------------------------
  // 6. Write output
  // ------------------------------------------------------------------
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "cannot open output\n";
    return;
  }

  auto wc = [](TCanvas *c) {
    c->Write();
    delete c;
  };

  auto pad_canvas = [](const std::string &n, const std::string &t) {
    TCanvas *c = new TCanvas(n.c_str(), t.c_str(), 1600, 1000);
    c->Divide(4, 3);
    return c;
  };

  // Writes per-paddle canvases + 3-panel summary canvas for a pair of 2D histogram sets
  auto write2d = [&](TDirectory *parent, const std::string &dname, const std::string &label,
                     std::vector<RH2> &h01, std::vector<RH2> &h23,
                     RH2 &h_sum01, RH2 &h_sum23, RH2 &h_tot) {
    TDirectory *d = parent->mkdir(dname.c_str());
    d->cd();
    TCanvas *c01 = pad_canvas("c_" + dname + "_000_001", label + " 000-001");
    for (int p = 0; p < N_PADDLES; ++p) { c01->cd(p + 1); h01[p]->DrawCopy("COLZ"); }
    wc(c01);
    TCanvas *c23 = pad_canvas("c_" + dname + "_100_101", label + " 100-101");
    for (int p = 0; p < N_PADDLES; ++p) { c23->cd(p + 1); h23[p]->DrawCopy("COLZ"); }
    wc(c23);
    TCanvas *cs = new TCanvas(("c_" + dname + "_sums").c_str(), (label + " sums").c_str(), 2400, 800);
    cs->Divide(3, 1);
    cs->cd(1); h_sum01->DrawCopy("COLZ");
    cs->cd(2); h_sum23->DrawCopy("COLZ");
    cs->cd(3); h_tot->DrawCopy("COLZ");
    wc(cs);
  };

  // ---- 1D punchthrough time ----
  {
    TDirectory *d = fout.mkdir("punchthrough_time");
    d->cd();
    TCanvas *c01 = pad_canvas("c_pt_000_001", "PT time 000-001 per paddle");
    for (int p = 0; p < N_PADDLES; ++p) {
      c01->cd(p + 1);
      h_pt_01[p]->DrawCopy();
    }
    wc(c01);
    TCanvas *c23 = pad_canvas("c_pt_100_101", "PT time 100-101 per paddle");
    for (int p = 0; p < N_PADDLES; ++p) {
      c23->cd(p + 1);
      h_pt_23[p]->DrawCopy();
    }
    wc(c23);
    TCanvas *cs = new TCanvas("c_pt_sums", "PT time sums", 1600, 600);
    cs->Divide(3, 1);
    cs->cd(1);
    h_pt_01_sum->DrawCopy();
    cs->cd(2);
    h_pt_23_sum->DrawCopy();
    cs->cd(3);
    h_pt_tot->DrawCopy();
    wc(cs);
  }

  // ---- 2D PT vs Edep ----
  {
    TDirectory *d = fout.mkdir("pt_vs_edep");
    write2d(d, "front", "PT vs Edep front", h_pt_edep_f_01, h_pt_edep_f_23,
            h_pt_edep_f_01_sum, h_pt_edep_f_23_sum, h_pt_edep_f_tot);
    write2d(d, "back",  "PT vs Edep back",  h_pt_edep_b_01, h_pt_edep_b_23,
            h_pt_edep_b_01_sum, h_pt_edep_b_23_sum, h_pt_edep_b_tot);
  }

  // ---- 2D PT vs EdepAmp ----
  {
    TDirectory *d = fout.mkdir("pt_vs_edepamp");
    write2d(d, "front", "PT vs EdepAmp front", h_pt_edepamp_f_01, h_pt_edepamp_f_23,
            h_pt_edepamp_f_01_sum, h_pt_edepamp_f_23_sum, h_pt_edepamp_f_tot);
    write2d(d, "back",  "PT vs EdepAmp back",  h_pt_edepamp_b_01, h_pt_edepamp_b_23,
            h_pt_edepamp_b_01_sum, h_pt_edepamp_b_23_sum, h_pt_edepamp_b_tot);
  }

  fout.Close();
  std::cout << "[lad_punch_through_cut] Done. Wrote " << out_file << "\n";
}
