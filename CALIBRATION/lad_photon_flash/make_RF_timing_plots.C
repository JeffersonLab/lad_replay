// make_RF_timing_plots.C
//
// RF-timing peak analysis, one output ROOT file for a whole list of runs.
//
// Files in the .dat list are grouped by run number, then ADJACENT runs are
// combined into groups of at least min_files files (default 15, adjustable via
// the 4th argument) so each group has enough statistics.  All files of a group
// are chained and analysed TOGETHER; there is one output point per group,
// plotted at the group's file-weighted mean run number.  For each group this
// macro builds the path-length-corrected, RF-corrected time-of-flight
// distribution (analogous to the "P_c_tof_corr_summary" canvas produced by
// lad_tof/lad_tof_fast.C, but using the goodhit_hit_tof_rfcorr_* branches
// instead of goodhit_hit_tof_*), then:
//   * PRIMARY: measures the phase of the 4 ns beam-bunch comb from the satellite
//     peaks in [-50, 10] ns via the first-harmonic (Rayleigh) estimator.  All
//     comb teeth share the main peak's phase mod 4 ns, so pooling them measures
//     the phase far more precisely (and more robustly at low stats) than fitting
//     the main peak.  -> comb phase phi vs run number.
//   * CROSS-CHECK: also fits a Gaussian to the main peak (tallest bin in
//     [10, 30] ns, fit over peak +/- 1.5 ns).  Its center, converted to a phase
//     (mod 4 ns), is overlaid on the SAME phase-vs-run canvas as a second set of
//     points, and its width sigma is plotted vs run on its own canvas.
//   * writes a per-group summary canvas (planes 000,001,100,101,200 + all-plane
//     total) showing each distribution with the main-peak fit and a text box
//     giving the harmonic phi, the fit phi, and the fit width.
//
// The run number is the first 5-digit token in the file name (e.g.
// LAD_COIN_22614_1_1_-1.root -> 22614).  Files with no 5-digit run number are
// skipped.  One point per group.
//
// Both spectrometers (P and H) are analysed.  Hits are taken with no tracking
// cut (the tracking-independent "all_hits" category of lad_tof_fast).
//
// Output ROOT file layout (<label> = run<NNNNN> or runs<first>-<last>):
//   per_group/<label>/{P,H}_c_tof_corr_summary  -- summary canvas per group
//   per_group/<label>/<spec>_tofcorr_<name>     -- the histograms
//   summary/{P,H}_c_peak_phase                  -- comb-phase-vs-run, 6 pads
//   summary/{P,H}_c_peak_width                  -- width-vs-run (sanity), 6 pads
//
// Usage:
//   root -l -b -q 'make_RF_timing_plots.C("files.dat","out.root")'
//   root -l -b -q 'make_RF_timing_plots.C("files.dat","out.root",8)'     // 8 threads
//   root -l -b -q 'make_RF_timing_plots.C("files.dat","out.root",8,20)'  // 20-file groups
// 3rd arg caps the RDataFrame worker threads (default 4; <=0 uses all cores).
// 4th arg is min_files per group (default 15; 1 disables combining).
// Each group (all its files chained) is processed with one RDataFrame under
// implicit MT, so the event loop is parallelised across cores while memory stays
// bounded to one group at a time.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cstdio>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

// =====================================================================
// Geometry / binning constants (kept identical to lad_tof_fast.C so the
// distribution is directly comparable to P_c_tof_corr_summary)
// =====================================================================
const int    N_PLANES = 5, N_PADDLES = 11, N_SPECS = 2;
const int    N_TARGETS = N_PLANES + 1;               // 5 planes + all-plane total
const double hodo_radii[N_PLANES] = {615., 655.6, 523., 563.6, 615.}; // cm
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const char specs[N_SPECS] = {'P', 'H'};

const int    NBINS_TCORR = 900;
const double XMIN_TCORR = -150., XMAX_TCORR = 300.;
const double FIT_HALF_WINDOW = 1.5;   // Gaussian fit range = peak +/- this (ns)
const double MIN_ENTRIES_FOR_FIT = 50.;
// The RF-corrected coincidence peak sits at ~20 ns in tof_rfcorr - L/c; search
// for the maximum only inside this band so the fit cannot latch onto a taller
// noise bin elsewhere in the (mostly flat) background.
const double PEAK_SEARCH_MIN = 10.;
const double PEAK_SEARCH_MAX = 30.;

// Primary observable: the common phase of the 4 ns beam-bunch comb, measured
// from the clean satellite (accidental) peaks that sit well away from the main
// bump.  All comb teeth share the main peak's phase mod RF_PERIOD, so pooling
// them gives the phase far more precisely than fitting the main peak alone.
// RF_PERIOD matches BG_PERIOD_NS used by lad_tof_fast's bgsub_tof.
const double RF_PERIOD    = 4.0;      // ns
const double PHASE_WIN_MIN =  -50.;   // satellite-peak window (exactly 15 periods,
const double PHASE_WIN_MAX =   10.;   // just below the main peak so the two techniques compare)

const char *DEFAULT_DAT_FILE =
  "../files/run-lists/all_C3_runlist_22745-23590.dat";
const char *DEFAULT_OUT_FILE = "files/root/RF_timing_plots_full.root";

// label for target index t: 0..4 -> plane name, 5 -> "total"
static std::string target_label(int t) {
  return (t < N_PLANES) ? std::string(plane_names[t]) : std::string("total");
}

// result of the per-distribution measurements
struct PeakFit {
  // First-harmonic phase of the 4 ns comb (the primary result).
  bool   phase_ok = false;
  double phase = 0., phase_err = 0.;     // comb phase mod RF_PERIOD, in (-T/2, T/2]
  // Main-peak Gaussian fit -- kept only as a width sanity check.
  bool   ok = false;
  double center = 0., center_err = 0.;   // Gaussian mean      (par 1)
  double height = 0., height_err = 0.;   // Gaussian amplitude (par 0)
  double width = 0., width_err = 0.;     // Gaussian sigma     (par 2)
};

// A batch of adjacent runs combined so it holds at least min_files files.
struct RunGroup {
  int       first = -1, last = -1;      // run-number span of the group
  long long run_sum = 0, fcount = 0;    // for the file-weighted mean run (x position)
  std::vector<std::string> files;       // all files in the group (chained together)
  double xrun() const { return fcount ? double(run_sum) / double(fcount) : double(first); }
};

// =====================================================================
// First-harmonic (Rayleigh) phase of a periodic comb, measured over [lo, hi].
//   C = sum c_i cos(2*pi*x_i/T),  S = sum c_i sin(2*pi*x_i/T)
//   phase   = (T/2pi) * atan2(S, C)                  in (-T/2, T/2]
//   phase_err = (T/2pi) * sqrt(N/2) / sqrt(C^2+S^2)  (N = total counts in window)
// This is the ML phase estimator for a sinusoidal modulation; it uses every
// count and degrades gracefully at low stats.  Precision is set by the comb's
// modulation depth (sqrt(C^2+S^2)), not the raw count.
// =====================================================================
static void harmonic_phase(const TH1D *h, double lo, double hi, double T,
                           double &phase, double &phase_err, bool &ok) {
  double C = 0., S = 0., N = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x < lo || x > hi) continue;
    const double c = h->GetBinContent(b);
    const double w = 2. * M_PI * x / T;
    C += c * std::cos(w);
    S += c * std::sin(w);
    N += c;
  }
  const double A = std::sqrt(C * C + S * S);
  if (N < MIN_ENTRIES_FOR_FIT || A <= 0.) { ok = false; phase = 0.; phase_err = 0.; return; }
  phase     = T / (2. * M_PI) * std::atan2(S, C);
  phase_err = T / (2. * M_PI) * std::sqrt(N / 2.) / A;
  ok = true;
}

// Wrap an absolute time into a comb phase in (-T/2, T/2] (e.g. the main-peak fit
// center -> its phase mod T, on the same footing as harmonic_phase).
static double wrap_to_period(double x, double T) {
  double p = std::fmod(x, T);
  if (p >   T / 2.) p -= T;
  if (p <= -T / 2.) p += T;
  return p;
}

// =====================================================================
// Extract the 5-digit run number from a file path (first such token).
// Returns -1 if none found.
// =====================================================================
static int parse_run_number(const std::string &path) {
  std::string base = gSystem->BaseName(path.c_str());
  std::smatch m;
  static const std::regex re("([0-9]{5})");
  if (std::regex_search(base, m, re)) return std::stoi(m[1].str());
  return -1;
}

// file-name stem (basename without .root), used only for log messages
static std::string file_stem(const std::string &path) {
  std::string base = gSystem->BaseName(path.c_str());
  size_t dot = base.rfind(".root");
  if (dot != std::string::npos) base = base.substr(0, dot);
  return base;
}

// =====================================================================
// Analyse one group: all files of the group (adjacent runs combined) are chained
// and analysed together (for statistics).  Builds the 12 corrected-tof
// histograms (6 per spec), measures the comb phase and the main-peak fit, writes
// the per-group summary canvases, returns the results.
// =====================================================================
static bool analyze_run(const std::vector<std::string> &paths, TDirectory *perrun,
                        const std::string &rlabel,
                        std::array<std::array<PeakFit, N_TARGETS>, N_SPECS> &fits) {
  using RVd = ROOT::VecOps::RVec<double>;

  // Keep only files that actually open and hold tree T.
  std::vector<std::string> good;
  for (const auto &p : paths) {
    TFile ftest(p.c_str(), "READ");
    if (ftest.IsZombie() || !ftest.Get("T"))
      std::cerr << "[make_RF_timing_plots] cannot open tree T in " << p << " -- skipping file\n";
    else
      good.push_back(p);
  }
  if (good.empty()) {
    std::cerr << "[make_RF_timing_plots] " << rlabel << ": no usable files -- skipping run\n";
    return false;
  }

  ROOT::RDataFrame rdf("T", good);   // all files for this run, chained
  ROOT::RDF::RNode df = rdf;

  // --- column definitions: per-plane corrected-tof sums + all-plane total ---
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    const std::string pfx = sp + ".ladhod.goodhit_";
    std::vector<std::string> plane_cols;
    for (int pl = 0; pl < N_PLANES; ++pl) {
      const std::string side = (pl % 2 == 0) ? "0" : "1";
      const double R = hodo_radii[pl];
      const double plv = static_cast<double>(pl);
      const bool excl = (pl == 2 || pl == 3);   // planes 100/101 drop paddles 1 and 9
      const std::string col = sp + "_tofcorr_p" + std::to_string(pl) + "_sum";
      df = df.Define(col,
        [R, plv, excl](const RVd &pln, const RVd &pdl, const RVd &yp, const RVd &tf) {
          RVd r;
          for (size_t i = 0; i < pln.size(); ++i) {
            if (pln[i] != plv) continue;
            if (excl && (pdl[i] == 1. || pdl[i] == 9.)) continue;
            const double dx  = 22. * (pdl[i] - 6.);
            const double p2d = std::sqrt(yp[i] * yp[i] + dx * dx);
            r.push_back(tf[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3);
          }
          return r;
        },
        {pfx + "plane_" + side, pfx + "paddle_" + side,
         pfx + "hit_ypos_" + side, pfx + "hit_tof_rfcorr_" + side});
      plane_cols.push_back(col);
    }
    df = df.Define(sp + "_tofcorr_total",
      [](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e) {
        auto r = ROOT::VecOps::Concatenate(a, b);
        r = ROOT::VecOps::Concatenate(r, c);
        r = ROOT::VecOps::Concatenate(r, d);
        return ROOT::VecOps::Concatenate(r, e);
      },
      plane_cols);
  }

  // --- book histograms (lazy) ---
  std::array<std::array<ROOT::RDF::RResultPtr<TH1D>, N_TARGETS>, N_SPECS> hbk;
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    for (int t = 0; t < N_TARGETS; ++t) {
      const std::string col = (t < N_PLANES)
        ? sp + "_tofcorr_p" + std::to_string(t) + "_sum"
        : sp + "_tofcorr_total";
      const std::string hn = sp + "_tofcorr_" + target_label(t) + "_" + rlabel;
      const std::string ht = sp + " tof corr " + target_label(t) +
                             " (" + rlabel + ");tof_{rfcorr}-L/c (ns);Counts";
      hbk[is][t] = df.Histo1D({hn.c_str(), ht.c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, col);
    }
  }

  // Trigger the (multi-threaded) event loop once for this run.
  (void)hbk[0][0]->GetEntries();

  // --- fit peaks + draw per-run summary canvases ---
  TDirectory *fdir = perrun->mkdir(rlabel.c_str());
  if (!fdir) fdir = perrun;   // fall back if the label already exists
  fdir->cd();

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    std::array<TH1D *, N_TARGETS> hc{};
    std::array<TF1 *, N_TARGETS> fc{};
    std::array<TPaveText *, N_TARGETS> pt{};

    for (int t = 0; t < N_TARGETS; ++t) {
      TH1D *h = static_cast<TH1D *>(hbk[is][t]->Clone());
      h->SetDirectory(nullptr);
      hc[t] = h;

      PeakFit pf;
      if (h->GetEntries() >= MIN_ENTRIES_FOR_FIT && h->Integral() > 0.) {
        // Find the tallest bin only within the search band, then restore the
        // full axis range for drawing.
        h->GetXaxis()->SetRangeUser(PEAK_SEARCH_MIN, PEAK_SEARCH_MAX);
        const int    mbin  = h->GetMaximumBin();
        const double xpeak = h->GetBinCenter(mbin);
        const double ypeak = h->GetBinContent(mbin);
        h->GetXaxis()->SetRange(0, 0);       // reset zoom -> full range
        TF1 *g = new TF1((sp + "_gaus_" + target_label(t) + "_" + rlabel).c_str(), "gaus",
                         xpeak - FIT_HALF_WINDOW, xpeak + FIT_HALF_WINDOW);
        g->SetParameters(ypeak, xpeak, 0.5);
        g->SetParLimits(2, 1e-3, 20.);       // keep sigma positive/bounded
        g->SetLineColor(kRed);
        const int status = h->Fit(g, "RQ");  // R: use g's range; Q: quiet
        fc[t] = g;
        if (status == 0) {
          pf.ok         = true;
          pf.center     = g->GetParameter(1);
          pf.center_err = g->GetParError(1);
          pf.height     = g->GetParameter(0);
          pf.height_err = g->GetParError(0);
          pf.width      = g->GetParameter(2);
          pf.width_err  = g->GetParError(2);
        }
      }

      // Primary result: comb phase from the satellite peaks (pooled, robust).
      harmonic_phase(h, PHASE_WIN_MIN, PHASE_WIN_MAX, RF_PERIOD,
                     pf.phase, pf.phase_err, pf.phase_ok);
      fits[is][t] = pf;

      // Stamp the results onto the pad so they travel with the file
      // (independent of the viewer's gStyle OptFit setting).
      TPaveText *box = new TPaveText(0.13, 0.70, 0.45, 0.89, "NDC");
      box->SetFillColor(0);
      box->SetBorderSize(1);
      box->SetTextAlign(12);
      box->SetTextSize(0.033);
      char line[128];
      if (pf.phase_ok) {
        std::snprintf(line, sizeof line, "comb #phi = %.3f #pm %.3f ns", pf.phase, pf.phase_err); box->AddText(line);
      } else {
        box->AddText("comb #phi: too few satellites");
      }
      if (pf.ok) {
        std::snprintf(line, sizeof line, "fit #phi = %.3f #pm %.3f ns",
                      wrap_to_period(pf.center, RF_PERIOD), pf.center_err); box->AddText(line);
        std::snprintf(line, sizeof line, "fit #sigma = %.3f #pm %.3f ns", pf.width, pf.width_err); box->AddText(line);
      } else {
        box->AddText("fit failed / too few entries");
      }
      pt[t] = box;
    }

    // Summary canvas: 3x2 -> planes 000,001,100,101,200 then all-plane total.
    TCanvas *cs = new TCanvas((sp + "_c_tof_corr_summary").c_str(),
                              (sp + " tof corr summary (" + rlabel + ")").c_str(), 1600, 1000);
    cs->Divide(3, 2);
    for (int t = 0; t < N_TARGETS; ++t) {
      cs->cd(t + 1);
      hc[t]->Draw();
      if (fc[t]) fc[t]->Draw("SAME");
      pt[t]->Draw();
    }
    fdir->cd();
    cs->Write();
    for (int t = 0; t < N_TARGETS; ++t) hc[t]->Write();  // keep the fitted histos too
    delete cs;
    for (int t = 0; t < N_TARGETS; ++t) { delete pt[t]; delete fc[t]; delete hc[t]; }
  }
  return true;
}

// =====================================================================
// Main
// =====================================================================
void make_RF_timing_plots(const char *dat_file = DEFAULT_DAT_FILE,
                          const char *out_file = DEFAULT_OUT_FILE,
                          int nthreads = 4,
                          int min_files = 15) {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1110);
  if (min_files < 1) min_files = 1;   // 1 -> one group per run (no combining)

  if (nthreads > 0) { ROOT::EnableImplicitMT(nthreads);
    std::cout << "[make_RF_timing_plots] implicit MT: " << nthreads << " threads\n"; }
  else              { ROOT::EnableImplicitMT();
    std::cout << "[make_RF_timing_plots] implicit MT: all cores\n"; }

  // ---- read the file list, grouping files by run number ----
  // std::map keeps runs sorted ascending (nice x-order for the graphs).  Files
  // whose name has no 5-digit run number cannot be combined and are skipped.
  std::map<int, std::vector<std::string>> runs;
  size_t nfiles = 0;
  {
    std::ifstream fin(dat_file);
    if (!fin.is_open()) { std::cerr << "cannot open list " << dat_file << "\n"; return; }
    std::string ln;
    while (std::getline(fin, ln)) {
      size_t a = ln.find_first_not_of(" \t\r\n");
      if (a == std::string::npos) continue;
      std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
      if (p.empty() || p[0] == '#') continue;
      const int run = parse_run_number(p);
      if (run < 0) {
        std::cerr << "[make_RF_timing_plots] WARNING: no 5-digit run number in "
                  << file_stem(p) << " -- skipping\n";
        continue;
      }
      runs[run].push_back(p);
      ++nfiles;
    }
  }
  if (runs.empty()) { std::cerr << "no usable files in " << dat_file << "\n"; return; }
  std::cout << "[make_RF_timing_plots] " << nfiles << " file(s) in "
            << runs.size() << " run(s) to process\n";

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) { std::cerr << "cannot open output " << out_file << "\n"; return; }
  TDirectory *pergroup = fout.mkdir("per_group");

  // ---- combine adjacent runs into groups of >= min_files files ----
  // Runs are already sorted ascending; walk them, accumulating files until the
  // group reaches min_files, then start a new group.  A trailing under-filled
  // group is merged back into the previous one so every group has >= min_files
  // (unless the whole list has fewer files than that).
  std::vector<RunGroup> groups;
  {
    RunGroup cur;
    for (const auto &kv : runs) {
      if (cur.first < 0) cur.first = kv.first;
      cur.last = kv.first;
      const long long n = (long long)kv.second.size();
      for (const auto &f : kv.second) cur.files.push_back(f);
      cur.run_sum += (long long)kv.first * n;
      cur.fcount  += n;
      if ((int)cur.files.size() >= min_files) { groups.push_back(std::move(cur)); cur = RunGroup(); }
    }
    if (!cur.files.empty()) {                       // trailing partial group
      if (!groups.empty() && (int)cur.files.size() < min_files) {
        RunGroup &prev = groups.back();             // merge it back into the previous group
        prev.last     = cur.last;
        prev.run_sum += cur.run_sum;
        prev.fcount  += cur.fcount;
        for (auto &f : cur.files) prev.files.push_back(std::move(f));
      } else {
        groups.push_back(std::move(cur));
      }
    }
  }
  std::cout << "[make_RF_timing_plots] combined into " << groups.size()
            << " group(s) of >= " << min_files << " file(s)\n";

  // per [spec][target] accumulators of (x=mean run, value, error).
  //   harmonic comb phase (primary technique): gxp / gyp / gep
  //   fit-derived phase (main-peak center mod T): gxf / gyf / gef
  //   fit width (sanity check):                  gxw / gyw / gew
  // Kept separate because the harmonic phase and the main-peak fit can
  // succeed/fail independently, so their point lists may differ.
  std::array<std::array<std::vector<double>, N_TARGETS>, N_SPECS> gxp, gyp, gep;
  std::array<std::array<std::vector<double>, N_TARGETS>, N_SPECS> gxf, gyf, gef;
  std::array<std::array<std::vector<double>, N_TARGETS>, N_SPECS> gxw, gyw, gew;

  // ---- process each group (all its files combined) ----
  for (size_t ig = 0; ig < groups.size(); ++ig) {
    const RunGroup &g = groups[ig];
    const std::string rlabel = (g.first == g.last)
      ? "run" + std::to_string(g.first)
      : "runs" + std::to_string(g.first) + "-" + std::to_string(g.last);
    const double xrun = g.xrun();
    std::cout << "[make_RF_timing_plots] (" << (ig + 1) << "/" << groups.size()
              << ") " << rlabel << "  (" << g.files.size() << " file(s), x=" << xrun << ")"
              << std::endl;

    std::array<std::array<PeakFit, N_TARGETS>, N_SPECS> fits;
    if (!analyze_run(g.files, pergroup, rlabel, fits)) continue;

    for (int is = 0; is < N_SPECS; ++is)
      for (int t = 0; t < N_TARGETS; ++t) {
        const PeakFit &pf = fits[is][t];
        if (pf.phase_ok) {
          gxp[is][t].push_back(xrun);
          gyp[is][t].push_back(pf.phase); gep[is][t].push_back(pf.phase_err);
        }
        if (pf.ok) {
          // fit center -> phase mod T; the constant shift leaves the error alone
          gxf[is][t].push_back(xrun);
          gyf[is][t].push_back(wrap_to_period(pf.center, RF_PERIOD));
          gef[is][t].push_back(pf.center_err);
          gxw[is][t].push_back(xrun);
          gyw[is][t].push_back(pf.width); gew[is][t].push_back(pf.width_err);
        }
      }
  }

  // ---- build peak-phase / peak-width vs run-number graphs ----
  TDirectory *sdir = fout.mkdir("summary");

  auto make_graph = [](const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<double> &ey, const std::string &name,
                       const std::string &title, int mstyle, int color) -> TGraphErrors * {
    const int n = static_cast<int>(x.size());
    std::vector<double> ex(n, 0.);
    TGraphErrors *g = (n > 0)
      ? new TGraphErrors(n, const_cast<double *>(x.data()), const_cast<double *>(y.data()),
                         ex.data(), const_cast<double *>(ey.data()))
      : new TGraphErrors();
    g->SetName(name.c_str());
    g->SetTitle(title.c_str());
    g->SetMarkerStyle(mstyle);
    g->SetMarkerSize(0.9);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    return g;
  };

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);

    TCanvas *cc = new TCanvas((sp + "_c_peak_phase").c_str(),
                              (sp + " comb phase vs run").c_str(), 1600, 1000);
    TCanvas *cw = new TCanvas((sp + "_c_peak_width").c_str(),
                              (sp + " peak width vs run").c_str(), 1600, 1000);
    cc->Divide(3, 2);
    cw->Divide(3, 2);

    // TMultiGraphs own their added graphs; delete the mg (not the graphs).
    std::vector<TMultiGraph *> mgs;
    std::vector<TLegend *>     legs;
    std::vector<TGraphErrors *> wgraphs;
    for (int t = 0; t < N_TARGETS; ++t) {
      const std::string lab = target_label(t);

      // Phase pad: two techniques overlaid.
      TGraphErrors *gh = make_graph(gxp[is][t], gyp[is][t], gep[is][t],   // harmonic (satellites)
        sp + "_peakphase_" + lab, "", 20, kBlue + 1);
      TGraphErrors *gf = make_graph(gxf[is][t], gyf[is][t], gef[is][t],   // main-peak fit
        sp + "_fitphase_" + lab, "", 24, kRed + 1);

      TMultiGraph *mg = new TMultiGraph();
      mg->SetName((sp + "_mgphase_" + lab).c_str());
      mg->SetTitle((sp + " comb phase " + lab + ";run number;comb phase #phi (ns, mod 4)").c_str());
      if (gh->GetN() > 0) mg->Add(gh, "P"); else { delete gh; gh = nullptr; }
      if (gf->GetN() > 0) mg->Add(gf, "P"); else { delete gf; gf = nullptr; }

      cc->cd(t + 1);
      if (mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetEntries() > 0) {
        mg->Draw("AP");
        TLegend *lg = new TLegend(0.50, 0.78, 0.88, 0.89);
        lg->SetTextSize(0.03);
        if (gh) lg->AddEntry(gh, "satellite comb (Rayleigh)", "P");
        if (gf) lg->AddEntry(gf, "main-peak fit", "P");
        lg->Draw();
        legs.push_back(lg);
      }
      mgs.push_back(mg);

      // Width pad: single technique (fit sigma).
      TGraphErrors *gw = make_graph(gxw[is][t], gyw[is][t], gew[is][t],
        sp + "_peakwidth_" + lab, sp + " peak width " + lab + ";run number;peak width #sigma (ns)",
        20, kBlue + 1);
      cw->cd(t + 1); if (gw->GetN() > 0) gw->Draw("AP");
      wgraphs.push_back(gw);
    }
    // Only the combined canvases are stored (individual graphs are not written).
    sdir->cd();
    cc->Write();
    cw->Write();
    delete cc;
    delete cw;
    for (auto *mg : mgs) delete mg;      // also deletes the phase graphs it owns
    for (auto *lg : legs) delete lg;
    for (auto *g : wgraphs) delete g;
  }

  fout.Close();
  std::cout << "[make_RF_timing_plots] Done. Wrote " << out_file << "\n";
}
