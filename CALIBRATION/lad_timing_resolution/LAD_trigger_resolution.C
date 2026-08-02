// LAD_trigger_resolution.C
//
// Trigger-timing plots for the LAD/SHMS/HMS trigger TDC signals.
//
// For every event with g.evtyp==1 it histograms the calibrated (_tdcTime) and
// raw (_tdcTimeRaw) TDC times of the trigger signals recorded in the SHMS and
// HMS trigger apparatuses (T.shms.* and T.hms.*).  All SHMS histograms are
// written into a "shms" directory and all HMS histograms into an "hms"
// directory of a single output ROOT file.
//
// Input can be either
//   * a text file (.dat/.txt) listing one ROOT-file path per line, or
//   * a path/glob passed straight to TChain::Add (wildcards allowed).
//
// Usage:
//   root -l -b -q 'LAD_trigger_resolution.C("files/my_runs.dat")'
//   root -l -b -q 'LAD_trigger_resolution.C("/path/to/LAD_*.root","trig_res.root")'
//

#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
// Default input ROOT file used when the macro is run with no arguments.
static const char *DEFAULT_INPUT =
    "/Users/lucasehinger/Library/CloudStorage/OneDrive-Personal/Documents/MIT/Research/LAD/LAD_coding/ROOTfiles/"
    "LAD_COIN/PRODUCTION/LAD_COIN_22614_1_1_-1.root";
static const char *DEFAULT_OUT_FILE = "LAD_trigger_resolution_output.root";
// SHMS triggers fire on evtyp==1, HMS triggers on evtyp==2, so each
// spectrometer's directory gets its own event-type cut.
static const char *SHMS_CUT = "g.evtyp==1";
static const char *HMS_CUT  = "g.evtyp==2";
static const int NBINS      = 500; // bins per histogram (range auto-set)

// Raw TDC channel -> ns conversion (t_{shms,hms}_tdcchanperns; CAEN 1190 LSB).
// Raw quantities are plotted as (raw * NS_PER_CH) so every axis is in ns.
static const double NS_PER_CH = 0.09766;

// Trigger signals to plot for each spectrometer.  All axes are in ns.
//   name       : signal name (branch is T.<prefix>.<name>_tdcTime[Raw])
//   raw        : also make the _tdcTimeRaw plot (converted to ns)
//   tmin,tmax  : fixed ns x-range for the _tdcTime plot (auto-ranged if tmin>=tmax)
//   rmin,rmax  : fixed x-range for the _tdcTimeRaw plot, given in RAW CHANNELS
//                (converted to ns internally; auto-ranged if rmin>=rmax)
struct TrigVar {
  std::string name;
  bool        raw;
  double      tmin, tmax;
  double      rmin, rmax;
};

// Difference plots (e.g. TrigLAD - pTRIG1) for each spectrometer.  Axes in ns.
//   name      : histogram name suffix
//   expr      : TTree::Draw expression
//   xmin,xmax : fixed x-range (auto-ranged if xmin>=xmax)
//   raw       : true  -> expr is in raw channels; scaled to ns and bounds are
//                        interpreted as raw channels
//               false -> expr is already in ns (_tdcTime); bounds are ns
struct DiffVar {
  std::string name;
  std::string expr;
  double      xmin, xmax;
  bool        raw;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Add every non-empty, non-comment line of a text file to the chain.
static int add_from_list(TChain *chain, const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Error: could not open file list '" << filename << "'\n";
    return 0;
  }
  int added = 0;
  std::string line;
  while (std::getline(infile, line)) {
    size_t a = line.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
      continue;
    std::string path = line.substr(a, line.find_last_not_of(" \t\r\n") - a + 1);
    if (path.empty() || path[0] == '#')
      continue;
    added += chain->Add(path.c_str());
  }
  return added;
}

// Book one histogram for T.<prefix>.<var> into the given directory, with the
// plotted quantity multiplied by `scale` (used to convert raw TDC channels to
// ns).  The bounds [xmin,xmax] are given in the variable's native units and are
// scaled the same way.  If xmax>xmin the range is fixed; otherwise TTree::Draw
// auto-determines it (only the bin count is fixed in the auto case).
static void book_trig_hist(TChain *chain, TDirectory *dir, const char *prefix, const char *cut,
                           const std::string &var, double xmin = 0., double xmax = 0., double scale = 1.0) {
  dir->cd();
  TString hname = Form("h_%s_%s", prefix, var.c_str());
  TString vexpr = (scale == 1.0) ? Form("T.%s.%s", prefix, var.c_str())
                                 : Form("T.%s.%s*%.10g", prefix, var.c_str(), scale);
  double lo = xmin * scale, hi = xmax * scale;
  TString expr;
  if (hi > lo)
    expr = Form("%s>>%s(%d,%g,%g)", vexpr.Data(), hname.Data(), NBINS, lo, hi);
  else
    expr = Form("%s>>%s(%d)", vexpr.Data(), hname.Data(), NBINS);

  chain->Draw(expr, cut, "goff");

  TH1 *h = dynamic_cast<TH1 *>(dir->Get(hname));
  if (!h) {
    std::cerr << "  Warning: no histogram produced for T." << prefix << "." << var
              << " (branch missing or no passing entries)\n";
    return;
  }

  h->SetTitle(Form("T.%s.%s;%s (ns);Counts", prefix, var.c_str(), var.c_str()));
  std::cout << "  T." << prefix << "." << var << ": " << (Long64_t)h->GetEntries() << " entries\n";
  h->Write();
}

// Fill and write all trigger histograms for one spectrometer.
static void fill_spectrometer(TChain *chain, TDirectory *dir, const char *prefix, const char *cut,
                              const std::vector<TrigVar> &vars) {
  std::cout << "Filling '" << prefix << "' trigger histograms (cut: " << cut << ")...\n";
  for (const TrigVar &v : vars) {
    book_trig_hist(chain, dir, prefix, cut, v.name + "_tdcTime", v.tmin, v.tmax); // already ns
    if (v.raw)
      book_trig_hist(chain, dir, prefix, cut, v.name + "_tdcTimeRaw", v.rmin, v.rmax, NS_PER_CH);
  }
}

// Book one raw-TDC difference histogram into the given directory. Fixed bounds
// send no-hit (zero) outliers to under/overflow, so no extra cut is needed
// beyond the directory's event-type cut.
static void book_diff_hist(TChain *chain, TDirectory *dir, const char *cut, const DiffVar &d) {
  dir->cd();
  double scale  = d.raw ? NS_PER_CH : 1.0;
  TString hname = Form("h_%s", d.name.c_str());
  TString vexpr = d.raw ? Form("(%s)*%.10g", d.expr.c_str(), NS_PER_CH) : TString(d.expr.c_str());
  double lo = d.xmin * scale, hi = d.xmax * scale;
  TString draw;
  if (hi > lo)
    draw = Form("%s>>%s(%d,%g,%g)", vexpr.Data(), hname.Data(), NBINS, lo, hi);
  else
    draw = Form("%s>>%s(%d)", vexpr.Data(), hname.Data(), NBINS);

  chain->Draw(draw, cut, "goff");

  TH1 *h = dynamic_cast<TH1 *>(dir->Get(hname));
  if (!h) {
    std::cerr << "  Warning: no histogram produced for " << d.expr << "\n";
    return;
  }
  h->SetTitle(Form("%s;%s (ns);Counts", d.name.c_str(), d.expr.c_str()));

  // Gaussian fit over the plotted range.  TH1::Fit only uses in-range bins
  // (1..NBINS), so the under/overflow entries are excluded automatically.
  double sigma = 0.;
  if (h->GetEntries() > 0) {
    h->Fit("gaus", "Q"); // quiet; full-axis range => no under/overflow bins
    if (TF1 *g = h->GetFunction("gaus")) {
      sigma            = g->GetParameter(2);
      double sigma_err = g->GetParError(2);
      double mu        = g->GetParameter(1);
      // Attach a label so the width is drawn (and saved) with the histogram.
      TPaveText *pt = new TPaveText(0.60, 0.70, 0.88, 0.88, "NDC");
      pt->SetBorderSize(1);
      pt->SetFillColor(0);
      pt->SetTextAlign(12);
      pt->AddText(Form("#mu = %.2f ns", mu));
      pt->AddText(Form("#sigma = %.3f #pm %.3f ns", sigma, sigma_err));
      pt->AddText(Form("FWHM = %.3f ns", 2.35482 * sigma));
      h->GetListOfFunctions()->Add(pt);
    }
  }

  std::cout << "  " << d.name << ": " << (Long64_t)h->GetEntries() << " entries, sigma = " << sigma << " ns\n";
  h->Write();
}

// Fill and write all difference histograms for one spectrometer.
static void fill_diffs(TChain *chain, TDirectory *dir, const char *cut, const std::vector<DiffVar> &diffs) {
  for (const DiffVar &d : diffs)
    book_diff_hist(chain, dir, cut, d);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void LAD_trigger_resolution(const char *input = DEFAULT_INPUT, const char *out_file = DEFAULT_OUT_FILE) {
  gROOT->SetBatch(kTRUE);

  // -------------------------------------------------------------------------
  // 1. Build the TChain from the requested ROOT files
  // -------------------------------------------------------------------------
  TChain *chain = new TChain("T");
  std::string in(input);
  bool is_list =
      (in.size() >= 4) && (in.compare(in.size() - 4, 4, ".dat") == 0 || in.compare(in.size() - 4, 4, ".txt") == 0);

  if (is_list)
    add_from_list(chain, in);
  else
    chain->Add(input); // single file or wildcard/glob

  Long64_t nentries = chain->GetEntries();
  std::cout << "[LAD_trigger_resolution] chained " << chain->GetNtrees() << " file(s), " << nentries << " entries\n";
  if (nentries == 0) {
    std::cerr << "Error: empty TChain -- check the input path/list.\n";
    return;
  }

  // -------------------------------------------------------------------------
  // 2. Variables to plot per spectrometer
  // -------------------------------------------------------------------------
  // SHMS trigger apparatus (T.shms.*): LAD triggers, SHMS pTRIG1/2, SHMS->LAD trigger.
  // Columns: {name, raw?, tmin, tmax, rmin, rmax}  (min>=max => auto-range)
  std::vector<TrigVar> shms_vars = {
      {"ltrig1LAD", true, 0, 0, 0, 0},
      {"ltrig2LAD", true, 0, 0, 0, 0},
      {"pTRIG1", true, 345, 355, 9100, 9700},
      {"pTRIG2", true, 333, 343, 9000, 9600},
      {"shmsTrigLAD", true, 1790, 1850, 15200, 16200}};

  // HMS trigger apparatus (T.hms.*): LAD triggers, HMS pTRIG3/4, HMS->LAD trigger.
  // Calibrated _tdcTime windows are similar to SHMS, but the raw TDC channels
  // are on different modules (~32500 for pTRIG3/4, ~15750 for hmsTrigLAD).
  std::vector<TrigVar> hms_vars = {
      {"ltrig1LAD", true, 0, 0, 0, 0},
      {"ltrig2LAD", true, 0, 0, 0, 0},
      {"pTRIG3", true, 345, 355, 31500, 33500},
      {"pTRIG4", true, 340, 350, 32000, 33000},
      {"hmsTrigLAD", true, 1790, 1900, 15200, 16200}};

  // Difference plots (TrigLAD - pTRIGn, pTRIGa - pTRIGb), all axes in ns.
  // Columns: {name, expr, xmin, xmax, raw}.  raw=true bounds are in channels.
  std::vector<DiffVar> shms_diffs = {
      // raw-TDC differences (channels -> ns)
      {"shms_TrigLAD_minus_pTRIG1_raw", "T.shms.shmsTrigLAD_tdcTimeRaw-T.shms.pTRIG1_tdcTimeRaw", 5000, 7000, true},
      {"shms_TrigLAD_minus_pTRIG2_raw", "T.shms.shmsTrigLAD_tdcTimeRaw-T.shms.pTRIG2_tdcTimeRaw", 5000, 7000, true},
      {"shms_pTRIG1_minus_pTRIG2_raw", "T.shms.pTRIG1_tdcTimeRaw-T.shms.pTRIG2_tdcTimeRaw", -200, 500, true},
      // calibrated _tdcTime differences (already ns).  These differ from the raw
      // differences above: the LAD-trigger and pTRIG channels have different TDC
      // reference times, which cancel only for pTRIG-pTRIG (not TrigLAD-pTRIG).
      {"shms_TrigLAD_minus_pTRIG1_time", "T.shms.shmsTrigLAD_tdcTime-T.shms.pTRIG1_tdcTime", 1420, 1540, false},
      {"shms_TrigLAD_minus_pTRIG2_time", "T.shms.shmsTrigLAD_tdcTime-T.shms.pTRIG2_tdcTime", 1420, 1540, false},
      {"shms_pTRIG1_minus_pTRIG2_time", "T.shms.pTRIG1_tdcTime-T.shms.pTRIG2_tdcTime", -20, 50, false}};

  std::vector<DiffVar> hms_diffs = {
      // raw-TDC differences (channels -> ns)
      {"hms_TrigLAD_minus_pTRIG3_raw", "T.hms.hmsTrigLAD_tdcTimeRaw-T.hms.pTRIG3_tdcTimeRaw", -18000, -15000, true},
      {"hms_TrigLAD_minus_pTRIG4_raw", "T.hms.hmsTrigLAD_tdcTimeRaw-T.hms.pTRIG4_tdcTimeRaw", -18000, -15000, true},
      {"hms_pTRIG3_minus_pTRIG4_raw", "T.hms.pTRIG3_tdcTimeRaw-T.hms.pTRIG4_tdcTimeRaw", -1000, 1000, true},
      // calibrated _tdcTime differences (already ns).  These differ from the raw
      // differences above: the LAD-trigger and pTRIG channels have different TDC
      // reference times, which cancel only for pTRIG-pTRIG (not TrigLAD-pTRIG).
      {"hms_TrigLAD_minus_pTRIG3_time", "T.hms.hmsTrigLAD_tdcTime-T.hms.pTRIG3_tdcTime", 1400, 1620, false},
      {"hms_TrigLAD_minus_pTRIG4_time", "T.hms.hmsTrigLAD_tdcTime-T.hms.pTRIG4_tdcTime", 1400, 1620, false},
      {"hms_pTRIG3_minus_pTRIG4_time", "T.hms.pTRIG3_tdcTime-T.hms.pTRIG4_tdcTime", -100, 100, false}};

  // -------------------------------------------------------------------------
  // 3. Output file with one directory per spectrometer
  // -------------------------------------------------------------------------
  TFile *fout = new TFile(out_file, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: could not create output file '" << out_file << "'\n";
    return;
  }

  TDirectory *d_shms = fout->mkdir("shms");
  TDirectory *d_hms  = fout->mkdir("hms");
  // Difference plots live in a "differences" sub-directory of each spectrometer.
  TDirectory *d_shms_diff = d_shms->mkdir("differences");
  TDirectory *d_hms_diff  = d_hms->mkdir("differences");

  fill_spectrometer(chain, d_shms, "shms", SHMS_CUT, shms_vars);
  fill_diffs(chain, d_shms_diff, SHMS_CUT, shms_diffs);
  fill_spectrometer(chain, d_hms, "hms", HMS_CUT, hms_vars);
  fill_diffs(chain, d_hms_diff, HMS_CUT, hms_diffs);

  // Histograms are Write()-ten into their directories inside book_trig_hist();
  // Close() persists the directory structure without duplicating them.
  fout->Close();
  std::cout << "[LAD_trigger_resolution] wrote " << out_file << "\n";
}
