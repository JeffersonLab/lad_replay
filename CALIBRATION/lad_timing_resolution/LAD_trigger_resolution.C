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
// A "lad_triggers" sub-directory of each compares the LAD triggers recorded in
// the TS crates (pTRIG5/6) with those in the LAD crate (ltrig1/2LAD), using
// double differences that cancel the per-module 1190 window jitter (see
// fill_lad_triggers).  Raw cross-crate differences, including the TrigLAD-pTRIGn
// ones in "differences", are dominated by that ~25 ns window jitter.
//
// Input can be either
//   * a text file (.dat/.txt) listing one ROOT-file path per line, or
//   * a path/glob passed straight to TChain::Add (wildcards allowed).
//
// Usage:
//   root -l -b -q 'LAD_trigger_resolution.C("files/my_runs.dat")'
//   root -l -b -q 'LAD_trigger_resolution.C("/path/to/LAD_*.root","trig_res.root")'
//   root -l -b -q 'LAD_trigger_resolution.C("in.root","out.root",true)'   // LAD-trigger section only
//
// The TTree::Draw sections re-read the input once per histogram (~25 passes);
// lad_only=true skips them and runs only the single-pass RDataFrame section.
//

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
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
// LAD-trigger timing (multi-hit channels, RDataFrame)
// ---------------------------------------------------------------------------
// Every 1190 latches its trigger window on its own 25 ns clock, so a raw
// difference between channels in DIFFERENT modules (e.g. LAD crate ROC9 slot 31
// vs SHMS ROC2 slot 19) carries ~25 ns of window jitter and says nothing about
// the trigger logic.  The window cancels only within one module.  To compare
// crates, each side is referenced to a copy of the electron trigger in its own
// module, giving a double difference
//     DD = (A_LADcrate - refA_LADcrate) - (B_TScrate - refB_TScrate)
// where A and B are the same LAD trigger seen in the two crates.  DD is flat
// in the window jitter; any width left is real relative jitter between the two
// trigger paths (e.g. a 4 ns clock re-latch would show up as a 4 ns box).
//
// TS inputs (LAD wiki): TRIG1 SHMS 3/4, TRIG2 SHMS EL-REAL, TRIG3 HMS 3/4,
// TRIG4 HMS EL-REAL, TRIG5 LAD & (SHMS||HMS EL-REAL), TRIG6 LAD single.
// LAD crate (ROC9 slot 31): ltrig1LAD (LAD+REFBAR), ltrig2LAD (LAD
// n-multiplicity), shmsTrigLAD, hmsTrigLAD.  Same channels appear in T.shms
// and T.hms; the LAD-crate ones carry no reference subtraction.
//
// The LAD-crate 1190 records BOTH edges of each logic pulse (pairs a fixed
// width apart, stable over runs 22565-23750).  Only leading edges are kept.
using RVecD = ROOT::RVec<double>;

static double lad_edge_width(const std::string &ch) { // raw channels; 0 = single-edge channel
  if (ch == "shmsTrigLAD") return 793;
  if (ch == "hmsTrigLAD") return 773;
  if (ch == "ltrig1LAD") return 560;
  if (ch == "ltrig2LAD") return 557;
  return 0;
}
static const double EDGE_TOL = 30; // raw channels (~3 ns)
// Difference histograms: one bin per TDC channel, centred on whole channels, +-2.5 us.
static const int DD_NBINS  = 2 * 25600;
static const double DD_LO  = (-25600 - 0.5) * NS_PER_CH;
static const double DD_HI  = (25600 - 0.5) * NS_PER_CH;

// Leading edges of a hit list, in ns.  A hit is a trailing edge if another hit
// sits `width` (+-EDGE_TOL) channels before it.
static RVecD leading_edges_ns(const RVecD &hits, double width) {
  RVecD out;
  for (double x : hits) {
    if (x <= 0)
      continue;
    bool trailing = false;
    if (width > 0)
      for (double y : hits)
        if (y > 0 && std::abs(x - y - width) < EDGE_TOL) {
          trailing = true;
          break;
        }
    if (!trailing)
      out.push_back(x * NS_PER_CH);
  }
  return out;
}

// Signal of one channel: prefix ("shms"/"hms") + channel name.
struct LadSig {
  std::string prefix, name;
  std::string branch() const { return "T." + prefix + "." + name + "_vecTdcTimeRaw"; }
  std::string column() const { return "le_" + prefix + "_" + name; }
};

// (A_i - refA) - (B_j - refB) for all leading-edge pairs.  An empty ref name
// means no reference.  fit=true fits the core with a Gaussian.
struct LadDiff {
  std::string name, title;
  LadSig a, b;
  LadSig refa, refb;
  bool fit;
};

// Fit the core of a peaked histogram: Gaussian in +-fit_half ns around the
// maximum (the LAD-trigger peaks have a one-sided shoulder, so the window is
// kept tight), plus the FWHM read straight from the bins.  With a finite
// pk_hint the maximum is searched only within +-1 ns of it (for low-statistics
// copies of a histogram).  Returns the peak position.
static double annotate_peak(TH1 *h, double pk_hint = NAN, double fit_half = 0.6, double show_half = 15.) {
  if (h->GetEntries() <= 0)
    return NAN;
  int bmax = h->GetMaximumBin();
  if (std::isfinite(pk_hint)) {
    bmax = h->FindBin(pk_hint);
    for (int b = h->FindBin(pk_hint - 1.); b <= h->FindBin(pk_hint + 1.); ++b)
      if (h->GetBinContent(b) > h->GetBinContent(bmax))
        bmax = b;
  }
  double pk  = h->GetBinCenter(bmax);
  double max = h->GetBinContent(bmax);
  int lo = bmax, hi = bmax;
  while (lo > 1 && h->GetBinContent(lo - 1) > max / 2)
    --lo;
  while (hi < h->GetNbinsX() && h->GetBinContent(hi + 1) > max / 2)
    ++hi;
  double fwhm = h->GetBinLowEdge(hi + 1) - h->GetBinLowEdge(lo);
  TF1 g("g_core", "gaus", pk - fit_half, pk + fit_half);
  g.SetParameters(max, pk, fit_half / 2);
  h->Fit(&g, "QR0");
  if (TF1 *f = h->GetFunction("g_core"))
    f->ResetBit(TF1::kNotDraw);
  TPaveText *pt = new TPaveText(0.60, 0.70, 0.88, 0.88, "NDC");
  pt->SetBorderSize(1);
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(Form("peak = %.2f ns", pk));
  pt->AddText(Form("core #sigma = %.3f #pm %.3f ns", g.GetParameter(2), g.GetParError(2)));
  pt->AddText(Form("FWHM (bins) = %.2f ns", fwhm));
  h->GetListOfFunctions()->Add(pt);
  h->GetXaxis()->SetRangeUser(pk - show_half, pk + show_half);
  std::cout << "  " << h->GetName() << ": " << (Long64_t)h->GetEntries() << " entries, peak " << pk
            << " ns, core sigma " << g.GetParameter(2) << " ns, FWHM " << fwhm << " ns\n";
  return pk;
}

// Book and write the LAD-trigger histograms of one spectrometer.
//   evcut : event-type cut, e.g. "g.evtyp==1"
//   sigs  : single-channel leading-edge time plots
//   diffs : differences / double differences
static void fill_lad_triggers(TChain *chain, TDirectory *dir, const char *evcut, const std::vector<LadSig> &sigs,
                              const std::vector<LadDiff> &diffs) {
  std::cout << "Filling LAD-trigger histograms in " << dir->GetPath() << " (cut: " << evcut << ")...\n";
  auto has = [&](const LadSig &s) { return s.name.empty() || chain->GetBranch(s.branch().c_str()) != nullptr; };

  ROOT::RDataFrame rdf(*chain);
  ROOT::RDF::RNode df = rdf.Filter(evcut);
  // Leading-edge (ns) column for every channel that is used.
  std::set<std::string> defined;
  auto define = [&](const LadSig &s) {
    if (s.name.empty() || defined.count(s.column()))
      return;
    double w = lad_edge_width(s.name);
    df = df.Define(s.column(), [w](const RVecD &h) { return leading_edges_ns(h, w); }, {s.branch()});
    defined.insert(s.column());
  };
  df = df.Define("no_ref", [] { return RVecD{0.}; });
  auto col = [](const LadSig &s) { return s.name.empty() ? std::string("no_ref") : s.column(); };

  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
  enum { kNoFit, kFit, kFitAtPrevPeak };
  std::vector<int> fits;
  for (const LadSig &s : sigs) {
    if (!has(s)) {
      std::cerr << "  Warning: missing branch " << s.branch() << "\n";
      continue;
    }
    define(s);
    TString n = Form("h_%s_%s_leading", s.prefix.c_str(), s.name.c_str());
    hs.push_back(df.Histo1D({n.Data(), Form("T.%s.%s leading edges;raw time (ns);Hits", s.prefix.c_str(), s.name.c_str()),
                             4000, 0., 4000.},
                            s.column()));
    fits.push_back(kNoFit);
  }
  for (const LadDiff &d : diffs) {
    if (!has(d.a) || !has(d.b) || !has(d.refa) || !has(d.refb)) {
      std::cerr << "  Warning: missing branch for " << d.name << ", skipped\n";
      continue;
    }
    for (const LadSig *s : {&d.a, &d.b, &d.refa, &d.refb})
      define(*s);
    // Reference = earliest leading edge (one pulse per event for the e-trigger copies).
    auto dd = [](const RVecD &a, const RVecD &b, const RVecD &ra, const RVecD &rb) {
      RVecD out;
      if (ra.empty() || rb.empty())
        return out;
      double r = ROOT::VecOps::Min(ra) - ROOT::VecOps::Min(rb);
      for (double x : a)
        for (double y : b)
          out.push_back(x - y - r);
      return out;
    };
    std::vector<std::string> cols = {d.a.column(), d.b.column(), col(d.refa), col(d.refb)};
    std::string vn = "v_" + d.name;
    auto dfd = df.Define(vn, dd, cols);
    // One bin per TDC channel over the full 1190 window (all values are whole channels,
    // so wider-than-channel bins alias); the display range is set on the peak afterwards.
    hs.push_back(dfd.Histo1D({Form("h_%s", d.name.c_str()), (d.title + ";(ns);Pairs").c_str(), DD_NBINS, DD_LO, DD_HI}, vn));
    fits.push_back(d.fit ? kFit : kNoFit);
    if (d.fit) { // same, only events with exactly one leading edge in A and B (no pairing ambiguity)
      auto one = [](const RVecD &a, const RVecD &b) { return a.size() == 1 && b.size() == 1; };
      hs.push_back(dfd.Filter(one, {d.a.column(), d.b.column()})
                       .Histo1D({Form("h_%s_single", d.name.c_str()),
                                 (d.title + " [1 leading edge in each];(ns);Events").c_str(), DD_NBINS, DD_LO, DD_HI},
                                vn));
      fits.push_back(kFitAtPrevPeak);
    }
  }

  dir->cd();
  double prev_peak = NAN;
  for (size_t i = 0; i < hs.size(); ++i) {
    TH1D *h = static_cast<TH1D *>(hs[i]->Clone());
    if (fits[i] == kFit)
      prev_peak = annotate_peak(h);
    else if (fits[i] == kFitAtPrevPeak)
      annotate_peak(h, prev_peak);
    else
      std::cout << "  " << h->GetName() << ": " << (Long64_t)h->GetEntries() << " entries\n";
    h->Write();
  }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void LAD_trigger_resolution(const char *input = DEFAULT_INPUT, const char *out_file = DEFAULT_OUT_FILE,
                            bool lad_only = false) {
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

  if (!lad_only) {
    fill_spectrometer(chain, d_shms, "shms", SHMS_CUT, shms_vars);
    fill_diffs(chain, d_shms_diff, SHMS_CUT, shms_diffs);
    fill_spectrometer(chain, d_hms, "hms", HMS_CUT, hms_vars);
    fill_diffs(chain, d_hms_diff, HMS_CUT, hms_diffs);
  }

  // -------------------------------------------------------------------------
  // 4. LAD triggers: TS crate vs LAD crate (see fill_lad_triggers above)
  // -------------------------------------------------------------------------
  // LAD triggers reach the TS only in some run periods (pTRIG5 e.g. 23105-23108
  // and 23727-23750, pTRIG6 most runs >= 23105); ltrig1LAD barely fires.
  const LadSig none{"", ""};
  // shmsTrigLAD/hmsTrigLAD are the SHMS/HMS 3/4 triggers, so the matching TS-crate
  // refs are pTRIG1/pTRIG3 (EL-REAL refs pTRIG2/pTRIG4 kept for comparison).
  // SHMS: LAD crate ref = shmsTrigLAD, ROC2 ref = pTRIG1 (same 1190 as pTRIG5/6).
  const LadSig s_l1{"shms", "ltrig1LAD"}, s_l2{"shms", "ltrig2LAD"}, s_sT{"shms", "shmsTrigLAD"};
  const LadSig s_p1{"shms", "pTRIG1"}, s_p2{"shms", "pTRIG2"}, s_p5{"shms", "pTRIG5"}, s_p6{"shms", "pTRIG6"};
  const LadSig s_pT2{"shms", "pT2"};
  std::vector<LadSig> shms_lad_sigs = {s_l1, s_l2, s_sT, s_p5, s_p6};
  std::vector<LadDiff> shms_lad_diffs = {
      // within one module: window cancels
      {"shms_ltrig2LAD_minus_shmsTrigLAD", "LAD crate: ltrig2LAD - shmsTrigLAD", s_l2, s_sT, none, none, false},
      {"shms_pTRIG6_minus_pTRIG2", "ROC2: pTRIG6 - pTRIG2", s_p6, s_p2, none, none, false},
      {"shms_pTRIG5_minus_pTRIG2", "ROC2: pTRIG5 - pTRIG2", s_p5, s_p2, none, none, true},
      // across crates, raw: dominated by the two 1190 windows (compare TrigLAD - pTRIGn above)
      {"shms_ltrig2LAD_minus_pTRIG6_raw", "LAD crate - ROC2 (raw): ltrig2LAD - pTRIG6", s_l2, s_p6, none, none, false},
      {"shms_shmsTrigLAD_minus_pTRIG5_raw", "LAD crate - ROC2 (raw): shmsTrigLAD - pTRIG5", s_sT, s_p5, none, none, false},
      // across crates, double difference: windows cancel, LAD trigger seen in both crates
      {"shms_DD_ltrig2_pTRIG6_ref_pTRIG1", "(ltrig2LAD - shmsTrigLAD) - (pTRIG6 - pTRIG1)", s_l2, s_p6, s_sT, s_p1, true},
      {"shms_DD_ltrig2_pTRIG6_ref_pTRIG2", "(ltrig2LAD - shmsTrigLAD) - (pTRIG6 - pTRIG2)", s_l2, s_p6, s_sT, s_p2, true},
      {"shms_DD_ltrig2_pTRIG6_ref_pT2", "(ltrig2LAD - shmsTrigLAD) - (pTRIG6 - pT2)", s_l2, s_p6, s_sT, s_pT2, true},
      {"shms_DD_ltrig1_pTRIG6_ref_pTRIG1", "(ltrig1LAD - shmsTrigLAD) - (pTRIG6 - pTRIG1)", s_l1, s_p6, s_sT, s_p1, true}};

  // HMS: LAD crate ref = hmsTrigLAD; TS inputs are recorded in ROC1 (T.hms.pTRIGn,
  // ref pTRIG3 = HMS 3/4) and in ROC2 (T.shms.pTRIGn).  hT2 sits in another
  // 1190 than T.hms.pTRIG6, so it is NOT a valid ref for a double difference.
  const LadSig h_l2{"hms", "ltrig2LAD"}, h_hT{"hms", "hmsTrigLAD"};
  const LadSig h_p3{"hms", "pTRIG3"}, h_p4{"hms", "pTRIG4"}, h_p5{"hms", "pTRIG5"}, h_p6{"hms", "pTRIG6"};
  const LadSig s_p3{"shms", "pTRIG3"}, s_p4{"shms", "pTRIG4"};
  std::vector<LadSig> hms_lad_sigs = {h_l2, h_hT, h_p5, h_p6};
  std::vector<LadDiff> hms_lad_diffs = {
      {"hms_ltrig2LAD_minus_hmsTrigLAD", "LAD crate: ltrig2LAD - hmsTrigLAD", h_l2, h_hT, none, none, false},
      {"hms_pTRIG6_minus_pTRIG4", "ROC1: pTRIG6 - pTRIG4", h_p6, h_p4, none, none, false},
      {"hms_ROC2_pTRIG5_minus_pTRIG4", "ROC2: pTRIG5 - pTRIG4", s_p5, s_p4, none, none, true},
      {"hms_ltrig2LAD_minus_pTRIG6_raw", "LAD crate - ROC1 (raw): ltrig2LAD - pTRIG6", h_l2, h_p6, none, none, false},
      {"hms_DD_ltrig2_pTRIG6_ROC1_ref_pTRIG3", "(ltrig2LAD - hmsTrigLAD) - (pTRIG6 - pTRIG3) [ROC1]", h_l2, h_p6, h_hT, h_p3,
       true},
      {"hms_DD_ltrig2_pTRIG6_ROC2_ref_pTRIG3", "(ltrig2LAD - hmsTrigLAD) - (pTRIG6 - pTRIG3) [ROC2]", h_l2, s_p6, h_hT, s_p3,
       true},
      {"hms_DD_ltrig2_pTRIG6_ROC1_ref_pTRIG4", "(ltrig2LAD - hmsTrigLAD) - (pTRIG6 - pTRIG4) [ROC1]", h_l2, h_p6, h_hT, h_p4,
       true},
      {"hms_DD_ltrig2_pTRIG6_ROC2_ref_pTRIG4", "(ltrig2LAD - hmsTrigLAD) - (pTRIG6 - pTRIG4) [ROC2]", h_l2, s_p6, h_hT, s_p4,
       true}};

  // Electron triggers only: events where both the SHMS and HMS triggers were
  // recorded (coincidences, any event type).  The HMS-SHMS trigger separation
  // seen in the LAD crate is compared with the one seen in ROC2/ROC1, so the
  // LAD trigger logic drops out.  Works for every run period, but statistics
  // are low (~30 events per 4 GB segment).
  const LadSig h_sT{"hms", "shmsTrigLAD"}, h_p1{"hms", "pTRIG1"}, h_p2{"hms", "pTRIG2"};
  std::vector<LadDiff> ee_diffs = {
      {"ee_DD_ROC2", "(hmsTrigLAD - shmsTrigLAD) - (pTRIG3 - pTRIG1) [ROC2]", h_hT, s_p3, h_sT, s_p1, true},
      {"ee_DD_ROC1", "(hmsTrigLAD - shmsTrigLAD) - (pTRIG3 - pTRIG1) [ROC1]", h_hT, h_p3, h_sT, h_p1, true},
      {"ee_DD_ROC2_elreal", "(hmsTrigLAD - shmsTrigLAD) - (pTRIG4 - pTRIG2) [ROC2]", h_hT, s_p4, h_sT, s_p2, true},
      {"ee_DD_ROC1_elreal", "(hmsTrigLAD - shmsTrigLAD) - (pTRIG4 - pTRIG2) [ROC1]", h_hT, h_p4, h_sT, h_p2, true}};

  fill_lad_triggers(chain, d_shms->mkdir("lad_triggers"), SHMS_CUT, shms_lad_sigs, shms_lad_diffs);
  fill_lad_triggers(chain, d_hms->mkdir("lad_triggers"), HMS_CUT, hms_lad_sigs, hms_lad_diffs);
  fill_lad_triggers(chain, fout->mkdir("both_etrig"), "g.evtyp>0", {}, ee_diffs);

  // Histograms are Write()-ten into their directories inside book_trig_hist();
  // Close() persists the directory structure without duplicating them.
  fout->Close();
  std::cout << "[LAD_trigger_resolution] wrote " << out_file << "\n";
}
