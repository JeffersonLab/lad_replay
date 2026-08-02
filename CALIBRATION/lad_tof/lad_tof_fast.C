// lad_tof_fast.C
//
// Multi-threaded LAD ToF analysis using RDataFrame.
//
// Branch convention (hcana):
//   X.ladhod.goodhit_<var>_0      -> arrays for planes 0,2,4  (X = P or H)
//   X.ladhod.goodhit_<var>_1      -> arrays for planes 1,3
//   X.ladhod.goodhit_chiSquare    -> per-hit chi-squared, same length as _0 and _1 arrays
//
// Plane numbering vs display name:
//   index 0 -> 000,  1 -> 001,  2 -> 100,  3 -> 101,  4 -> 200
//
// Sum plots exclude plane 100 (idx 2) and plane 101 (idx 3) paddles 1 and 9.
//
// Tracking variants (per-hit chiSquare branches, X = P or H):
//   standard         -> X.ladhod.goodhit_chiSquare
//   xz               -> X.ladhod.goodhit_chiSquare_xz
//   noTrackVertex    -> X.ladhod.goodhit_chiSquare_noTrackVertex
//   noTrackVertex_xz -> X.ladhod.goodhit_chiSquare_noTrackVertex_xz
// Only variants present in the data (for both spectrometers) are produced; the
// standard tracking is always assumed present.
//
// Each output directory contains:
//   all_hits              -- no cut (tracking-independent)
//   <variant>/has_track   -- that variant's chiSquare < 100 per hit
//   <variant>/no_track    -- that variant's chiSquare >= 100 per hit
//
// Usage:
//   root -l -b -q 'lad_tof_fast.C("input.dat","out.root")'
//   root -l -b -q 'lad_tof_fast.C("input.dat","out.root",8)'  // 8 MT threads
// The 3rd arg caps RDataFrame worker threads (default 4); each thread holds its
// own copy of every histogram, so lower it if the job is OOM-killed, raise it if
// there is spare RAM.  <=0 uses all cores.
//

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TROOT.h>
#include <cmath>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
#if __has_include(<ROOT/RDFHelpers.hxx>)
#include <ROOT/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#elif __has_include(<ROOT/RDF/RDFHelpers.hxx>)
#include <ROOT/RDF/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#endif
#endif

#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// =====================================================================
// Histogram binning
// =====================================================================
const int NBINS_TOF=200;      const double XMIN_TOF=0.,      XMAX_TOF=100.;
const int NBINS_YPOS=200;     const double XMIN_YPOS=-100.,  XMAX_YPOS=100.;
const int NBINS_HITTIME=200;  const double XMIN_HITTIME=1550.,XMAX_HITTIME=2050.;
const int NBINS_EDEP=200;     const double XMIN_EDEP=0.,     XMAX_EDEP=100.;
const int NBINS_EDEP_AMP=300; const double XMIN_EDEP_AMP=0., XMAX_EDEP_AMP=300;
const int NBINS_PT=150;       const double XMIN_PT=-5.,      XMAX_PT=10.; // punch throughs
const int NBINS_TCORR=900;    const double XMIN_TCORR=-150.,  XMAX_TCORR=300.;

// N_TRACKS is the max # of tracking variants; N_CATS = 1 (all_hits) + 2*N_TRACKS
// (has_track/no_track per variant).  Both are compile-time capacities for the
// histogram-holder arrays; the actually-used counts (ntracks, ncats) are runtime.
const int N_PLANES=5, N_PADDLES=11, N_SPECS=2;
const int N_TRACKS=4, N_CATS=1+2*N_TRACKS;
const double BG_PERIOD_NS  = 4.0; // N: period of the repeating background in ns
const int    PROTON_REBIN  = 15;  // rebin factor applied before sideband-subtracted efficiency
const double hodo_radii[N_PLANES]={615.,655.6,523.,563.6,615.}; // cm
const char* const plane_names[N_PLANES]={"000","001","100","101","200"};
const std::array<char,N_SPECS> specs={'P','H'};

const char *DEFAULT_DAT_FILE="../files/run-lists/all_C3_runlist_22745-23590.dat";
// const char *DEFAULT_DAT_FILE="../files/run-lists/all_C3_runlist_good_only.dat";
const char *DEFAULT_OUT_FILE="files/root_fast/timing_C3_22745-23590_PH_all_5k.root";

// =====================================================================

struct VarConfig { std::string name; int nbins; double xmin,xmax; };

// Background-subtract a corrected-tof histogram.
// The first 10*period_ns of the x-range is used as background.
// That window is folded modulo period_ns (10 full periods), averaged,
// then the repeating pattern is subtracted from the entire x-range.
TH1D* bgsub_tof(const TH1D* h, double period_ns = BG_PERIOD_NS) {
  TH1D* out = (TH1D*)h->Clone((std::string(h->GetName())+"_bgsub").c_str());
  out->SetTitle((std::string(h->GetTitle())+" (bgsub)").c_str());

  int    nbins   = h->GetNbinsX();
  double xmin    = h->GetXaxis()->GetXmin();
  double bw      = h->GetBinWidth(1);
  double bg_end  = xmin + 10. * period_ns;
  double n_per   = 10.;

  int N_tmpl = std::max(1, (int)std::round(period_ns / bw));
  std::vector<double> tmpl(N_tmpl, 0.);

  for (int b = 1; b <= nbins; ++b) {
    double x = h->GetBinCenter(b);
    if (x >= bg_end) continue;
    double rel = x - xmin;
    int ti = (int)(std::fmod(rel, period_ns) / bw);
    if (ti < 0 || ti >= N_tmpl) continue;
    tmpl[ti] += h->GetBinContent(b);
  }
  for (int i = 0; i < N_tmpl; ++i) tmpl[i] /= n_per;

  for (int b = 1; b <= nbins; ++b) {
    double x   = h->GetBinCenter(b);
    double rel = x - xmin;
    int ti = ((int)(std::fmod(rel, period_ns) / bw) % N_tmpl + N_tmpl) % N_tmpl;
    out->SetBinContent(b, h->GetBinContent(b) - tmpl[ti]);
  }
  return out;
}

// Sideband-subtract a histogram by computing the mean bin content
// in [xmin, xmin+sideband_ns] and subtracting that constant from all bins.
TH1D* flat_bgsub(const TH1D* h, double sideband_ns = 50.) {
  TH1D* out = (TH1D*)h->Clone((std::string(h->GetName())+"_sb").c_str());
  out->SetTitle((std::string(h->GetTitle())+" (sb-sub)").c_str());
  double bg_end = h->GetXaxis()->GetXmin() + sideband_ns;
  int n_sb = 0; double sum_sb = 0.;
  for(int b = 1; b <= h->GetNbinsX(); ++b){
    if(h->GetBinCenter(b) >= bg_end) break;
    sum_sb += h->GetBinContent(b); ++n_sb;
  }
  double bg = (n_sb > 0) ? sum_sb / n_sb : 0.;
  for(int b = 1; b <= h->GetNbinsX(); ++b)
    out->SetBinContent(b, h->GetBinContent(b) - bg);
  return out;
}

// A tracking variant: 'dir' is its output sub-directory name; 'tsuf' is the
// suffix appended to "goodhit_chiSquare" (empty for the standard/full tracking).
struct Track { std::string dir, tsuf; };

// chi_mode: 0=all hits, 1=chiSquare<100, 2=chiSquare>=100.
// 'tsuf' selects which tracking variant's chiSquare column drives the cut.
struct Cat { std::string suf, dir; int chi_mode; std::string tsuf; };

void lad_tof_fast(const char *dat_file=DEFAULT_DAT_FILE, const char *out_file=DEFAULT_OUT_FILE,
                  int nthreads=100) {

  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  // RDataFrame clones every booked histogram per worker thread, so peak memory
  // scales with the thread count.  With 4 tracking variants (9 categories) the
  // full result set is large; cap the threads to avoid the OOM killer (which
  // shows up as a silent SIGKILL / exit 9) on memory-limited farm nodes.
  // Pass nthreads<=0 to use all available cores.
  if (nthreads>0) { ROOT::EnableImplicitMT(nthreads);
    std::cout<<"[lad_tof_fast] implicit MT: "<<nthreads<<" threads\n"; }
  else            { ROOT::EnableImplicitMT();
    std::cout<<"[lad_tof_fast] implicit MT: all cores\n"; }

  // ---------------------------------------------------------------
  // 1. TChain
  // ---------------------------------------------------------------
  TChain chain("T");
  {
    std::ifstream fin(dat_file);
    if (!fin.is_open()) { std::cerr<<"cannot open "<<dat_file<<"\n"; return; }
    std::string ln;
    while (std::getline(fin,ln)) {
      size_t a=ln.find_first_not_of(" \t\r\n");
      if (a==std::string::npos) continue;
      std::string p=ln.substr(a,ln.find_last_not_of(" \t\r\n")-a+1);
      if (p.empty()||p[0]=='#') continue;
      chain.Add(p.c_str());
    }
  }
  std::cout<<"[lad_tof_fast] entries: "<<chain.GetEntries()<<"\n";
  if (!chain.GetEntries()) { std::cerr<<"empty chain\n"; return; }

  // ---------------------------------------------------------------
  // 1b. Detect which tracking variants are present in the data.
  //     Each variant t adds branches X.ladhod.goodhit_chiSquare<tsuf>.  A variant
  //     is used only if present for BOTH spectrometers (so every alias below
  //     resolves); the standard tracking is always assumed present.
  // ---------------------------------------------------------------
  chain.LoadTree(0);
  auto has_branch=[&chain](const std::string& n){ return chain.GetBranch(n.c_str())!=nullptr; };
  std::vector<Track> tracks;
  tracks.push_back({"standard",""});
  {
    const std::array<Track,3> variants={{
      {"xz","_xz"},{"noTrackVertex","_noTrackVertex"},{"noTrackVertex_xz","_noTrackVertex_xz"}}};
    for (const auto& v:variants) {
      bool ok=true;
      for (int is=0;is<N_SPECS;++is)
        if(!has_branch(std::string(1,specs[is])+".ladhod.goodhit_chiSquare"+v.tsuf)) ok=false;
      if(ok) tracks.push_back(v);
      else    std::cout<<"[lad_tof_fast] tracking variant '"<<v.dir<<"' absent in data; skipping\n";
    }
  }
  const int ntracks=(int)tracks.size();
  std::cout<<"[lad_tof_fast] tracking variants ("<<ntracks<<"):";
  for(const auto&t:tracks) std::cout<<" "<<t.dir;
  std::cout<<"\n";

  // ---------------------------------------------------------------
  // 2. RDataFrame + aliases
  // ---------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df=rdf;

  ROOT::RDF::RResultPtr<ULong64_t> progress_count;
#ifdef LAD_HAS_RDF_PROGRESSBAR
  ROOT::RDF::Experimental::AddProgressBar(df);
#else
  {
    const ULong64_t total=chain.GetEntries();
    const ULong64_t step=std::max<ULong64_t>(1ULL,total/200ULL);
    auto t0=std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());
    progress_count=df.Count();
    progress_count.OnPartialResult(step,[total,t0](ULong64_t n){
      double f=total?double(n)/double(total):0.;
      double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-*t0).count();
      std::fprintf(stderr,"\r[lad_tof_fast] %6.2f%% (%llu/%llu) %.1fs eta %.1fs   ",
                   100.*f,(unsigned long long)n,(unsigned long long)total,s,f>0?s*(1./f-1.):0.);
      std::fflush(stderr);
    });
  }
#endif

  for (int is=0;is<N_SPECS;++is) {
    const std::string sp(1,specs[is]);
    const std::string pfx=sp+".ladhod.goodhit_";
    for (const std::string &s:{"0","1"}) {
      df=df.Alias(sp+"_plane_"   +s, pfx+"plane_"      +s);
      df=df.Alias(sp+"_paddle_"  +s, pfx+"paddle_"     +s);
      df=df.Alias(sp+"_hittime_" +s, pfx+"hittime_"    +s);
      df=df.Alias(sp+"_edep_"    +s, pfx+"hitedep_"    +s);
      df=df.Alias(sp+"_edep_amp_"+s, pfx+"hitedep_amp_"+s);
      df=df.Alias(sp+"_tof_"     +s, pfx+"hit_tof_"    +s);
      df=df.Alias(sp+"_ypos_"    +s, pfx+"hit_ypos_"   +s);
    }
    for (const auto& tk:tracks)
      df=df.Alias(sp+"_chiSquare"+tk.tsuf, pfx+"chiSquare"+tk.tsuf);
    df=df.Alias(sp+"_isProton_0",  pfx+"isProton_0");
    df=df.Alias(sp+"_isProton_1",  pfx+"isProton_1");
  }

  // ---------------------------------------------------------------
  // 3. Categories and variables
  // ---------------------------------------------------------------
  // One "all_hits" category (tracking-independent) plus has_track/no_track for
  // each detected tracking variant, each living in its own <variant>/ sub-dir.
  // For the standard variant (tsuf="") the column suffixes are "_has"/"_no",
  // matching the original single-tracking column names.
  std::vector<Cat> cats;
  cats.push_back({"", "all_hits", 0, ""});
  for (const auto& tk:tracks) {
    cats.push_back({tk.tsuf+"_has", tk.dir+"/has_track", 1, tk.tsuf});
    cats.push_back({tk.tsuf+"_no",  tk.dir+"/no_track",  2, tk.tsuf});
  }
  const int ncats=(int)cats.size();

  const std::vector<VarConfig> vars={
    {"tof",     NBINS_TOF,     XMIN_TOF,     XMAX_TOF},
    {"ypos",    NBINS_YPOS,    XMIN_YPOS,    XMAX_YPOS},
    {"hittime", NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME},
    {"edep",    NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP},
    {"edep_amp",NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP},
  };

  auto chi_jit=[](const std::string &chi, int mode)->std::string{
    if (mode==1) return " && "+chi+" < 100";
    if (mode==2) return " && "+chi+" >= 100";
    return "";
  };

  // ---------------------------------------------------------------
  // 4. Column definitions — loop over both spectrometers
  // ---------------------------------------------------------------
  for (int is=0;is<N_SPECS;++is) {
    const std::string sp(1,specs[is]);

    // ------------------------------------------------------------------
    // 4a. Per-plane/paddle var columns
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const std::string &cs=cat.suf;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      const std::string cx=chi_jit(chi,cat.chi_mode);
      for (int plane=0;plane<N_PLANES;++plane) {
        const std::string side=(plane%2==0)?"0":"1";
        const std::string pc=sp+"_plane_"+side, pdc=sp+"_paddle_"+side;
        for (int paddle=0;paddle<N_PADDLES;++paddle) {
          const std::string mask="("+pc+"=="+std::to_string(plane)+" && "+pdc+"=="+std::to_string(paddle)+cx+")";
          for (const auto &v:vars)
            df=df.Define(sp+"_"+v.name+"_p"+std::to_string(plane)+"_b"+std::to_string(paddle)+cs,
                         sp+"_"+v.name+"_"+side+"["+mask+"]");
        }
        // Per-plane sum: exclude paddles 1 and 9 for planes 100 (idx 2) and 101 (idx 3)
        const std::string excl=(plane==2||plane==3)?" && "+pdc+"!=1 && "+pdc+"!=9":"";
        const std::string mpl="("+pc+"=="+std::to_string(plane)+cx+excl+")";
        for (const auto &v:vars)
          df=df.Define(sp+"_"+v.name+"_p"+std::to_string(plane)+"_sum"+cs,
                       sp+"_"+v.name+"_"+side+"["+mpl+"]");
      }
      // _total: concatenate per-plane sums so paddle exclusions are inherited
      for (const auto &v:vars) {
        std::vector<std::string> srcs;
        for(int p=0;p<N_PLANES;++p)
          srcs.push_back(sp+"_"+v.name+"_p"+std::to_string(p)+"_sum"+cs);
        df=df.Define(sp+"_"+v.name+"_total"+cs,
          [](const ROOT::VecOps::RVec<double>&v0,const ROOT::VecOps::RVec<double>&v1,
             const ROOT::VecOps::RVec<double>&v2,const ROOT::VecOps::RVec<double>&v3,
             const ROOT::VecOps::RVec<double>&v4){
            auto r=ROOT::VecOps::Concatenate(v0,v1); r=ROOT::VecOps::Concatenate(r,v2);
            r=ROOT::VecOps::Concatenate(r,v3); return ROOT::VecOps::Concatenate(r,v4);
          },srcs);
      }
    }

    // ------------------------------------------------------------------
    // 4b. edep vs punchthrough — split into 000-001 (_01) and 100-101 (_23) pairs
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const int cm=cat.chi_mode; const std::string &cs=cat.suf;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      for (int paddle=0;paddle<N_PADDLES;++paddle) {
        const double pv=paddle; const std::string ps=std::to_string(paddle);
        using RVd=ROOT::VecOps::RVec<double>;
        using Rpd=ROOT::VecOps::RVec<std::pair<double,double>>;

        // --- 000-001 pair (pl0==0, pl1==1) ---
        df=df.Define(sp+"_edep_0_vs_pt_01_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd e0,RVd pl1,RVd pd1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=0.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=1.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({e0[i],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_edep_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edep_1_vs_pt_01_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd pl1,RVd pd1,RVd e1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=0.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=1.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({e1[j],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",sp+"_edep_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edepamp_0_vs_pt_01_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd ea,RVd pl1,RVd pd1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=0.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=1.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({ea[i],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_edep_amp_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edepamp_1_vs_pt_01_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd pl1,RVd pd1,RVd ea,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=0.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=1.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({ea[j],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",sp+"_edep_amp_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        // --- 100-101 pair (pl0==2, pl1==3) ---
        df=df.Define(sp+"_edep_0_vs_pt_23_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd e0,RVd pl1,RVd pd1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=2.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=3.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({e0[i],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_edep_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edep_1_vs_pt_23_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd pl1,RVd pd1,RVd e1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=2.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=3.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({e1[j],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",sp+"_edep_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edepamp_0_vs_pt_23_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd ea,RVd pl1,RVd pd1,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=2.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=3.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({ea[i],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_edep_amp_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        df=df.Define(sp+"_edepamp_1_vs_pt_23_b"+ps+cs,
          [pv,cm](RVd pl0,RVd pd0,RVd pl1,RVd pd1,RVd ea,RVd h0,RVd h1,RVd chi){
            Rpd r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=2.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=3.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back({ea[j],h1[j]-h0[i]});
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",sp+"_edep_amp_1",sp+"_hittime_0",sp+"_hittime_1",chi});

        // Split pair columns into x(PT) and y(edep) for Histo2D
        auto sx=[](const Rpd&v){RVd x; for(auto&p:v)x.push_back(p.second); return x;};
        auto sy=[](const Rpd&v){RVd y; for(auto&p:v)y.push_back(p.first);  return y;};
        for(const std::string &pp:{"_01","_23"}){
          df=df.Define(sp+"_pt_0_b"     +pp+"_"+ps+cs,sx,{sp+"_edep_0_vs_pt"   +pp+"_b"+ps+cs});
          df=df.Define(sp+"_edep_0_b_pt"+pp+"_"+ps+cs,sy,{sp+"_edep_0_vs_pt"   +pp+"_b"+ps+cs});
          df=df.Define(sp+"_pt_1_b"     +pp+"_"+ps+cs,sx,{sp+"_edep_1_vs_pt"   +pp+"_b"+ps+cs});
          df=df.Define(sp+"_edep_1_b_pt"+pp+"_"+ps+cs,sy,{sp+"_edep_1_vs_pt"   +pp+"_b"+ps+cs});
          df=df.Define(sp+"_ptamp_0_b"      +pp+"_"+ps+cs,sx,{sp+"_edepamp_0_vs_pt"+pp+"_b"+ps+cs});
          df=df.Define(sp+"_edepamp_0_b_pt" +pp+"_"+ps+cs,sy,{sp+"_edepamp_0_vs_pt"+pp+"_b"+ps+cs});
          df=df.Define(sp+"_ptamp_1_b"      +pp+"_"+ps+cs,sx,{sp+"_edepamp_1_vs_pt"+pp+"_b"+ps+cs});
          df=df.Define(sp+"_edepamp_1_b_pt" +pp+"_"+ps+cs,sy,{sp+"_edepamp_1_vs_pt"+pp+"_b"+ps+cs});
        }
      }
    }

    // ------------------------------------------------------------------
    // 4c. tof vs edep 2D pair columns (per-paddle, no plane filtering needed)
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const int cm=cat.chi_mode; const std::string &cs=cat.suf;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      for (int paddle=0;paddle<N_PADDLES;++paddle) {
        const double pv=paddle; const std::string ps=std::to_string(paddle);

        auto mk_te=[&](const std::string &col,const std::string &pd_src,
                        const std::string &tof_src,const std::string &edep_src){
          df=df.Define(col,[pv,cm](ROOT::VecOps::RVec<double> pd,ROOT::VecOps::RVec<double> t,
                                    ROOT::VecOps::RVec<double> e, ROOT::VecOps::RVec<double> chi){
            ROOT::VecOps::RVec<std::pair<double,double>> r;
            for(size_t i=0;i<pd.size();++i){
              if(pd[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              r.push_back({t[i],e[i]});
            } return r;
          },{pd_src,tof_src,edep_src,chi});
        };
        mk_te(sp+"_tof_edep_0_b"   +ps+cs,sp+"_paddle_0",sp+"_tof_0",sp+"_edep_0");
        mk_te(sp+"_tof_edep_1_b"   +ps+cs,sp+"_paddle_1",sp+"_tof_1",sp+"_edep_1");
        mk_te(sp+"_tof_edepamp_0_b"+ps+cs,sp+"_paddle_0",sp+"_tof_0",sp+"_edep_amp_0");
        mk_te(sp+"_tof_edepamp_1_b"+ps+cs,sp+"_paddle_1",sp+"_tof_1",sp+"_edep_amp_1");

        auto sx=[](const ROOT::VecOps::RVec<std::pair<double,double>>&v){
          ROOT::VecOps::RVec<double> x; for(auto&p:v)x.push_back(p.first);  return x;};
        auto sy=[](const ROOT::VecOps::RVec<std::pair<double,double>>&v){
          ROOT::VecOps::RVec<double> y; for(auto&p:v)y.push_back(p.second); return y;};
        df=df.Define(sp+"_tof_0_b"        +ps+cs,sx,{sp+"_tof_edep_0_b"   +ps+cs});
        df=df.Define(sp+"_edep_0_b_tof"   +ps+cs,sy,{sp+"_tof_edep_0_b"   +ps+cs});
        df=df.Define(sp+"_tof_1_b"        +ps+cs,sx,{sp+"_tof_edep_1_b"   +ps+cs});
        df=df.Define(sp+"_edep_1_b_tof"   +ps+cs,sy,{sp+"_tof_edep_1_b"   +ps+cs});
        df=df.Define(sp+"_tofamp_0_b"     +ps+cs,sx,{sp+"_tof_edepamp_0_b"+ps+cs});
        df=df.Define(sp+"_edepamp_0_b_tof"+ps+cs,sy,{sp+"_tof_edepamp_0_b"+ps+cs});
        df=df.Define(sp+"_tofamp_1_b"     +ps+cs,sx,{sp+"_tof_edepamp_1_b"+ps+cs});
        df=df.Define(sp+"_edepamp_1_b_tof"+ps+cs,sy,{sp+"_tof_edepamp_1_b"+ps+cs});
      }
    }

    // ------------------------------------------------------------------
    // 4d. Punchthrough per-paddle and sums
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const int cm=cat.chi_mode; const std::string &cs=cat.suf;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      for (int paddle=0;paddle<N_PADDLES;++paddle) {
        const double pv=paddle; const std::string ps=std::to_string(paddle);

        df=df.Define(sp+"_punchthrough_01_b"+ps+cs,
          [pv,cm](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,ROOT::VecOps::RVec<double> h0,
                  ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,ROOT::VecOps::RVec<double> h1,
                  ROOT::VecOps::RVec<double> chi){
            ROOT::VecOps::RVec<double> r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=0.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=1.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back(h1[j]-h0[i]);
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_hittime_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_1",chi});

        df=df.Define(sp+"_punchthrough_23_b"+ps+cs,
          [pv,cm](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,ROOT::VecOps::RVec<double> h0,
                  ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,ROOT::VecOps::RVec<double> h1,
                  ROOT::VecOps::RVec<double> chi){
            ROOT::VecOps::RVec<double> r;
            for(size_t i=0;i<pl0.size();++i){
              if(pl0[i]!=2.||pd0[i]!=pv) continue;
              if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=3.||pd1[j]!=pv) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                r.push_back(h1[j]-h0[i]);
              }
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_hittime_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_1",chi});
      }

      df=df.Define(sp+"_punchthrough_01_sum"+cs,
        [cm](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,ROOT::VecOps::RVec<double> h0,
             ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,ROOT::VecOps::RVec<double> h1,
             ROOT::VecOps::RVec<double> chi){
          ROOT::VecOps::RVec<double> r;
          for(size_t i=0;i<pl0.size();++i){
            if(pl0[i]!=0.) continue;
            if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
            for(size_t j=0;j<pl1.size();++j){
              if(pl1[j]!=1.||pd1[j]!=pd0[i]) continue;
              if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
              r.push_back(h1[j]-h0[i]);
            }
          } return r;
        },{sp+"_plane_0",sp+"_paddle_0",sp+"_hittime_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_1",chi});

      // 23 sum: exclude paddles 1 and 9 (planes 100/101)
      df=df.Define(sp+"_punchthrough_23_sum"+cs,
        [cm](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,ROOT::VecOps::RVec<double> h0,
             ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,ROOT::VecOps::RVec<double> h1,
             ROOT::VecOps::RVec<double> chi){
          ROOT::VecOps::RVec<double> r;
          for(size_t i=0;i<pl0.size();++i){
            if(pl0[i]!=2.) continue;
            if(pd0[i]==1.||pd0[i]==9.) continue; // exclude bad paddles
            if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
            for(size_t j=0;j<pl1.size();++j){
              if(pl1[j]!=3.||pd1[j]!=pd0[i]) continue;
              if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
              r.push_back(h1[j]-h0[i]);
            }
          } return r;
        },{sp+"_plane_0",sp+"_paddle_0",sp+"_hittime_0",sp+"_plane_1",sp+"_paddle_1",sp+"_hittime_1",chi});

      df=df.Define(sp+"_punchthrough_sum_total"+cs,
        [](const ROOT::VecOps::RVec<double>&a,const ROOT::VecOps::RVec<double>&b){
          return ROOT::VecOps::Concatenate(a,b);
        },{sp+"_punchthrough_01_sum"+cs,sp+"_punchthrough_23_sum"+cs});
    }

    // ------------------------------------------------------------------
    // 4e. Front veto per-paddle and sums
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const int cm=cat.chi_mode; const std::string &cs=cat.suf;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      for (int paddle=0;paddle<N_PADDLES;++paddle) {
        const double pav=static_cast<double>(paddle); const std::string ps=std::to_string(paddle);

        auto mk_fv=[&](const std::string &col,double bk,double fr,const std::string &val_src){
          df=df.Define(col,
            [bk,fr,pav,cm](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,
                            ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,
                            ROOT::VecOps::RVec<double> val,ROOT::VecOps::RVec<double> chi){
              ROOT::VecOps::RVec<double> r;
              for(size_t j=0;j<pl1.size();++j){
                if(pl1[j]!=bk||pd1[j]!=pav) continue;
                if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
                bool veto=false;
                for(size_t i=0;i<pl0.size();++i){
                  if(pl0[i]!=fr||pd0[i]!=pav) continue;
                  if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
                  veto=true; break;
                }
                if(!veto) r.push_back(val[j]);
              } return r;
            },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",val_src,chi});
        };
        mk_fv(sp+"_hittime_fv_01_b"+ps+cs,1.,0.,sp+"_hittime_1");
        mk_fv(sp+"_tof_fv_01_b"    +ps+cs,1.,0.,sp+"_tof_1");
        mk_fv(sp+"_hittime_fv_23_b"+ps+cs,3.,2.,sp+"_hittime_1");
        mk_fv(sp+"_tof_fv_23_b"    +ps+cs,3.,2.,sp+"_tof_1");
      }

      // Sum lambdas: for 23 pair also exclude paddles 1 and 9
      auto mk_fv_sum=[&](const std::string &col,double bk,double fr,
                          const std::string &val_src,bool excl_pad){
        df=df.Define(col,
          [bk,fr,cm,excl_pad](ROOT::VecOps::RVec<double> pl0,ROOT::VecOps::RVec<double> pd0,
                               ROOT::VecOps::RVec<double> pl1,ROOT::VecOps::RVec<double> pd1,
                               ROOT::VecOps::RVec<double> val,ROOT::VecOps::RVec<double> chi){
            ROOT::VecOps::RVec<double> r;
            for(size_t j=0;j<pl1.size();++j){
              if(pl1[j]!=bk) continue;
              if(excl_pad&&(pd1[j]==1.||pd1[j]==9.)) continue;
              if(cm==1&&chi[j]>=100.) continue; if(cm==2&&chi[j]<100.) continue;
              double pa=pd1[j]; bool veto=false;
              for(size_t i=0;i<pl0.size();++i){
                if(pl0[i]!=fr||pd0[i]!=pa) continue;
                if(cm==1&&chi[i]>=100.) continue; if(cm==2&&chi[i]<100.) continue;
                veto=true; break;
              }
              if(!veto) r.push_back(val[j]);
            } return r;
          },{sp+"_plane_0",sp+"_paddle_0",sp+"_plane_1",sp+"_paddle_1",val_src,chi});
      };
      mk_fv_sum(sp+"_hittime_fv_01_sum"+cs,1.,0.,sp+"_hittime_1",false);
      mk_fv_sum(sp+"_tof_fv_01_sum"    +cs,1.,0.,sp+"_tof_1",    false);
      mk_fv_sum(sp+"_hittime_fv_23_sum"+cs,3.,2.,sp+"_hittime_1", true);
      mk_fv_sum(sp+"_tof_fv_23_sum"    +cs,3.,2.,sp+"_tof_1",     true);

      df=df.Define(sp+"_hittime_fv_total"+cs,
        [](const ROOT::VecOps::RVec<double>&a,const ROOT::VecOps::RVec<double>&b){
          return ROOT::VecOps::Concatenate(a,b);
        },{sp+"_hittime_fv_01_sum"+cs,sp+"_hittime_fv_23_sum"+cs});
      df=df.Define(sp+"_tof_fv_total"+cs,
        [](const ROOT::VecOps::RVec<double>&a,const ROOT::VecOps::RVec<double>&b){
          return ROOT::VecOps::Concatenate(a,b);
        },{sp+"_tof_fv_01_sum"+cs,sp+"_tof_fv_23_sum"+cs});
    }

    // ------------------------------------------------------------------
    // 4f. Path-length corrected tof (per-paddle and per-plane sum)
    //     Sum for planes 2/3 excludes paddles 1 and 9
    // ------------------------------------------------------------------
    for (const auto &cat:cats) {
      const std::string &cs=cat.suf;
      const int cm2=cat.chi_mode;
      const std::string chi=sp+"_chiSquare"+cat.tsuf;
      for (int plane=0;plane<N_PLANES;++plane) {
        const std::string side=(plane%2==0)?"0":"1";
        const double pl_val=static_cast<double>(plane);
        const double R=hodo_radii[plane];
        const bool excl_pad=(plane==2||plane==3);
        for (int paddle=0;paddle<N_PADDLES;++paddle) {
          const double dx=22.*(static_cast<double>(paddle)-6.);
          const std::string tc=sp+"_tof_p" +std::to_string(plane)+"_b"+std::to_string(paddle)+cs;
          const std::string yc=sp+"_ypos_p"+std::to_string(plane)+"_b"+std::to_string(paddle)+cs;
          df=df.Define(sp+"_tof_corrected_p"+std::to_string(plane)+"_b"+std::to_string(paddle)+cs,
            [dx,R](const ROOT::VecOps::RVec<double>&tv,const ROOT::VecOps::RVec<double>&yv){
              ROOT::VecOps::RVec<double> r; r.reserve(tv.size());
              for(size_t i=0;i<tv.size();++i){double p2d=std::sqrt(yv[i]*yv[i]+dx*dx);
                r.push_back(tv[i]-std::sqrt(p2d*p2d+R*R)/100./0.3);} return r;
            },{tc,yc});
        }
        df=df.Define(sp+"_tof_corrected_p"+std::to_string(plane)+"_sum"+cs,
          [pl_val,R,cm2,excl_pad](const ROOT::VecOps::RVec<double>&plv,const ROOT::VecOps::RVec<double>&pdv,
                                   const ROOT::VecOps::RVec<double>&yv, const ROOT::VecOps::RVec<double>&tv,
                                   const ROOT::VecOps::RVec<double>&chi){
            ROOT::VecOps::RVec<double> r;
            for(size_t i=0;i<plv.size();++i){
              if(plv[i]!=pl_val) continue;
              if(excl_pad&&(pdv[i]==1.||pdv[i]==9.)) continue;
              if(cm2==1&&chi[i]>=100.) continue; if(cm2==2&&chi[i]<100.) continue;
              double dx=22.*(pdv[i]-6.); double p2d=std::sqrt(yv[i]*yv[i]+dx*dx);
              r.push_back(tv[i]-std::sqrt(p2d*p2d+R*R)/100./0.3);} return r;
          },{sp+"_plane_"+side,sp+"_paddle_"+side,sp+"_ypos_"+side,sp+"_tof_"+side,chi});
      }

      auto cat5=[&](const std::string &pfx)->std::vector<std::string>{
        std::vector<std::string> v;
        for(int p=0;p<N_PLANES;++p) v.push_back(pfx+std::to_string(p)+"_sum"+cs);
        return v;};
      df=df.Define(sp+"_tof_corrected_total"+cs,
        [](const ROOT::VecOps::RVec<double>&v0,const ROOT::VecOps::RVec<double>&v1,
           const ROOT::VecOps::RVec<double>&v2,const ROOT::VecOps::RVec<double>&v3,
           const ROOT::VecOps::RVec<double>&v4){
          auto r=ROOT::VecOps::Concatenate(v0,v1); r=ROOT::VecOps::Concatenate(r,v2);
          r=ROOT::VecOps::Concatenate(r,v3); return ROOT::VecOps::Concatenate(r,v4);
        },cat5(sp+"_tof_corrected_p"));
    }
    // ------------------------------------------------------------------
    // 4g. Proton-tagged corrected tof (side-1 hits with isProton_1==1)
    //     Planes 001 (idx 1) and 101 (idx 3) only; radius looked up from hodo_radii.
    // ------------------------------------------------------------------
    {
      using RVd=ROOT::VecOps::RVec<double>;
      auto mk_proton=[&](const std::string& col, bool req_track, const std::string& chicol){
        df=df.Define(col,
          [req_track](const RVd& pl1,const RVd& pd1,const RVd& y1,const RVd& t1,
                      const RVd& ip1,const RVd& chi){
            RVd r;
            for(size_t i=0;i<pl1.size();++i){
              if(ip1[i]!=1.) continue;
              if(req_track&&chi[i]>=100.) continue;
              int pi=(int)std::round(pl1[i]);
              if(pi!=1&&pi!=3) continue; // only planes 001 and 101
              double R=hodo_radii[pi];
              double dx=22.*(pd1[i]-6.);
              double p2d=std::sqrt(y1[i]*y1[i]+dx*dx);
              r.push_back(t1[i]-std::sqrt(p2d*p2d+R*R)/100./0.3);
            } return r;
          },{sp+"_plane_1",sp+"_paddle_1",sp+"_ypos_1",sp+"_tof_1",sp+"_isProton_1",chicol});
      };
      // proton-only is tracking-independent; proton+track once per tracking variant.
      mk_proton(sp+"_tof_corr_proton", false, sp+"_chiSquare");
      for(const auto& tk:tracks)
        mk_proton(sp+"_tof_corr_proton_track"+tk.tsuf, true, sp+"_chiSquare"+tk.tsuf);

      // Per-paddle and per-plane sums for planes 001 (idx 1) and 101 (idx 3)
      for(int pp=0;pp<2;++pp){
        const int   pi_c =(pp==0)?1:3;
        const double R_c  =hodo_radii[pi_c];
        const bool excl_c=(pi_c==3); // exclude paddles 1,9 for plane 101

        for(int paddle=0;paddle<N_PADDLES;++paddle){
          const double pv_c=static_cast<double>(paddle);
          const std::string ps=std::to_string(paddle);
          auto mk_p=[&](const std::string& col,bool req_track,const std::string& chicol){
            df=df.Define(col,
              [pi_c,R_c,pv_c,req_track](const RVd& pl1,const RVd& pd1,const RVd& y1,
                                         const RVd& t1, const RVd& ip1,const RVd& chi){
                RVd r;
                for(size_t i=0;i<pl1.size();++i){
                  if(ip1[i]!=1.) continue;
                  if((int)std::round(pl1[i])!=pi_c) continue;
                  if(pd1[i]!=pv_c) continue;
                  if(req_track&&chi[i]>=100.) continue;
                  double dx=22.*(pd1[i]-6.);
                  double p2d=std::sqrt(y1[i]*y1[i]+dx*dx);
                  r.push_back(t1[i]-std::sqrt(p2d*p2d+R_c*R_c)/100./0.3);
                } return r;
              },{sp+"_plane_1",sp+"_paddle_1",sp+"_ypos_1",sp+"_tof_1",sp+"_isProton_1",chicol});
          };
          mk_p(sp+"_tof_corr_proton_p"+std::to_string(pi_c)+"_b"+ps, false, sp+"_chiSquare");
          for(const auto& tk:tracks)
            mk_p(sp+"_tof_corr_proton_track_p"+std::to_string(pi_c)+"_b"+ps+tk.tsuf, true, sp+"_chiSquare"+tk.tsuf);
        }

        // Per-plane sum (plane 101 excludes paddles 1 and 9)
        auto mk_sum=[&](const std::string& col,bool req_track,const std::string& chicol){
          df=df.Define(col,
            [pi_c,R_c,excl_c,req_track](const RVd& pl1,const RVd& pd1,const RVd& y1,
                                         const RVd& t1, const RVd& ip1,const RVd& chi){
              RVd r;
              for(size_t i=0;i<pl1.size();++i){
                if(ip1[i]!=1.) continue;
                if((int)std::round(pl1[i])!=pi_c) continue;
                if(excl_c&&(pd1[i]==1.||pd1[i]==9.)) continue;
                if(req_track&&chi[i]>=100.) continue;
                double dx=22.*(pd1[i]-6.);
                double p2d=std::sqrt(y1[i]*y1[i]+dx*dx);
                r.push_back(t1[i]-std::sqrt(p2d*p2d+R_c*R_c)/100./0.3);
              } return r;
            },{sp+"_plane_1",sp+"_paddle_1",sp+"_ypos_1",sp+"_tof_1",sp+"_isProton_1",chicol});
        };
        mk_sum(sp+"_tof_corr_proton_p"+std::to_string(pi_c)+"_sum", false, sp+"_chiSquare");
        for(const auto& tk:tracks)
          mk_sum(sp+"_tof_corr_proton_track_p"+std::to_string(pi_c)+"_sum"+tk.tsuf, true, sp+"_chiSquare"+tk.tsuf);
      }
    }
  } // end spec loop (column definitions)

  // ===============================================================
  // 5. Histogram booking
  // ===============================================================
  using RH1=ROOT::RDF::RResultPtr<::TH1D>;
  using RH2=ROOT::RDF::RResultPtr<::TH2D>;

  std::array<std::array<std::map<std::string,std::vector<std::vector<RH1>>>,N_CATS>,N_SPECS> hpad;
  std::array<std::array<std::map<std::string,std::vector<RH1>>,N_CATS>,N_SPECS>              hsum;
  std::array<std::array<std::map<std::string,RH1>,N_CATS>,N_SPECS>                           htot;
  std::array<std::array<std::vector<RH1>,N_CATS>,N_SPECS> h_pt_01_pad,h_pt_23_pad;
  std::array<std::array<RH1,N_CATS>,N_SPECS>              h_pt_01_sum,h_pt_23_sum,h_pt_tot;
  std::array<std::array<std::vector<RH1>,N_CATS>,N_SPECS> h_ht_fv_01,h_ht_fv_23,h_tof_fv_01,h_tof_fv_23;
  std::array<std::array<RH1,N_CATS>,N_SPECS> h_ht_fv_01_sum,h_ht_fv_23_sum,h_ht_fv_tot;
  std::array<std::array<RH1,N_CATS>,N_SPECS> h_tof_fv_01_sum,h_tof_fv_23_sum,h_tof_fv_tot;
  std::array<std::array<std::vector<std::vector<RH1>>,N_CATS>,N_SPECS> h_tof_corr;
  std::array<std::array<std::vector<RH1>,N_CATS>,N_SPECS>              h_tof_corr_sum;
  std::array<std::array<RH1,N_CATS>,N_SPECS>                           h_tof_corr_tot;
  // edep vs PT split by plane pair
  std::array<std::array<std::map<std::string,std::vector<RH2>>,N_CATS>,N_SPECS> h_edep_vs_pt_01,h_edep_vs_pt_23;
  std::array<std::array<std::map<std::string,std::vector<RH2>>,N_CATS>,N_SPECS> h_tof_vs_edep;
  // Proton-tagged corrected tof (per spec, no category split — cuts are built-in).
  // Proton-only histos have no tracking dimension; proton+track histos carry an
  // extra [track] index (one booking per detected tracking variant).
  std::array<RH1,N_SPECS>                      h_proton_tof;
  std::array<std::array<RH1,N_TRACKS>,N_SPECS> h_proton_track_tof;
  // Per-paddle [spec][pp (0=001,1=101)][paddle] and per-plane sum [spec][pp];
  // track variants insert a [track] index after [spec].
  std::array<std::array<std::vector<RH1>,2>,N_SPECS>                       h_proton_pad;
  std::array<std::array<std::array<std::vector<RH1>,2>,N_TRACKS>,N_SPECS>  h_proton_track_pad;
  std::array<std::array<RH1,2>,N_SPECS>                                    h_proton_psum;
  std::array<std::array<std::array<RH1,2>,N_TRACKS>,N_SPECS>               h_proton_track_psum;

  for (int is=0;is<N_SPECS;++is) {
    const std::string sp(1,specs[is]);
    for (int ic=0;ic<ncats;++ic) {
      const std::string &cs=cats[ic].suf;
      auto bk1=[&](const std::string &col,const std::string &ttl,int nb,double lo,double hi)->RH1{
        return df.Histo1D({col.c_str(),ttl.c_str(),nb,lo,hi},col);};

      // Vars
      for(const auto&v:vars){
        hpad[is][ic][v.name].resize(N_PLANES);
        hsum[is][ic][v.name].resize(N_PLANES);
        for(int p=0;p<N_PLANES;++p) hpad[is][ic][v.name][p].resize(N_PADDLES);
      }
      for(const auto&v:vars){
        for(int pl=0;pl<N_PLANES;++pl){
          for(int pa=0;pa<N_PADDLES;++pa){
            const std::string c=sp+"_"+v.name+"_p"+std::to_string(pl)+"_b"+std::to_string(pa)+cs;
            const std::string t=sp+" "+v.name+" pl"+plane_names[pl]+" pd"+std::to_string(pa)+";"+v.name+";Counts";
            hpad[is][ic][v.name][pl][pa]=bk1(c,t,v.nbins,v.xmin,v.xmax);
          }
          const std::string sc=sp+"_"+v.name+"_p"+std::to_string(pl)+"_sum"+cs;
          const std::string st=sp+" "+v.name+" plane "+plane_names[pl]+";"+v.name+";Counts";
          hsum[is][ic][v.name][pl]=bk1(sc,st,v.nbins,v.xmin,v.xmax);
        }
        const std::string tc=sp+"_"+v.name+"_total"+cs;
        htot[is][ic][v.name]=bk1(tc,sp+" "+v.name+" all planes;"+v.name+";Counts",v.nbins,v.xmin,v.xmax);
      }

      // Punchthrough
      h_pt_01_pad[is][ic].resize(N_PADDLES); h_pt_23_pad[is][ic].resize(N_PADDLES);
      for(int pa=0;pa<N_PADDLES;++pa){
        const std::string ps=std::to_string(pa);
        h_pt_01_pad[is][ic][pa]=bk1(sp+"_punchthrough_01_b"+ps+cs,
          sp+" PT 000-001 pd"+ps+";PT(ns);Counts",NBINS_PT,XMIN_PT,XMAX_PT);
        h_pt_23_pad[is][ic][pa]=bk1(sp+"_punchthrough_23_b"+ps+cs,
          sp+" PT 100-101 pd"+ps+";PT(ns);Counts",NBINS_PT,XMIN_PT,XMAX_PT);
      }
      h_pt_01_sum[is][ic]=bk1(sp+"_punchthrough_01_sum"+cs,sp+" PT 000-001 sum;PT(ns);Counts",NBINS_PT,XMIN_PT,XMAX_PT);
      h_pt_23_sum[is][ic]=bk1(sp+"_punchthrough_23_sum"+cs,sp+" PT 100-101 sum;PT(ns);Counts",NBINS_PT,XMIN_PT,XMAX_PT);
      h_pt_tot[is][ic]   =bk1(sp+"_punchthrough_sum_total"+cs,sp+" PT total;PT(ns);Counts",    NBINS_PT,XMIN_PT,XMAX_PT);

      // Front veto
      h_ht_fv_01[is][ic].resize(N_PADDLES);  h_ht_fv_23[is][ic].resize(N_PADDLES);
      h_tof_fv_01[is][ic].resize(N_PADDLES); h_tof_fv_23[is][ic].resize(N_PADDLES);
      for(int pa=0;pa<N_PADDLES;++pa){
        const std::string ps=std::to_string(pa);
        h_ht_fv_01[is][ic][pa] =bk1(sp+"_hittime_fv_01_b"+ps+cs,sp+" ht fv 001 pd"+ps+";ht(ns);Counts", NBINS_HITTIME,XMIN_HITTIME,XMAX_HITTIME);
        h_tof_fv_01[is][ic][pa]=bk1(sp+"_tof_fv_01_b"+ps+cs,    sp+" tof fv 001 pd"+ps+";tof(ns);Counts",NBINS_TOF,XMIN_TOF,XMAX_TOF);
        h_ht_fv_23[is][ic][pa] =bk1(sp+"_hittime_fv_23_b"+ps+cs,sp+" ht fv 101 pd"+ps+";ht(ns);Counts", NBINS_HITTIME,XMIN_HITTIME,XMAX_HITTIME);
        h_tof_fv_23[is][ic][pa]=bk1(sp+"_tof_fv_23_b"+ps+cs,    sp+" tof fv 101 pd"+ps+";tof(ns);Counts",NBINS_TOF,XMIN_TOF,XMAX_TOF);
      }
      h_ht_fv_01_sum[is][ic] =bk1(sp+"_hittime_fv_01_sum"+cs,sp+" ht fv 001 sum;ht(ns);Counts",  NBINS_HITTIME,XMIN_HITTIME,XMAX_HITTIME);
      h_tof_fv_01_sum[is][ic]=bk1(sp+"_tof_fv_01_sum"    +cs,sp+" tof fv 001 sum;tof(ns);Counts",NBINS_TOF,XMIN_TOF,XMAX_TOF);
      h_ht_fv_23_sum[is][ic] =bk1(sp+"_hittime_fv_23_sum"+cs,sp+" ht fv 101 sum;ht(ns);Counts",  NBINS_HITTIME,XMIN_HITTIME,XMAX_HITTIME);
      h_tof_fv_23_sum[is][ic]=bk1(sp+"_tof_fv_23_sum"    +cs,sp+" tof fv 101 sum;tof(ns);Counts",NBINS_TOF,XMIN_TOF,XMAX_TOF);
      h_ht_fv_tot[is][ic]    =bk1(sp+"_hittime_fv_total"  +cs,sp+" ht fv total;ht(ns);Counts",   NBINS_HITTIME,XMIN_HITTIME,XMAX_HITTIME);
      h_tof_fv_tot[is][ic]   =bk1(sp+"_tof_fv_total"      +cs,sp+" tof fv total;tof(ns);Counts", NBINS_TOF,XMIN_TOF,XMAX_TOF);

      // Corrected tof
      h_tof_corr[is][ic].assign(N_PLANES,std::vector<RH1>(N_PADDLES));
      h_tof_corr_sum[is][ic].resize(N_PLANES);
      for(int pl=0;pl<N_PLANES;++pl){
        for(int pa=0;pa<N_PADDLES;++pa){
          const std::string ctc=sp+"_tof_corrected_p"+std::to_string(pl)+"_b"+std::to_string(pa)+cs;
          h_tof_corr[is][ic][pl][pa]=bk1(ctc,
            sp+" tof corr "+plane_names[pl]+" pd"+std::to_string(pa)+";tof-L/c(ns);Counts",
            NBINS_TCORR,XMIN_TCORR,XMAX_TCORR);
        }
        const std::string cs2=sp+"_tof_corrected_p"+std::to_string(pl)+"_sum"+cs;
        h_tof_corr_sum[is][ic][pl]=bk1(cs2,
          sp+" tof corr "+plane_names[pl]+" sum;tof-L/c(ns);Counts",NBINS_TCORR,XMIN_TCORR,XMAX_TCORR);
      }
      h_tof_corr_tot[is][ic]=bk1(sp+"_tof_corrected_total"+cs,
        sp+" tof corr total;tof-L/c(ns);Counts",NBINS_TCORR,XMIN_TCORR,XMAX_TCORR);

      // 2D edep vs PT — separate for 000-001 and 100-101
      for(const std::string&k:{"edep_0","edep_1","edep_amp_0","edep_amp_1"}){
        h_edep_vs_pt_01[is][ic][k].resize(N_PADDLES);
        h_edep_vs_pt_23[is][ic][k].resize(N_PADDLES);
      }
      for(int pa=0;pa<N_PADDLES;++pa){
        const std::string ps=std::to_string(pa);
        auto bk2pt=[&](std::map<std::string,std::vector<RH2>> &hmap,
                        const std::string &pp, const std::string&key,
                        const std::string&xc,const std::string&yc,
                        int nx,double lx,double hx,int ny,double ly,double hy){
          const std::string hn=sp+"_h_"+key+"_vs_pt"+pp+"_"+ps+cs;
          const std::string pairname=(pp=="_01")?"000-001":"100-101";
          hmap[key][pa]=df.Histo2D(
            {hn.c_str(),(sp+" "+key+" vs PT "+pairname+" pd"+ps+";PT(ns);"+key).c_str(),nx,lx,hx,ny,ly,hy},xc,yc);};
        bk2pt(h_edep_vs_pt_01[is][ic],"_01","edep_0",    sp+"_pt_0_b_01_"   +ps+cs,sp+"_edep_0_b_pt_01_"   +ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2pt(h_edep_vs_pt_01[is][ic],"_01","edep_1",    sp+"_pt_1_b_01_"   +ps+cs,sp+"_edep_1_b_pt_01_"   +ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2pt(h_edep_vs_pt_01[is][ic],"_01","edep_amp_0",sp+"_ptamp_0_b_01_"+ps+cs,sp+"_edepamp_0_b_pt_01_"+ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
        bk2pt(h_edep_vs_pt_01[is][ic],"_01","edep_amp_1",sp+"_ptamp_1_b_01_"+ps+cs,sp+"_edepamp_1_b_pt_01_"+ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
        bk2pt(h_edep_vs_pt_23[is][ic],"_23","edep_0",    sp+"_pt_0_b_23_"   +ps+cs,sp+"_edep_0_b_pt_23_"   +ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2pt(h_edep_vs_pt_23[is][ic],"_23","edep_1",    sp+"_pt_1_b_23_"   +ps+cs,sp+"_edep_1_b_pt_23_"   +ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2pt(h_edep_vs_pt_23[is][ic],"_23","edep_amp_0",sp+"_ptamp_0_b_23_"+ps+cs,sp+"_edepamp_0_b_pt_23_"+ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
        bk2pt(h_edep_vs_pt_23[is][ic],"_23","edep_amp_1",sp+"_ptamp_1_b_23_"+ps+cs,sp+"_edepamp_1_b_pt_23_"+ps+cs,NBINS_PT,XMIN_PT,XMAX_PT,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
      }

      // 2D tof vs edep
      for(const std::string&k:{"edep_0","edep_1","edep_amp_0","edep_amp_1"})
        h_tof_vs_edep[is][ic][k].resize(N_PADDLES);
      for(int pa=0;pa<N_PADDLES;++pa){
        const std::string ps=std::to_string(pa);
        auto bk2te=[&](const std::string&key,const std::string&xc,const std::string&yc,
                        int nx,double lx,double hx,int ny,double ly,double hy){
          const std::string hn=sp+"_h_tof_"+key+"_"+ps+cs;
          h_tof_vs_edep[is][ic][key][pa]=df.Histo2D(
            {hn.c_str(),(sp+" tof vs "+key+" pd"+ps+";tof(ns);"+key).c_str(),nx,lx,hx,ny,ly,hy},xc,yc);};
        bk2te("edep_0",    sp+"_tof_0_b"   +ps+cs,sp+"_edep_0_b_tof"   +ps+cs,NBINS_TOF,XMIN_TOF,XMAX_TOF,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2te("edep_1",    sp+"_tof_1_b"   +ps+cs,sp+"_edep_1_b_tof"   +ps+cs,NBINS_TOF,XMIN_TOF,XMAX_TOF,NBINS_EDEP,    XMIN_EDEP,    XMAX_EDEP);
        bk2te("edep_amp_0",sp+"_tofamp_0_b"+ps+cs,sp+"_edepamp_0_b_tof"+ps+cs,NBINS_TOF,XMIN_TOF,XMAX_TOF,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
        bk2te("edep_amp_1",sp+"_tofamp_1_b"+ps+cs,sp+"_edepamp_1_b_tof"+ps+cs,NBINS_TOF,XMIN_TOF,XMAX_TOF,NBINS_EDEP_AMP,XMIN_EDEP_AMP,XMAX_EDEP_AMP);
      }
    }
  } // end spec loop (histogram booking)

  // Proton histogram booking (outside cat loop — no category split).
  // Proton-only booked once; proton+track booked once per tracking variant.
  for(int is=0;is<N_SPECS;++is){
    const std::string sp(1,specs[is]);
    auto bkp=[&](const std::string& col,const std::string& ttl)->RH1{
      return df.Histo1D({col.c_str(),ttl.c_str(),NBINS_TCORR,XMIN_TCORR,XMAX_TCORR},col);};
    h_proton_tof[is]=bkp(sp+"_tof_corr_proton", sp+" tof corr proton;tof-L/c(ns);Counts");
    for(int it=0;it<ntracks;++it){
      const std::string& td=tracks[it].dir; const std::string& tu=tracks[it].tsuf;
      h_proton_track_tof[is][it]=bkp(sp+"_tof_corr_proton_track"+tu,
        sp+" tof corr proton+track ["+td+"];tof-L/c(ns);Counts");
    }
    for(int pp=0;pp<2;++pp){
      const int pi=(pp==0)?1:3;
      const std::string pn=plane_names[pi];
      h_proton_pad[is][pp].resize(N_PADDLES);
      for(int it=0;it<ntracks;++it) h_proton_track_pad[is][it][pp].resize(N_PADDLES);
      for(int pa=0;pa<N_PADDLES;++pa){
        const std::string ps=std::to_string(pa);
        const std::string cp=sp+"_tof_corr_proton_p"+std::to_string(pi)+"_b"+ps;
        h_proton_pad[is][pp][pa]=bkp(cp,sp+" tof corr proton "+pn+" pd"+ps+";tof-L/c(ns);Counts");
        for(int it=0;it<ntracks;++it){
          const std::string& td=tracks[it].dir; const std::string& tu=tracks[it].tsuf;
          const std::string ct=sp+"_tof_corr_proton_track_p"+std::to_string(pi)+"_b"+ps+tu;
          h_proton_track_pad[is][it][pp][pa]=bkp(ct,
            sp+" tof corr proton+track ["+td+"] "+pn+" pd"+ps+";tof-L/c(ns);Counts");
        }
      }
      const std::string cs2=sp+"_tof_corr_proton_p"+std::to_string(pi)+"_sum";
      h_proton_psum[is][pp]=bkp(cs2,sp+" tof corr proton "+pn+" sum;tof-L/c(ns);Counts");
      for(int it=0;it<ntracks;++it){
        const std::string& td=tracks[it].dir; const std::string& tu=tracks[it].tsuf;
        const std::string ts2=sp+"_tof_corr_proton_track_p"+std::to_string(pi)+"_sum"+tu;
        h_proton_track_psum[is][it][pp]=bkp(ts2,
          sp+" tof corr proton+track ["+td+"] "+pn+" sum;tof-L/c(ns);Counts");
      }
    }
  }

  // ===============================================================
  // 6. Trigger event loop
  // ===============================================================
  std::cout<<"[lad_tof_fast] Running event loop...\n";
  (void)htot[0][0][vars.front().name]->GetEntries();
#ifndef LAD_HAS_RDF_PROGRESSBAR
  std::fprintf(stderr,"\n");
#endif

  // ===============================================================
  // 7. Write output
  // ===============================================================
  TFile fout(out_file,"RECREATE");
  if(fout.IsZombie()){std::cerr<<"cannot open output\n";return;}

  auto wc=[](TCanvas*c){c->Write();delete c;};

  // Create (or fetch) a possibly-nested "a/b/c" directory under base.
  auto mkpath=[](TDirectory* base,const std::string& path)->TDirectory*{
    TDirectory* dcur=base; std::string seg; std::istringstream ss(path);
    while(std::getline(ss,seg,'/')){
      if(seg.empty()) continue;
      TDirectory* nd=dcur->GetDirectory(seg.c_str());
      if(!nd) nd=dcur->mkdir(seg.c_str());
      dcur=nd;
    }
    return dcur?dcur:base;
  };

  for (int is=0;is<N_SPECS;++is) {
    const std::string sp(1,specs[is]);
    TDirectory* sdir=fout.mkdir(sp.c_str());

    // Var directories
    for(const auto&v:vars){
      TDirectory*vdir=sdir->mkdir(v.name.c_str());
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(vdir,cats[ic].dir); cdir->cd();
        for(int pl=0;pl<N_PLANES;++pl){
          TCanvas*c=new TCanvas((sp+"_c_"+v.name+"_pl"+plane_names[pl]).c_str(),
                                (sp+" "+v.name+" plane "+plane_names[pl]).c_str(),1600,1000);
          c->Divide(4,3);
          for(int pa=0;pa<N_PADDLES;++pa){c->cd(pa+1);hpad[is][ic][v.name][pl][pa]->DrawCopy();}
          wc(c);
        }
        TCanvas*cs=new TCanvas((sp+"_c_"+v.name+"_summary").c_str(),(sp+" "+v.name+" summary").c_str(),1600,1000);
        cs->Divide(3,2);
        for(int pl=0;pl<N_PLANES;++pl){cs->cd(pl+1);hsum[is][ic][v.name][pl]->DrawCopy();}
        cs->cd(6);htot[is][ic][v.name]->DrawCopy();wc(cs);
      }
      sdir->cd();
    }

    // Punchthrough
    {
      TDirectory*d=sdir->mkdir("punchthrough");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        TCanvas*c01=new TCanvas((sp+"_c_pt_000_001").c_str(),(sp+" PT 000-001 per paddle").c_str(),1600,1000);
        c01->Divide(4,3);for(int p=0;p<N_PADDLES;++p){c01->cd(p+1);h_pt_01_pad[is][ic][p]->DrawCopy();}wc(c01);
        TCanvas*c23=new TCanvas((sp+"_c_pt_100_101").c_str(),(sp+" PT 100-101 per paddle").c_str(),1600,1000);
        c23->Divide(4,3);for(int p=0;p<N_PADDLES;++p){c23->cd(p+1);h_pt_23_pad[is][ic][p]->DrawCopy();}wc(c23);
        TCanvas*cs=new TCanvas((sp+"_c_pt_sums").c_str(),(sp+" PT sums").c_str(),1600,600);cs->Divide(3,1);
        cs->cd(1);h_pt_01_sum[is][ic]->DrawCopy();cs->cd(2);h_pt_23_sum[is][ic]->DrawCopy();cs->cd(3);h_pt_tot[is][ic]->DrawCopy();wc(cs);
      }
      sdir->cd();
    }

    // edep vs punchthrough — two canvases per key (000-001 and 100-101)
    {
      TDirectory*d=sdir->mkdir("edep_vs_punchthrough");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        for(const auto&key:{"edep_0","edep_1","edep_amp_0","edep_amp_1"}){
          std::string k(key);
          TCanvas*c01=new TCanvas((sp+"_c_"+k+"_vs_pt_000_001").c_str(),(sp+" "+k+" vs PT 000-001").c_str(),1600,1000);
          c01->Divide(4,3);for(int p=0;p<N_PADDLES;++p){c01->cd(p+1);h_edep_vs_pt_01[is][ic][key][p]->DrawCopy("COLZ");}wc(c01);
          TCanvas*c23=new TCanvas((sp+"_c_"+k+"_vs_pt_100_101").c_str(),(sp+" "+k+" vs PT 100-101").c_str(),1600,1000);
          c23->Divide(4,3);for(int p=0;p<N_PADDLES;++p){c23->cd(p+1);h_edep_vs_pt_23[is][ic][key][p]->DrawCopy("COLZ");}wc(c23);
        }
      }
      sdir->cd();
    }

    // tof vs edep
    {
      TDirectory*d=sdir->mkdir("tof_vs_edep");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        for(const auto&key:{"edep_0","edep_1","edep_amp_0","edep_amp_1"}){
          TCanvas*c=new TCanvas((sp+"_c_tof_vs_"+std::string(key)).c_str(),(sp+" tof vs "+std::string(key)).c_str(),1600,1000);
          c->Divide(4,3);for(int p=0;p<N_PADDLES;++p){c->cd(p+1);h_tof_vs_edep[is][ic][key][p]->DrawCopy("COLZ");}wc(c);
        }
      }
      sdir->cd();
    }

    // Front veto
    {
      TDirectory*d=sdir->mkdir("front_veto");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        TCanvas*cht01=new TCanvas((sp+"_c_ht_fv_001").c_str(),(sp+" ht fv 001 per paddle").c_str(),1600,1000);
        cht01->Divide(4,3);for(int p=0;p<N_PADDLES;++p){cht01->cd(p+1);h_ht_fv_01[is][ic][p]->DrawCopy();}wc(cht01);
        TCanvas*cht23=new TCanvas((sp+"_c_ht_fv_101").c_str(),(sp+" ht fv 101 per paddle").c_str(),1600,1000);
        cht23->Divide(4,3);for(int p=0;p<N_PADDLES;++p){cht23->cd(p+1);h_ht_fv_23[is][ic][p]->DrawCopy();}wc(cht23);
        TCanvas*chts=new TCanvas((sp+"_c_ht_fv_sums").c_str(),(sp+" ht fv sums").c_str(),1600,600);chts->Divide(3,1);
        chts->cd(1);h_ht_fv_01_sum[is][ic]->DrawCopy();chts->cd(2);h_ht_fv_23_sum[is][ic]->DrawCopy();chts->cd(3);h_ht_fv_tot[is][ic]->DrawCopy();wc(chts);
        TCanvas*ctof01=new TCanvas((sp+"_c_tof_fv_001").c_str(),(sp+" tof fv 001 per paddle").c_str(),1600,1000);
        ctof01->Divide(4,3);for(int p=0;p<N_PADDLES;++p){ctof01->cd(p+1);h_tof_fv_01[is][ic][p]->DrawCopy();}wc(ctof01);
        TCanvas*ctof23=new TCanvas((sp+"_c_tof_fv_101").c_str(),(sp+" tof fv 101 per paddle").c_str(),1600,1000);
        ctof23->Divide(4,3);for(int p=0;p<N_PADDLES;++p){ctof23->cd(p+1);h_tof_fv_23[is][ic][p]->DrawCopy();}wc(ctof23);
        TCanvas*ctofs=new TCanvas((sp+"_c_tof_fv_sums").c_str(),(sp+" tof fv sums").c_str(),1600,600);ctofs->Divide(3,1);
        ctofs->cd(1);h_tof_fv_01_sum[is][ic]->DrawCopy();ctofs->cd(2);h_tof_fv_23_sum[is][ic]->DrawCopy();ctofs->cd(3);h_tof_fv_tot[is][ic]->DrawCopy();wc(ctofs);
      }
      sdir->cd();
    }

    // Corrected tof + background-subtracted
    {
      TDirectory*d=sdir->mkdir("corrected");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        for(int pl=0;pl<N_PLANES;++pl){
          TCanvas*ct=new TCanvas((sp+"_c_tof_corr_"+plane_names[pl]).c_str(),
                                  (sp+" tof corr "+plane_names[pl]).c_str(),1600,1000);
          ct->Divide(4,3);for(int p=0;p<N_PADDLES;++p){ct->cd(p+1);h_tof_corr[is][ic][pl][p]->DrawCopy();}wc(ct);
        }
        TCanvas*cts=new TCanvas((sp+"_c_tof_corr_summary").c_str(),(sp+" tof corr summary").c_str(),1600,1000);
        cts->Divide(3,2);
        for(int pl=0;pl<N_PLANES;++pl){cts->cd(pl+1);h_tof_corr_sum[is][ic][pl]->DrawCopy();}
        cts->cd(6);h_tof_corr_tot[is][ic]->DrawCopy();wc(cts);
      }
      sdir->cd();
    }

    // Background-subtracted corrected tof — separate top-level directory
    {
      TDirectory*d=sdir->mkdir("corrected_bgsub");
      for(int ic=0;ic<ncats;++ic){
        TDirectory*cdir=mkpath(d,cats[ic].dir);cdir->cd();
        for(int pl=0;pl<N_PLANES;++pl){
          TCanvas*cbg=new TCanvas((sp+"_c_tof_corr_bgsub_"+plane_names[pl]).c_str(),
                                   (sp+" tof corr bgsub "+plane_names[pl]).c_str(),1600,1000);
          cbg->Divide(4,3);
          for(int pa=0;pa<N_PADDLES;++pa){
            cbg->cd(pa+1);
            TH1D* hb=bgsub_tof(h_tof_corr[is][ic][pl][pa].GetPtr());
            hb->DrawCopy(); delete hb;
          } wc(cbg);
        }
        TCanvas*cbgs=new TCanvas((sp+"_c_tof_corr_bgsub_summary").c_str(),(sp+" tof corr bgsub summary").c_str(),1600,1000);
        cbgs->Divide(3,2);
        for(int pl=0;pl<N_PLANES;++pl){
          cbgs->cd(pl+1);
          TH1D* hb=bgsub_tof(h_tof_corr_sum[is][ic][pl].GetPtr());
          hb->DrawCopy(); delete hb;
        }
        cbgs->cd(6);
        {TH1D* hb=bgsub_tof(h_tof_corr_tot[is][ic].GetPtr());hb->DrawCopy();delete hb;}
        wc(cbgs);
      }
      sdir->cd();
    }

    // Proton-tagged tof: one sub-directory per tracking variant.  The proton-only
    // (no-track) histos are shared across variants; each variant folder draws them
    // alongside its own proton+track histos so the ratio is self-contained.
    {
      TDirectory*d=sdir->mkdir("proton_id");
      for(int it=0; it<ntracks; ++it){
        const std::string& tu=tracks[it].tsuf;
        TDirectory*td=mkpath(d,tracks[it].dir); td->cd();

        TCanvas*c=new TCanvas((sp+"_c_proton_tof").c_str(),(sp+" proton tof").c_str(),1800,600);
        c->Divide(3,1);
        c->cd(1); h_proton_tof[is]->DrawCopy();
        c->cd(2); h_proton_track_tof[is][it]->DrawCopy();
        c->cd(3);
        TH1D* hratio=(TH1D*)h_proton_track_tof[is][it]->Clone((sp+"_proton_track_ratio"+tu).c_str());
        hratio->SetTitle((sp+" proton track/total;tof-L/c(ns);ratio").c_str());
        hratio->Divide(h_proton_tof[is].GetPtr());
        hratio->DrawCopy(); delete hratio;
        wc(c);

        // Sideband-subtracted efficiency: (track-bg)/(total-bg)
        // Background = mean bin content in first 50 ns of x-range
        TCanvas*csb=new TCanvas((sp+"_c_proton_tof_sb").c_str(),(sp+" proton tof sb-sub").c_str(),1800,600);
        csb->Divide(3,1);
        TH1D* hp_rb=(TH1D*)h_proton_tof[is]->Clone((sp+"_proton_rb"+tu).c_str());       hp_rb->Rebin(PROTON_REBIN);
        TH1D* ht_rb=(TH1D*)h_proton_track_tof[is][it]->Clone((sp+"_proton_track_rb"+tu).c_str()); ht_rb->Rebin(PROTON_REBIN);
        TH1D* hpb=flat_bgsub(hp_rb);
        TH1D* htb=flat_bgsub(ht_rb);
        delete hp_rb; delete ht_rb;
        TH1D* hratio_sb=(TH1D*)htb->Clone((sp+"_proton_track_ratio_sb"+tu).c_str());
        hratio_sb->SetTitle((sp+" proton (track-bg)/(total-bg);tof-L/c(ns);ratio").c_str());
        hratio_sb->Divide(hpb);
        csb->cd(1); hpb->DrawCopy();
        csb->cd(2); htb->DrawCopy();
        csb->cd(3); hratio_sb->DrawCopy();
        wc(csb);
        delete hpb; delete htb; delete hratio_sb;

        // Per-paddle canvases for planes 001 and 101
        for(int pp=0;pp<2;++pp){
          const int pi=(pp==0)?1:3;
          const std::string pn=plane_names[pi];

          TCanvas*cp=new TCanvas((sp+"_c_proton_"+pn).c_str(),(sp+" proton "+pn+" per paddle").c_str(),1600,1000);
          cp->Divide(4,3);
          for(int pa=0;pa<N_PADDLES;++pa){cp->cd(pa+1);h_proton_pad[is][pp][pa]->DrawCopy();}
          wc(cp);

          TCanvas*ct=new TCanvas((sp+"_c_proton_track_"+pn).c_str(),(sp+" proton+track "+pn+" per paddle").c_str(),1600,1000);
          ct->Divide(4,3);
          for(int pa=0;pa<N_PADDLES;++pa){ct->cd(pa+1);h_proton_track_pad[is][it][pp][pa]->DrawCopy();}
          wc(ct);

          TCanvas*cr=new TCanvas((sp+"_c_proton_ratio_"+pn).c_str(),(sp+" proton ratio "+pn+" per paddle").c_str(),1600,1000);
          cr->Divide(4,3);
          for(int pa=0;pa<N_PADDLES;++pa){
            cr->cd(pa+1);
            TH1D* hr=(TH1D*)h_proton_track_pad[is][it][pp][pa]->Clone(
              (sp+"_proton_ratio_p"+std::to_string(pi)+"_b"+std::to_string(pa)+tu).c_str());
            hr->Divide(h_proton_pad[is][pp][pa].GetPtr());
            hr->DrawCopy(); delete hr;
          }
          wc(cr);
        }

        // Summary canvases: 3 pads (plane 001, plane 101, total)
        {
          TCanvas*csp2=new TCanvas((sp+"_c_proton_sum").c_str(),(sp+" proton plane sums").c_str(),1800,600);
          csp2->Divide(3,1);
          csp2->cd(1);h_proton_psum[is][0]->DrawCopy();
          csp2->cd(2);h_proton_psum[is][1]->DrawCopy();
          csp2->cd(3);h_proton_tof[is]->DrawCopy();
          wc(csp2);

          TCanvas*cst2=new TCanvas((sp+"_c_proton_track_sum").c_str(),(sp+" proton+track plane sums").c_str(),1800,600);
          cst2->Divide(3,1);
          cst2->cd(1);h_proton_track_psum[is][it][0]->DrawCopy();
          cst2->cd(2);h_proton_track_psum[is][it][1]->DrawCopy();
          cst2->cd(3);h_proton_track_tof[is][it]->DrawCopy();
          wc(cst2);

          TCanvas*csr2=new TCanvas((sp+"_c_proton_ratio_sum").c_str(),(sp+" proton ratio plane sums").c_str(),1800,600);
          csr2->Divide(3,1);
          for(int pp=0;pp<2;++pp){
            csr2->cd(pp+1);
            TH1D* hr=(TH1D*)h_proton_track_psum[is][it][pp]->Clone((sp+"_proton_ratio_sum_"+std::to_string(pp)+tu).c_str());
            hr->Divide(h_proton_psum[is][pp].GetPtr()); hr->DrawCopy(); delete hr;
          }
          csr2->cd(3);
          {TH1D* hr=(TH1D*)h_proton_track_tof[is][it]->Clone((sp+"_proton_ratio_total"+tu).c_str());
           hr->Divide(h_proton_tof[is].GetPtr()); hr->DrawCopy(); delete hr;}
          wc(csr2);
        }
      } // end tracking-variant loop
      sdir->cd();
    }

    fout.cd();
  } // end spec loop (output)

  fout.Close();
  std::cout<<"[lad_tof_fast] Done. Wrote "<<out_file<<"\n";
}
