#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>

using namespace std;
double mod_value = 4.00801; // RF period in ns

struct RunHistograms {
  int run_num = 0;
  bool valid = false;
  std::unique_ptr<TH1D> h_P;
  std::unique_ptr<TH1D> h_H;
  std::unique_ptr<TH1D> h_mod_P;
  std::unique_ptr<TH1D> h_mod_H;
};

RunHistograms processRun(int run_num) {
  RunHistograms result;
  result.run_num = run_num;

  auto makeEmpty = [&](const char *name, const char *title, int bins, double low, double high) {
    std::unique_ptr<TH1D> h(new TH1D(name, title, bins, low, high));
    h->SetDirectory(nullptr);
    return h;
  };

  auto makeEmptySet = [&]() {
    result.h_P = makeEmpty(Form("h_P_%d", run_num), "RF Offset P", 500, 0, 4000);
    result.h_H = makeEmpty(Form("h_H_%d", run_num), "RF Offset H", 500, 0, 4000);
    result.h_mod_P = makeEmpty(Form("h_mod_P_%d", run_num), "RF Offset Mod P", 500, 0, mod_value);
    result.h_mod_H = makeEmpty(Form("h_mod_H_%d", run_num), "RF Offset Mod H", 500, 0, mod_value);
    result.valid = true;
  };

  try {
    ROOT::RDataFrame df("T",
                        Form("/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                             "LAD_COIN_%d_0_0_50000.root",
                             run_num));

    auto dfP = df.Filter("g.evtyp==1").Define("P_rf", "P.ladkin.t_vertex-T.shms.pRF_tdcTime");
    auto dfH = df.Filter("g.evtyp==2").Define("H_rf", "H.ladkin.t_vertex-T.hms.hRF_tdcTime");

    // Helper lambda to get center 98th percentile range
    auto getPercentileRange = [&](TH1D *h, double &low, double &high) {
      double sum = h->Integral();
      if (sum <= 0) {
        low = 0;
        high = 4000;
        return;
      }
      double cumsum = 0;
      int low_bin = 0, high_bin = h->GetNbinsX();

      // Find the 1st percentile
      for (int i = 1; i <= h->GetNbinsX(); i++) {
        cumsum += h->GetBinContent(i);
        if (cumsum > sum * 0.01) {
          low_bin = i;
          break;
        }
      }

      // Find the 99th percentile
      cumsum = 0;
      for (int i = h->GetNbinsX(); i >= 1; i--) {
        cumsum += h->GetBinContent(i);
        if (cumsum > sum * 0.01) {
          high_bin = i;
          break;
        }
      }

      low  = h->GetBinCenter(low_bin);
      high = h->GetBinCenter(high_bin);
    };

    auto h_P_temp = dfP.Histo1D({Form("h_P_temp_%d", run_num), "RF Offset P Temp", 1000, 0, 4000}, "P_rf");
    double p_low, p_high;
    getPercentileRange(h_P_temp.GetPtr(), p_low, p_high);

    auto h_H_temp = dfH.Histo1D({Form("h_H_temp_%d", run_num), "RF Offset H Temp", 1000, 0, 4000}, "H_rf");
    double h_low, h_high;
    getPercentileRange(h_H_temp.GetPtr(), h_low, h_high);

    auto h_P = dfP.Histo1D({Form("h_P_%d", run_num), "RF Offset P", 500, p_low, p_high}, "P_rf");
    auto h_H = dfH.Histo1D({Form("h_H_%d", run_num), "RF Offset H", 500, h_low, h_high}, "H_rf");

    auto dfPmod = dfP.Define("P_rf_mod", Form("fmod(P_rf, %.5f)", mod_value));
    auto dfHmod = dfH.Define("H_rf_mod", Form("fmod(H_rf, %.5f)", mod_value));

    auto h_mod_P = dfPmod.Histo1D({Form("h_mod_P_%d", run_num), "RF Offset Mod P", 500, 0, mod_value}, "P_rf_mod");
    auto h_mod_H = dfHmod.Histo1D({Form("h_mod_H_%d", run_num), "RF Offset Mod H", 500, 0, mod_value}, "H_rf_mod");

    result.h_P.reset((TH1D *)h_P->Clone());
    result.h_H.reset((TH1D *)h_H->Clone());
    result.h_mod_P.reset((TH1D *)h_mod_P->Clone());
    result.h_mod_H.reset((TH1D *)h_mod_H->Clone());

    result.h_P->SetDirectory(nullptr);
    result.h_H->SetDirectory(nullptr);
    result.h_mod_P->SetDirectory(nullptr);
    result.h_mod_H->SetDirectory(nullptr);

    result.valid = true;
    return result;
  } catch (...) {
    makeEmptySet();
    return result;
  }
}

void get_RF_offset_plot(const char* run_list_file = "run_numbers_all.txt") {

  gROOT->SetBatch(kTRUE);
  ROOT::EnableThreadSafety();
  TH1::AddDirectory(kFALSE);

  // Read run numbers from file
  vector<int> run_nums;
  ifstream infile(run_list_file);
  if (!infile.is_open()) {
    cerr << "Error: Could not open run list file: " << run_list_file << endl;
    cerr << "Please create a file with one run number per line." << endl;
    return;
  }
  
  int run_num;
  while (infile >> run_num) {
    run_nums.push_back(run_num);
  }
  infile.close();
  
  if (run_nums.empty()) {
    cerr << "Error: No run numbers found in " << run_list_file << endl;
    return;
  }
  
  cout << "Read " << run_nums.size() << " run numbers from " << run_list_file << endl;
  std::sort(run_nums.begin(), run_nums.end());

  const unsigned int max_threads = 20;
  unsigned int thread_count = std::min<unsigned int>(max_threads,
                                                     std::max<unsigned int>(1, std::thread::hardware_concurrency()));
  ROOT::EnableImplicitMT(thread_count);

  int total_runs = run_nums.size();
  std::vector<RunHistograms> results;
  results.reserve(total_runs);

  int current = 0;
  for (const auto &run_num : run_nums) {
    RunHistograms res = processRun(run_num);
    results.push_back(std::move(res));
    current++;
    double percent = (100.0 * current) / total_runs;
    cout << "\rProgress: " << current << "/" << total_runs << " (" << fixed << setprecision(1) << percent << "%) " << flush;
  }
  cout << endl; // Final newline after progress completes

  std::sort(results.begin(), results.end(), [](const RunHistograms &a, const RunHistograms &b) {
    return a.run_num < b.run_num;
  });

  TFile *out = new TFile("rf_offset_histograms_all.root", "RECREATE");

  TDirectory *dir_h_P = out->mkdir("h_P");
  TDirectory *dir_h_mod_P = out->mkdir("h_mod_P");
  TDirectory *dir_h_H = out->mkdir("h_H");
  TDirectory *dir_h_mod_H = out->mkdir("h_mod_H");

  for (const auto &res : results) {
    if (!res.valid) {
      continue;
    }

    if (res.h_P) {
      dir_h_P->cd();
      res.h_P->Write();
      out->cd();
    }
    if (res.h_mod_P) {
      dir_h_mod_P->cd();
      res.h_mod_P->Write();
      out->cd();
    }
    if (res.h_H) {
      dir_h_H->cd();
      res.h_H->Write();
      out->cd();
    }
    if (res.h_mod_H) {
      dir_h_mod_H->cd();
      res.h_mod_H->Write();
      out->cd();
    }
  }

  out->Close();
}
