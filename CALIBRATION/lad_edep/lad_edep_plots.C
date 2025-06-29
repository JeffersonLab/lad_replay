// Lucas Ehinger
// General LAD hodo plotting script
// Multithreaded, to speed up the processing
// Not great when it comes to memory management (requires lots of memory). No errors or major memory leaks, just
// duplication of memory across threads since root doesn't like mutexes. Multi-threading with ROOT (which is inherently
// not thread-safe) is hard, and I wasn't smart enough to do it elegantly, but this works.
#include </usr/lib/gcc/x86_64-redhat-linux/11/include/omp.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

const int MAX_DATA     = 100;
const int MAX_DATA_GEM = 1000;
const int CACHE_SIZE   = 50 * 1024 * 1024; // 50 MB cache size
const int MAX_EVTS_PER_THREAD = 500000;

struct hist_params {
  int NBINS;
  double MIN;
  double MAX;
};

const hist_params time_params    = {100, 1700.0, 2000.0};
const hist_params dt_params      = {200, -2, 8};
const hist_params adc_amp_params = {200, 0.0, 400.0};
const hist_params adc_int_params = {200, 0.0, 40.0};

const double TDC2NS = 0.09766; // TDC to ns conversion factor
const double ADC2NS = 0.0625;  // ADC to ns conversion factor

const int MIN_TDC_TIME = 1600;
const int MAX_TDC_TIME = 2000;

const char spec_prefix = 'H'; // Spectrometer to replay

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

double proton_cut_time[2][N_PADDLES] = {{3.5, 3, 5.5, 3, 3, 5, 4, 3.5, 3.5, 2.5, 3},
                                        {3, 0, 4, 5, 5, 2, 4.5, 4.5, 4, 0, 4}};

const double janky_diff_time_calib[2][N_PADDLES] = {{-0.25, 0.4, -1.6, 0.5, 0.7, -1.5, 0, 0.6, 0.3, 1, 0.7},
                                                    {1.1, 0.2, -0.1, -0.4, -0.9, 2, -0.7, -0.2, 0.25, 2, -0.3}};

TLegend *make_legend(TCanvas *c) {
  // Create and draw a new legend
  TLegend *legend = new TLegend(0.7, 0.6, 0.9, 0.9);
  int nhists      = c->GetListOfPrimitives()->GetSize();
  for (int i = 0; i < nhists; ++i) {
    TObject *obj = c->GetListOfPrimitives()->At(i);
    if (obj && obj->InheritsFrom(TH1::Class())) {
      // Only add if not already in the legend
      bool alreadyInLegend = false;
      for (int j = 0; j < legend->GetListOfPrimitives()->GetSize(); ++j) {
        TObject *legObj = legend->GetListOfPrimitives()->At(j);
        if (legObj && legObj->InheritsFrom("TLegendEntry")) {
          TLegendEntry *entry = (TLegendEntry *)legObj;
          if (entry->GetObject() == obj) {
            alreadyInLegend = true;
            break;
          }
        }
      }
      if (!alreadyInLegend) {
        legend->AddEntry(obj, obj->GetTitle(), "l");
      }
    }
  }
  return legend;
}

template <typename HistType>
void write_to_canvas_plane(HistType *hist_arr[N_PLANES][N_PADDLES], TFile *file, TString dir, TString var_name,
                           bool include_comp_plt = false) {
  // Create and navigate to the directory
  file->mkdir(dir);
  file->cd(dir);
  // Create a canvas for each plane
  TH1F *hist_plane_sum[N_PLANES] = {nullptr};
  for (int plane = 0; plane < N_PLANES; ++plane) {
    // Create canvas for top bars
    TCanvas *c1 = new TCanvas(Form("c_%s_plane_%s", var_name.Data(), plane_names[plane].c_str()),
                              Form("%s Plane %s", var_name.Data(), plane_names[plane].c_str()), 800, 600);
    c1->Divide(4, 3); // Divide canvas into subpads
    for (int bar = 0; bar < N_PADDLES; ++bar) {
      c1->cd(bar + 1);
      hist_arr[plane][bar]->Draw();
    }
    c1->Write();
    delete c1;

    // Sum all bars for this plane into a single histogram
    if (!hist_plane_sum[plane]) {
      hist_plane_sum[plane] = (TH1F *)hist_arr[plane][0]->Clone(
          Form("h_%s_plane_%s_all_bars", var_name.Data(), plane_names[plane].c_str()));
      hist_plane_sum[plane]->Reset();
      hist_plane_sum[plane]->SetTitle(Form("%s Plane %s All Bars", var_name.Data(), plane_names[plane].c_str()));
    }
    for (int bar = 0; bar < N_PADDLES; ++bar) {
      hist_plane_sum[plane]->Add(hist_arr[plane][bar]);
    }

    // Write the summed histogram for this plane
    // hist_plane_sum[plane]->Write();
  }

  // Create a canvas to overlay all planes
  TCanvas *c_all_planes =
      new TCanvas(Form("c_%s_all_planes", var_name.Data()), Form("%s All Planes", var_name.Data()), 1200, 800);
  c_all_planes->Divide(2, 3); // Divide canvas into subpads
  for (int plane = 0; plane < N_PLANES; ++plane) {
    c_all_planes->cd(plane + 1);
    hist_plane_sum[plane]->SetLineColor(plane + 1); // Assign different colors for each plane
    hist_plane_sum[plane]->Draw("HIST");
  }
  c_all_planes->Write();
  delete c_all_planes;

  // Sum all plane histograms into one final histogram
  TH1F *hist_final_sum = (TH1F *)hist_plane_sum[0]->Clone(Form("h_%s_final_sum", var_name.Data()));
  hist_final_sum->Reset();
  hist_final_sum->SetTitle(Form("%s All Planes", var_name.Data()));
  for (int plane = 0; plane < N_PLANES; ++plane) {
    hist_final_sum->Add(hist_plane_sum[plane]);
  }

  // Write the final summed histogram
  hist_final_sum->Write();

  delete hist_final_sum;
  for (int plane = 0; plane < N_PLANES; ++plane) {
    delete hist_plane_sum[plane];
  }
  return;
}

template <typename T> void add_branch(TChain *tree, const char *branch_name, T *branch_data) {
  // Add a branch to the tree
  tree->SetBranchAddress(branch_name, branch_data);
  tree->AddBranchToCache(branch_name, kTRUE);
  tree->SetBranchStatus(branch_name, 1); // Enable the branch
}

void process_chunk(int i_thread, int start, int end, std::vector<TString> &fileNames,
                   map<string, TH1 *[N_PLANES][N_PADDLES]> &hist_map) {

  TChain *T = new TChain("T");
  for (const auto &fileName : fileNames) {
    T->Add(fileName);
  }

  if (T->GetEntries() == 0) {
    std::cerr << "Error: Cannot open the ROOT files or no entries in the TChain!" << std::endl;
    return;
  }

  // Define arrays to hold the data
  // Double_t trk_d0[MAX_DATA_GEM];
  // Double_t trk_d0_good[MAX_DATA_GEM];
  // Double_t trk_projz[MAX_DATA_GEM];
  // Double_t trk_xloc[2][MAX_DATA_GEM], trk_yloc[2][MAX_DATA_GEM];
  // Double_t kin_trackID_0[MAX_DATA_GEM], kin_trackID_1[MAX_DATA_GEM];
  // Double_t kin_plane_0[MAX_DATA_GEM], kin_plane_1[MAX_DATA_GEM];
  // Double_t kin_paddle_0[MAX_DATA_GEM], kin_paddle_1[MAX_DATA_GEM];
  // Double_t kin_hittime_0[MAX_DATA_GEM], kin_hittime_1[MAX_DATA_GEM];
  // Double_t kin_hit_tof_0[MAX_DATA_GEM], kin_hit_tof_1[MAX_DATA_GEM];
  // Double_t kin_hittheta_0[MAX_DATA_GEM], kin_hittheta_1[MAX_DATA_GEM];
  // Double_t kin_hitphi_0[MAX_DATA_GEM], kin_hitphi_1[MAX_DATA_GEM];
  // Double_t kin_hitedep_0[MAX_DATA_GEM], kin_hitedep_1[MAX_DATA_GEM];
  // Double_t kin_dTrkHoriz_0[MAX_DATA_GEM], kin_dTrkHoriz_1[MAX_DATA_GEM];
  // Double_t kin_dTrkVert_0[MAX_DATA_GEM], kin_dTrkVert_1[MAX_DATA_GEM];
  // Int_t nTracks, nGoodHits;

  // Hodo-level data

  Double_t fullhit_time_avg[N_PLANES][MAX_DATA], fullhit_adc_avg[N_PLANES][MAX_DATA],
      fullhit_paddle[N_PLANES][MAX_DATA], fullhit_adc_amp_avg[N_PLANES][MAX_DATA], fullhit_ypos[N_PLANES][MAX_DATA];
  Int_t fullhit_n[N_PLANES];

  Double_t hodo_start_time, pTRIG1, pTRIG2, pTRIG3, pTRIG4, evtyp;

  // add_branch(T, Form("Ndata.%c.gem.trk.d0", spec_prefix), &nTracks);
  // add_branch(T, Form("Ndata.%c.ladkin.goodhit_trackid_0", spec_prefix), &nGoodHits);
  // add_branch(T, Form("%c.gem.trk.d0", spec_prefix), &trk_d0);
  // add_branch(T, Form("%c.gem.trk.d0_good", spec_prefix), &trk_d0_good);
  // add_branch(T, Form("%c.gem.trk.projz", spec_prefix), &trk_projz);
  // add_branch(T, Form("%c.gem.trk.x1_local", spec_prefix), &trk_xloc[0]);
  // add_branch(T, Form("%c.gem.trk.x2_local", spec_prefix), &trk_xloc[1]);
  // add_branch(T, Form("%c.gem.trk.y1_local", spec_prefix), &trk_yloc[0]);
  // add_branch(T, Form("%c.gem.trk.y2_local", spec_prefix), &trk_yloc[1]);
  // add_branch(T, Form("%c.ladkin.goodhit_trackid_0", spec_prefix), &kin_trackID_0);
  // add_branch(T, Form("%c.ladkin.goodhit_trackid_1", spec_prefix), &kin_trackID_1);
  // add_branch(T, Form("%c.ladkin.goodhit_plane_0", spec_prefix), &kin_plane_0);
  // add_branch(T, Form("%c.ladkin.goodhit_plane_1", spec_prefix), &kin_plane_1);
  // add_branch(T, Form("%c.ladkin.goodhit_paddle_0", spec_prefix), &kin_paddle_0);
  // add_branch(T, Form("%c.ladkin.goodhit_paddle_1", spec_prefix), &kin_paddle_1);
  // add_branch(T, Form("%c.ladkin.goodhit_hittime_0", spec_prefix), &kin_hittime_0);
  // add_branch(T, Form("%c.ladkin.goodhit_hittime_1", spec_prefix), &kin_hittime_1);
  // add_branch(T, Form("%c.ladkin.goodhit_tof_0", spec_prefix), &kin_hit_tof_0);
  // add_branch(T, Form("%c.ladkin.goodhit_tof_1", spec_prefix), &kin_hit_tof_1);
  // add_branch(T, Form("%c.ladkin.goodhit_hittheta_0", spec_prefix), &kin_hittheta_0);
  // add_branch(T, Form("%c.ladkin.goodhit_hittheta_1", spec_prefix), &kin_hittheta_1);
  // add_branch(T, Form("%c.ladkin.goodhit_hitphi_0", spec_prefix), &kin_hitphi_0);
  // add_branch(T, Form("%c.ladkin.goodhit_hitphi_1", spec_prefix), &kin_hitphi_1);
  // add_branch(T, Form("%c.ladkin.goodhit_hitedep_0", spec_prefix), &kin_hitedep_0);
  // add_branch(T, Form("%c.ladkin.goodhit_hitedep_1", spec_prefix), &kin_hitedep_1);
  // add_branch(T, Form("%c.ladkin.goodhit_dTrkHoriz_0", spec_prefix), &kin_dTrkHoriz_0);
  // add_branch(T, Form("%c.ladkin.goodhit_dTrkHoriz_1", spec_prefix), &kin_dTrkHoriz_1);
  // add_branch(T, Form("%c.ladkin.goodhit_dTrkVert_0", spec_prefix), &kin_dTrkVert_0);
  // add_branch(T, Form("%c.ladkin.goodhit_dTrkVert_1", spec_prefix), &kin_dTrkVert_1);

  T->SetCacheSize(CACHE_SIZE);
  T->SetBranchStatus("*", 0); // Disable all branches initially
  add_branch(T, Form("%c.hod.starttime", spec_prefix), &hodo_start_time);
  add_branch(T, "T.hms.pTRIG1_tdcTime", &pTRIG1);
  add_branch(T, "T.hms.pTRIG2_tdcTime", &pTRIG2);
  add_branch(T, "T.hms.pTRIG3_tdcTime", &pTRIG3);
  add_branch(T, "T.hms.pTRIG4_tdcTime", &pTRIG4);
  add_branch(T, "g.evtyp", &evtyp);

  for (int plane = 0; plane < N_PLANES; ++plane) {
    add_branch(T, Form("%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()), &fullhit_time_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitEdep", spec_prefix, plane_names[plane].c_str()), &fullhit_adc_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitEdepAmp", spec_prefix, plane_names[plane].c_str()),
               &fullhit_adc_amp_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitPaddleNum", spec_prefix, plane_names[plane].c_str()),
               &fullhit_paddle[plane]);
    add_branch(T, Form("Ndata.%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()), &fullhit_n[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitPos", spec_prefix, plane_names[plane].c_str()), &fullhit_ypos[plane]);
  }

  Double_t hodo_hit_time[N_PLANES][N_PADDLES];
  Double_t hodo_hit_edep[N_PLANES][N_PADDLES];
  Double_t hodo_hit_adc_amp[N_PLANES][N_PADDLES];
  Double_t hodo_hit_ypos[N_PLANES][N_PADDLES];
  ////////////////////////////////////////////////////
  // Start Event Loop

  for (int i = start; i < end; ++i) {
    T->GetEntry(i);
    if ((spec_prefix == 'H' && int(evtyp) != 2) || (spec_prefix == 'P' && int(evtyp) != 1))
      continue; // Only process events with the correct evtyp based on spec_prefix

    // clear has_hodo_hit array

    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        hodo_hit_time[plane][bar]    = -999;
        hodo_hit_edep[plane][bar]    = -999;
        hodo_hit_adc_amp[plane][bar] = -999;
        hodo_hit_ypos[plane][bar]    = -999;
      }
    }

    // Loop over the hits and fill the arrays
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int i_top = 0; i_top < fullhit_n[plane]; ++i_top) {
        int bar = int(fullhit_paddle[plane][i_top]) - 1;
        hist_map["h_time_avg"][plane][bar]->Fill(fullhit_time_avg[plane][i_top]);
        hist_map["h_ADC_Int"][plane][bar]->Fill(fullhit_adc_avg[plane][i_top]);
        hist_map["h_ADC_Amp"][plane][bar]->Fill(fullhit_adc_amp_avg[plane][i_top]);
        hist_map["h_TDC_vs_ADC_Int"][plane][bar]->Fill(fullhit_time_avg[plane][i_top], fullhit_adc_avg[plane][i_top]);
        hist_map["h_TDC_vs_ADC_Amp"][plane][bar]->Fill(fullhit_time_avg[plane][i_top],
                                                       fullhit_adc_amp_avg[plane][i_top]);

        if (bar < 0 || bar >= N_PADDLES)
          continue; // Skip invalid bars
        hodo_hit_time[plane][bar]    = fullhit_time_avg[plane][i_top];
        hodo_hit_edep[plane][bar]    = fullhit_adc_avg[plane][i_top];
        hodo_hit_adc_amp[plane][bar] = fullhit_adc_amp_avg[plane][i_top];
        hodo_hit_ypos[plane][bar]    = fullhit_ypos[plane][i_top];
      }
    }

    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        int matching_plane;
        switch (plane) {
        case 0:
          matching_plane = 1; // Plane 001
          break;
        case 1:
          matching_plane = 0; // Plane 000
          break;
        case 2:
          matching_plane = 3; // Plane 001
          break;
        case 3:
          matching_plane = 2; // Plane 100
          break;
        case 4:
          matching_plane = 4; // Plane 200
          break;
        default:
          continue; // Skip invalid planes
        }
        if (hodo_hit_time[matching_plane][bar] != -999 && hodo_hit_time[plane][bar] != -999) {
          double tdc_diff = hodo_hit_time[matching_plane][bar] - hodo_hit_time[plane][bar];
          if (matching_plane < plane) {
            tdc_diff = -tdc_diff; // Reverse sign for the second plane
          }
          // if (plane / 2 < 2)
          //   tdc_diff += janky_diff_time_calib[plane / 2][bar];
          // Fill the histograms with the TDC difference
          hist_map["h_time_avg_punchthrough"][plane][bar]->Fill(hodo_hit_time[plane][bar]);
          hist_map["h_TDC_Diff_bar"][plane][bar]->Fill(tdc_diff);
          hist_map["h_TDC_Diff_vs_ADC_Int"][plane][bar]->Fill(tdc_diff, hodo_hit_edep[plane][bar]);
          hist_map["h_TDC_Diff_vs_ADC_Amp"][plane][bar]->Fill(tdc_diff, hodo_hit_adc_amp[plane][bar]);
          // Fill the Front vs Back ADC Integral and Amplitude histograms
          hist_map["h_Front_Back_ADC_Int"][plane][bar]->Fill(hodo_hit_edep[plane][bar],
                                                             hodo_hit_edep[matching_plane][bar]);
          hist_map["h_Front_Back_ADC_Amp"][plane][bar]->Fill(hodo_hit_adc_amp[plane][bar],
                                                             hodo_hit_adc_amp[matching_plane][bar]);
          hist_map["h_HitYPos"][plane][bar]->Fill(hodo_hit_ypos[plane][bar]);

          if (tdc_diff > proton_cut_time[plane / 2][bar]) {
            hist_map["h_time_avg_protons"][plane][bar]->Fill(hodo_hit_time[plane][bar]);
            // hist_map["h_TDC_Diff_vs_ADC_Int_protons"][plane][bar]->Fill(tdc_diff, hodo_hit_edep[plane][bar]);
            // hist_map["h_TDC_Diff_vs_ADC_Amp_protons"][plane][bar]->Fill(tdc_diff, hodo_hit_adc_amp[plane][bar]);
          }
        }
      }
    }

    // Loop over the hits and fill the histograms
    // Print the status as a percentage
    if (i % ((end - start) / 100) == 0 && i_thread == 0) {
      std::cout << "\rProcessing: " << int((i - start) * 100.0 / (end - start)) << "% completed." << std::flush;
    }
  } // End Event Loop

  return;
}

int lad_edep_plots() {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);
  ROOT::EnableThreadSafety();

  // Open multiple ROOT files
  std::vector<TString> fileNames = {
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22591_0_2_-1.root"};
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22609_0_6_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22610_0_6_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22611_0_6_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22613_0_6_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22615_0_6_-1.root",
      "/cache/hallc/c-lad/analysis/ehingerl/online_v1/LAD_COIN_23105_0_6_-1.root",
      "/cache/hallc/c-lad/analysis/ehingerl/online_v1/LAD_COIN_23106_0_6_-1.root",
      "/cache/hallc/c-lad/analysis/ehingerl/online_v1/LAD_COIN_23107_0_6_-1.root",
      "/cache/hallc/c-lad/analysis/ehingerl/online_v1/LAD_COIN_23108_0_6_-1.root",
  };
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22616_0_21_-1.root"};
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22382_0_21_-1_1.root"};
  TString outputFileName = Form("lad_edep_plots_new_timing_%c.root", spec_prefix);

  // Create a TChain to combine the trees from multiple files
  TChain *chain = new TChain("T");
  for (const auto &fileName : fileNames) {
    chain->Add(fileName);
  }

  if (chain->GetEntries() == 0) {
    std::cerr << "Error: Cannot open the ROOT files or no entries in the TChain!" << std::endl;
    delete chain;
    return -1;
  }

  // Get the TTree
  TTree *T = chain;
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    delete chain;
    return -1;
  }

  // Number of entries in the TTree
  int nEntries = T->GetEntries();

  // Number of threads to use
  int numThreads = std::thread::hardware_concurrency();
  int chunkSize  = nEntries / numThreads;

  // Adjust the number of threads if the chunk size is too small
  if (chunkSize < MAX_EVTS_PER_THREAD) {
    numThreads = std::max(1, nEntries / MAX_EVTS_PER_THREAD);
    chunkSize  = nEntries / numThreads;
  }

  std::vector<map<string, TH1 *[N_PLANES][N_PADDLES]>> hist_map_vec(numThreads);

  for (int thread = 0; thread < numThreads; ++thread) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        hist_map_vec[thread]["h_time_avg"][plane][bar] =
            new TH1D(Form("h_TDC_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar), time_params.NBINS,
                     time_params.MIN, time_params.MAX);
        hist_map_vec[thread]["h_time_avg"][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
        hist_map_vec[thread]["h_time_avg"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_time_avg_punchthrough"][plane][bar] =
            new TH1D(Form("h_TDC_punchthrough_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Punchthrough Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     time_params.NBINS, time_params.MIN, time_params.MAX);
        hist_map_vec[thread]["h_time_avg_punchthrough"][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
        hist_map_vec[thread]["h_time_avg_punchthrough"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_time_avg_protons"][plane][bar] =
            new TH1D(Form("h_TDC_protons_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Protons Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     time_params.NBINS, time_params.MIN, time_params.MAX);
        hist_map_vec[thread]["h_time_avg_protons"][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
        hist_map_vec[thread]["h_time_avg_protons"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_TDC_Diff_bar"][plane][bar] =
            new TH1D(Form("h_TDC_Diff_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Diff Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     dt_params.NBINS, dt_params.MIN, dt_params.MAX);
        hist_map_vec[thread]["h_TDC_Diff_bar"][plane][bar]->GetXaxis()->SetTitle("Time Difference (ns)");
        hist_map_vec[thread]["h_TDC_Diff_bar"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_ADC_Int"][plane][bar] =
            new TH1D(Form("h_ADC_Int_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("ADC Int Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     adc_int_params.NBINS, adc_int_params.MIN, adc_int_params.MAX);
        hist_map_vec[thread]["h_ADC_Int"][plane][bar]->GetXaxis()->SetTitle("ADC Integral (pC)");
        hist_map_vec[thread]["h_ADC_Int"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_ADC_Amp"][plane][bar] =
            new TH1D(Form("h_ADC_amp_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("ADC Amp Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     adc_amp_params.NBINS, adc_amp_params.MIN, adc_amp_params.MAX);
        hist_map_vec[thread]["h_ADC_Amp"][plane][bar]->GetXaxis()->SetTitle("ADC Amplitude (mV)");
        hist_map_vec[thread]["h_ADC_Amp"][plane][bar]->GetYaxis()->SetTitle("Counts");

        hist_map_vec[thread]["h_TDC_vs_ADC_Int"][plane][bar] =
            new TH2D(Form("h_TDC_vs_ADC_Int_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC vs ADC Int Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     time_params.NBINS, time_params.MIN, time_params.MAX, adc_int_params.NBINS, adc_int_params.MIN,
                     adc_int_params.MAX);
        hist_map_vec[thread]["h_TDC_vs_ADC_Int"][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
        hist_map_vec[thread]["h_TDC_vs_ADC_Int"][plane][bar]->GetYaxis()->SetTitle("ADC Integral (pC)");

        hist_map_vec[thread]["h_TDC_vs_ADC_Amp"][plane][bar] =
            new TH2D(Form("h_TDC_vs_ADC_Amp_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC vs ADC Amp Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     time_params.NBINS, time_params.MIN, time_params.MAX, adc_amp_params.NBINS, adc_amp_params.MIN,
                     adc_amp_params.MAX);
        hist_map_vec[thread]["h_TDC_vs_ADC_Amp"][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
        hist_map_vec[thread]["h_TDC_vs_ADC_Amp"][plane][bar]->GetYaxis()->SetTitle("ADC Amplitude (mV)");

        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Amp"][plane][bar] =
            new TH2D(Form("h_TDC_Diff_vs_ADC_amp_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Diff vs ADC Amp Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     dt_params.NBINS, dt_params.MIN, dt_params.MAX, adc_amp_params.NBINS, adc_amp_params.MIN,
                     adc_amp_params.MAX);
        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Amp"][plane][bar]->GetXaxis()->SetTitle("Time Difference (ns)");
        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Amp"][plane][bar]->GetYaxis()->SetTitle("ADC Amplitude (mV)");

        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Int"][plane][bar] =
            new TH2D(Form("h_TDC_Diff_vs_ADC_Int_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("TDC Diff vs ADC Int Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     dt_params.NBINS, dt_params.MIN, dt_params.MAX, adc_int_params.NBINS, adc_int_params.MIN,
                     adc_int_params.MAX);
        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Int"][plane][bar]->GetXaxis()->SetTitle("Time Difference (ns)");
        hist_map_vec[thread]["h_TDC_Diff_vs_ADC_Int"][plane][bar]->GetYaxis()->SetTitle("ADC Integral (pC)");

        // Histogram for Front vs Back ADC Integral
        hist_map_vec[thread]["h_Front_Back_ADC_Int"][plane][bar] =
            new TH2D(Form("h_Front_Back_ADC_Int_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("Front vs Back ADC Int Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     adc_int_params.NBINS, adc_int_params.MIN, adc_int_params.MAX, adc_int_params.NBINS,
                     adc_int_params.MIN, adc_int_params.MAX);
        hist_map_vec[thread]["h_Front_Back_ADC_Int"][plane][bar]->GetXaxis()->SetTitle("Front ADC Integral (pC)");
        hist_map_vec[thread]["h_Front_Back_ADC_Int"][plane][bar]->GetYaxis()->SetTitle("Back ADC Integral (pC)");

        // Histogram for Front vs Back ADC Amplitude
        hist_map_vec[thread]["h_Front_Back_ADC_Amp"][plane][bar] =
            new TH2D(Form("h_Front_Back_ADC_Amp_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
                     Form("Front vs Back ADC Amp Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar),
                     adc_amp_params.NBINS, adc_amp_params.MIN, adc_amp_params.MAX, adc_amp_params.NBINS,
                     adc_amp_params.MIN, adc_amp_params.MAX);
        hist_map_vec[thread]["h_Front_Back_ADC_Amp"][plane][bar]->GetXaxis()->SetTitle("Front ADC Amplitude (mV)");
        hist_map_vec[thread]["h_Front_Back_ADC_Amp"][plane][bar]->GetYaxis()->SetTitle("Back ADC Amplitude (mV)");

        hist_map_vec[thread]["h_HitYPos"][plane][bar] = new TH1D(
            Form("h_HitYPos_thread_%d_plane_%s_bar_%d", thread, plane_names[plane].c_str(), bar),
            Form("Hit Y Position Thread %d Plane %s Bar %d", thread, plane_names[plane].c_str(), bar), 100, -250, 250);
        hist_map_vec[thread]["h_HitYPos"][plane][bar]->GetXaxis()->SetTitle("Y Position (cm)");
        hist_map_vec[thread]["h_HitYPos"][plane][bar]->GetYaxis()->SetTitle("Counts");
      }
    }
  }

  // Create an output ROOT file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return -1;
  }

  // start threads
  cout << "Starting " << numThreads << " threads..." << endl;
  std::vector<std::thread> threads;
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    int start = i_thread * chunkSize;
    int end   = (i_thread == numThreads - 1) ? nEntries : start + chunkSize;
    threads.emplace_back(process_chunk, i_thread, start, end, ref(fileNames), ref(hist_map_vec[i_thread]));
  }
  // Wait for all threads to finish
  for (auto &thread : threads) {
    thread.join();
  }
  std::cout << "\rProcessing: 100% completed. \nMerging histograms" << std::endl;
  // Merge histograms from all threads by adding to the first thread
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    for (const auto &hist_pair : hist_map_vec[i_thread]) {

      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int bar = 0; bar < N_PADDLES; ++bar) {
          if (i_thread == 0) {
            (hist_map_vec[0][hist_pair.first][plane][bar])
                ->SetTitle(TString(hist_map_vec[0][hist_pair.first][plane][bar]->GetTitle())
                               .ReplaceAll(Form(" Thread %d", i_thread), ""));
            (hist_map_vec[0][hist_pair.first][plane][bar])
                ->SetName(TString(hist_map_vec[0][hist_pair.first][plane][bar]->GetName())
                              .ReplaceAll(Form("_thread_%d", i_thread), ""));
          } else {
            (hist_map_vec[0][hist_pair.first][plane][bar])->Add(hist_pair.second[plane][bar]);
          }
        }
      }
    }
  }

  cout << "Writing histograms to file..." << endl;
  // Write histograms to the output file
  outputFile->cd();

  write_to_canvas_plane(hist_map_vec[0]["h_time_avg"], outputFile, "KIN/TDC", "TDC");
  write_to_canvas_plane(hist_map_vec[0]["h_time_avg_punchthrough"], outputFile, "KIN/TDC_PUNCHTHROUGH",
                        "TDC_PUNCHTHROUGH");
  write_to_canvas_plane(hist_map_vec[0]["h_TDC_Diff_bar"], outputFile, "KIN/TDC_DIFF", "TDC_DIFF");
  write_to_canvas_plane(hist_map_vec[0]["h_ADC_Int"], outputFile, "KIN/ADC_INT", "ADC_INT");
  write_to_canvas_plane(hist_map_vec[0]["h_ADC_Amp"], outputFile, "KIN/ADC_AMP", "ADC_AMP");
  write_to_canvas_plane(hist_map_vec[0]["h_TDC_vs_ADC_Int"], outputFile, "KIN/TDC_VS_ADC_INT", "TDC_VS_ADC_INT");
  write_to_canvas_plane(hist_map_vec[0]["h_TDC_vs_ADC_Amp"], outputFile, "KIN/TDC_VS_ADC_AMP", "TDC_VS_ADC_AMP");
  write_to_canvas_plane(hist_map_vec[0]["h_TDC_Diff_vs_ADC_Amp"], outputFile, "KIN/TDC_DIFF_VS_ADC_AMP",
                        "TDC_DIFF_VS_ADC_AMP");
  write_to_canvas_plane(hist_map_vec[0]["h_TDC_Diff_vs_ADC_Int"], outputFile, "KIN/TDC_DIFF_VS_ADC_INT",
                        "TDC_DIFF_VS_ADC_INT");
  write_to_canvas_plane(hist_map_vec[0]["h_Front_Back_ADC_Int"], outputFile, "KIN/FRONT_BACK_ADC_INT",
                        "FRONT_BACK_ADC_INT");
  write_to_canvas_plane(hist_map_vec[0]["h_Front_Back_ADC_Amp"], outputFile, "KIN/FRONT_BACK_ADC_AMP",
                        "FRONT_BACK_ADC_AMP");
  write_to_canvas_plane(hist_map_vec[0]["h_time_avg_protons"], outputFile, "KIN/TDC_PROTONS", "TDC_PROTONS");

  write_to_canvas_plane(hist_map_vec[0]["h_HitYPos"], outputFile, "KIN/HIT_Y_POS", "HIT_Y_POS");

  outputFile->Close();
  delete outputFile;

  std::cout << "Done! All histograms written to file. Cleaning up..." << std::endl;

  std::_Exit(0);
  // It feels overly brutal to do this, but deallocating memory takes hours, and it reclaimed when the program exits
  // anyway. Other things attempted below failed, and took hours to finish.
  //
  // std::exit(0);

  // gROOT->Reset();
  // return;

  // // Clean up
  // for (int thread = 0; thread < numThreads; ++thread) {
  //   for (auto &hist_pair : hist_map_vec[thread]) {
  //     for (int i_trk_cut = 0; i_trk_cut < nHitCuts; ++i_trk_cut) {
  //       for (int plane = 0; plane < N_PLANES; ++plane) {
  //         for (int bar = 0; bar < N_PADDLES; ++bar) {
  //           delete hist_pair.second[i_trk_cut][plane][bar];
  //         }
  //       }
  //     }
  //   }
  // }

  // cout << "Done!" << endl;
  // // Close the output file

  // return;
}
