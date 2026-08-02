// Lucas Ehinger
// General LAD hodo plotting script
// Single-threaded version
// Not great when it comes to memory management (requires lots of memory). No errors or major memory leaks, just
// duplication of memory across threads since root doesn't like mutexes. Multi-threading with ROOT (which is inherently
// not thread-safe) is hard, and I wasn't smart enough to do it elegantly, but this works.
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const int MAX_DATA   = 50;
const int CACHE_SIZE = 100 * 1024 * 1024; // 100 MB cache size

struct hist_params {
  int NBINS;
  double MIN;
  double MAX;
};

const hist_params time_params      = {70, 1725.0, 1825.0};
const hist_params tof_params       = {140, 20, 150};
const hist_params tof_per_m_params = {200, -10, 10};
const hist_params tof_params_full  = {100, -50, 50};

const double TDC2NS = 0.09766; // TDC to ns conversion factor
const double ADC2NS = 0.0625;  // ADC to ns conversion factor

const int MINT_EVTS_PER_THREAD = 5000000;
bool process_by_file           = true;

const double edep_cut     = 70;  // mV
const double electron_cut = 0.8; // normalized electron energy cut for cal.etottracknorm

const int N_PLANES                 = 5;
const double hodo_angles[N_PLANES] = {150, 150, 127, 127, 104};               // Angles for each plane in deg
const double hodo_radii[N_PLANES]  = {615, 615 + 40.6, 523, 523 + 40.6, 615}; // Radii for each plane in cm

struct hit_cut {
  int n_bar_tolerance;
  bool require_hit;
  bool reject_antihit;
  string cut_name;
  bool include_comp_plt;
};
const int MIN_TDC_TIME = 1600;
const int MAX_TDC_TIME = 2000;

const int nHitCuts         = 7;
hit_cut hit_cuts[nHitCuts] = {{5, false, false, "All_Hits", true},
                              {0, true, false, "Matching_Hit_Tol_0", true},
                              {1, true, false, "Matching_Hit_Tol_1", false},
                              {2, true, false, "Matching_Hit_Tol_2", false},
                              {0, false, true, "Anti-Matching_Hit_Tol_0", true},
                              {1, false, true, "Anti-Matching_Hit_Tol_1", false},
                              {2, false, true, "Anti-Matching_Hit_Tol_2", false}};
const char spec_prefix     = 'P'; // Spectrometer to replay

const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

std::vector<TString> get_file_names(const std::string &filename) {
  std::vector<TString> fileNames;
  std::ifstream infile(filename);
  std::string line;
  while (std::getline(infile, line)) {
    // Skip empty lines and lines starting with '#'
    if (line.empty() || line[0] == '#')
      continue;
    fileNames.push_back(TString(line));
  }
  return fileNames;
}

template <typename T> void add_branch(TTree *tree, const char *branch_name, T *branch_data) {
  // Add a branch to the tree
  tree->SetBranchStatus(branch_name, 1); // Enable the branch
  tree->SetBranchAddress(branch_name, branch_data);
  tree->AddBranchToCache(branch_name, kTRUE);
}

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
  include_comp_plt = false;
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

  if (!include_comp_plt) {
    // Clean up
    delete hist_final_sum;
    for (int plane = 0; plane < N_PLANES; ++plane) {
      delete hist_plane_sum[plane];
    }
    return;
  }

  // Create a directory one level above 'dir' (i.e., sibling to 'dir')
  TString parentDir = dir;
  Ssiz_t slashPos   = parentDir.Last('/');
  if (slashPos != kNPOS) {
    parentDir.Remove(slashPos); // Remove last component to get parent
  } else {
    parentDir = ""; // No parent, use root
  }
  TString compDir = parentDir.IsNull() ? "COMP" : parentDir + "/COMP";
  if (!file->GetDirectory(compDir)) {
    file->mkdir(compDir);
  }
  file->cd(compDir);

  // For each plane, create a comparison canvas and draw all bar histograms
  for (int plane = 0; plane < N_PLANES; ++plane) {
    TCanvas *c_comp_plane = nullptr;
    TString cname         = Form("c_comp_plane_%s", plane_names[plane].c_str());
    if (gROOT->FindObject(cname)) {
      c_comp_plane = (TCanvas *)gROOT->FindObject(cname);
    } else {
      c_comp_plane = new TCanvas(cname, Form("Comparison Plane %s", plane_names[plane].c_str()), 1200, 800);
    }
    c_comp_plane->cd();
    // Only draw the sum over all bars for this plane
    // Get the number of objects already on the canvas
    hist_plane_sum[plane]->SetLineWidth(3);
    // Only draw if nothing is already on the canvas
    if (c_comp_plane->GetListOfPrimitives()->GetSize() == 0) {
      hist_plane_sum[plane]->Draw("HIST");
    } else {
      hist_plane_sum[plane]->Draw("HIST SAME");
    }

    // Adjust y-axis limits for all histograms on the canvas by grabbing them from the canvas
    double max_y = 0;
    double min_y = 0;
    int nhists   = c_comp_plane->GetListOfPrimitives()->GetSize();
    // Find the maximum y value among all histograms on the canvas
    for (int i = 0; i < nhists; ++i) {
      TObject *obj = c_comp_plane->GetListOfPrimitives()->At(i);
      if (obj && obj->InheritsFrom(TH1::Class())) {
        double this_max = ((TH1 *)obj)->GetMaximum();
        double this_min = ((TH1 *)obj)->GetMinimum();
        if (this_max > max_y)
          max_y = this_max;
        if (this_min < min_y)
          min_y = this_min;
      }
    }
    max_y *= 1.1; // Add some padding
    min_y *= 0.9; // Add some padding

    // Set the maximum for all histograms and redraw them
    for (int i = 0; i < nhists; ++i) {
      TObject *obj = c_comp_plane->GetListOfPrimitives()->At(i);
      if (obj && obj->InheritsFrom(TH1::Class())) {
        ((TH1 *)obj)->SetMaximum(max_y);
        ((TH1 *)obj)->SetMinimum(min_y);
        ((TH1 *)obj)->SetLineColor(i + 1); // Assign different colors for each histogram
        ((TH1 *)obj)->SetLineWidth(2);
        ((TH1 *)obj)->Draw(i == 0 ? "HIST" : "HIST SAME");
      }
    }
    // Remove any old legend
    TList *prims = c_comp_plane->GetListOfPrimitives();
    for (int i = prims->GetSize() - 1; i >= 0; --i) {
      TObject *obj = prims->At(i);
      if (obj && obj->InheritsFrom(TLegend::Class())) {
        prims->Remove(obj);
        delete obj;
      }
    }
    // Create and draw a new legend
    TLegend *legend = make_legend(c_comp_plane);
    legend->Draw();
    // Before writing, check if a canvas with the same name exists in the directory and delete it
    TObject *oldCanvas = gDirectory->Get(c_comp_plane->GetName());
    if (oldCanvas && oldCanvas != c_comp_plane) {
      gDirectory->Delete(Form("%s;*", c_comp_plane->GetName())); // Remove all versions with this name
    }
    c_comp_plane->Write();
    delete c_comp_plane;
    delete legend;
  }

  TCanvas *c_comp_all = nullptr;
  TString c_name_all  = "c_COMP_all_planes";
  if (gROOT->FindObject(c_name_all)) {
    c_comp_all = (TCanvas *)gROOT->FindObject(c_name_all);
  } else {
    c_comp_all = new TCanvas(c_name_all, Form("COMP All Planes"), 1200, 800);
  }
  c_comp_all->cd();
  hist_final_sum->SetLineWidth(3);
  // Only draw if nothing is already on the canvas
  if (c_comp_all->GetListOfPrimitives()->GetSize() == 0) {
    hist_final_sum->Draw("HIST");
  } else {
    hist_final_sum->Draw("HIST SAME");
  }
  // Adjust y-axis limits for all histograms on the canvas by grabbing them from the canvas
  double max_y = 0;
  double min_y = 0;
  int nhists   = c_comp_all->GetListOfPrimitives()->GetSize();
  // Find the maximum y value among all histograms on the canvas
  for (int i = 0; i < nhists; ++i) {
    TObject *obj = c_comp_all->GetListOfPrimitives()->At(i);
    if (obj && obj->InheritsFrom(TH1::Class())) {
      double this_max = ((TH1 *)obj)->GetMaximum();
      double this_min = ((TH1 *)obj)->GetMinimum();
      if (this_max > max_y)
        max_y = this_max;
      if (this_min < min_y)
        min_y = this_min;
    }
  }
  max_y *= 1.1; // Add some padding
  min_y *= 0.9; // Add some padding
  // Set the maximum for all histograms and redraw them
  for (int i = 0; i < nhists; ++i) {
    TObject *obj = c_comp_all->GetListOfPrimitives()->At(i);
    if (obj && obj->InheritsFrom(TH1::Class())) {
      ((TH1 *)obj)->SetMaximum(max_y);
      ((TH1 *)obj)->SetMinimum(min_y);
      ((TH1 *)obj)->SetLineColor(i + 1); // Assign different colors for each histogram
      ((TH1 *)obj)->SetLineWidth(2);
      ((TH1 *)obj)->Draw(i == 0 ? "HIST" : "HIST SAME");
    }
  }
  // Remove any old legend
  TList *prims = c_comp_all->GetListOfPrimitives();
  for (int i = prims->GetSize() - 1; i >= 0; --i) {
    TObject *obj = prims->At(i);
    if (obj && obj->InheritsFrom(TLegend::Class())) {
      prims->Remove(obj);
      delete obj;
    }
  }
  // Create and draw a new legend
  TLegend *legend = make_legend(c_comp_all);
  legend->Draw();

  // Before writing, check if a canvas with the same name exists in the directory and delete it
  TObject *oldCanvasAll = gDirectory->Get(c_comp_all->GetName());
  if (oldCanvasAll && oldCanvasAll != c_comp_all) {
    gDirectory->Delete(Form("%s;*", c_comp_all->GetName())); // Remove all versions with this name
  }
  c_comp_all->Write();

  // Clean up
  delete hist_final_sum;
  delete legend;
  delete c_comp_all;
  for (int plane = 0; plane < N_PLANES; ++plane) {
    delete hist_plane_sum[plane];
  }
}

void make_matching_hit_histo(int fullhit_n[N_PLANES], Double_t fullhit_paddle[N_PLANES][MAX_DATA],
                             bool has_hodo_hit[N_PLANES][N_PADDLES], Double_t hit_times[N_PLANES][MAX_DATA],
                             TH1 *hist[nHitCuts][N_PLANES][N_PADDLES], int min_time, int max_time) {
  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      for (int i_hit_cut = 0; i_hit_cut < nHitCuts; i_hit_cut++) {
        // Avg hits
        for (int i_hit = 0; i_hit < fullhit_n[plane]; ++i_hit) {
          if (hit_times[plane][i_hit] < min_time || hit_times[plane][i_hit] > max_time)
            continue; // Skip hits outside the time range

          int bar         = int(fullhit_paddle[plane][i_hit]) - 1;
          int match_plane = -1;
          switch (plane) {
          case 0:
            match_plane = 1;
            break;
          case 1:
            match_plane = 0;
            break;
          case 2:
            match_plane = 3;
            break;
          case 3:
            match_plane = 2;
            break;
          case 4:
            break;
          default:
            break;
          }
          if (plane == 4 && (hit_cuts[i_hit_cut].require_hit || hit_cuts[i_hit_cut].reject_antihit))
            continue; // Skip plane 4 if require_hit or reject_antihit is set

          if (bar < 0 || bar >= N_PADDLES)
            continue; // Skip invalid bars
          bool is_good_hit = true;

          if (hit_cuts[i_hit_cut].require_hit) {
            is_good_hit = false;
            for (int i_bar = bar - hit_cuts[i_hit_cut].n_bar_tolerance;
                 i_bar <= bar + hit_cuts[i_hit_cut].n_bar_tolerance; ++i_bar) {
              if (i_bar < 0 || i_bar >= N_PADDLES)
                continue; // Skip invalid bars
              if (has_hodo_hit[match_plane][i_bar]) {
                is_good_hit = true;
                break;
              }
            }
          }

          if (hit_cuts[i_hit_cut].reject_antihit) {
            if (!hit_cuts[i_hit_cut].require_hit)
              is_good_hit = true;

            for (int i_bar = bar - hit_cuts[i_hit_cut].n_bar_tolerance;
                 i_bar <= bar + hit_cuts[i_hit_cut].n_bar_tolerance; ++i_bar) {
              if (i_bar < 0 || i_bar >= N_PADDLES)
                continue; // Skip invalid bars
              if (has_hodo_hit[match_plane][i_bar]) {
                is_good_hit = false;
                break;
              }
            }
          }

          if (is_good_hit) {
            hist[i_hit_cut][plane][bar]->Fill(hit_times[plane][i_hit]);
          }
        }
      }
    }
  }
}

int lad_tof_plots_single() {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);

  // Open multiple ROOT files
  // std::vector<TString> fileNames = {
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22591_0_2_-1.root"};

  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22572_0_6_2000000.root"};
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22382_0_21_-1.root",
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22382_0_21_-1_1.root"};
  std::vector<TString> fileNames = get_file_names("../files/run-lists/all_C3_runlist_22745-23590.dat");
  // std::vector<TString> fileNames = get_file_names("../files/run-lists/all_C3_runlist_23700.dat");
  // std::vector<TString> fileNames = get_file_names("../files/run-lists/all_C3_runlist.dat");
  // std::vector<TString> fileNames = get_file_names("files/LD2_setting1_new2.dat");

  TString outputFileName = Form("files/root/updated_timing_plots_single_C3_22745-23590_%c.root", spec_prefix);
  // TString outputFileName =
  //     Form("files/raw_hodo_timing_plots/raw_hodo_timing_plots_LD2_setting1_new2_%c.root", spec_prefix);

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
  // nEntries     = 10000;

  // Number of threads to use
  int numThreads = 1;

  std::vector<map<string, TH1 *[nHitCuts][N_PLANES][N_PADDLES]>> hist_map_vec(numThreads);

  for (int thread = 0; thread < numThreads; ++thread) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        for (int i_hit_cut = 0; i_hit_cut < nHitCuts; i_hit_cut++) {
          hist_map_vec[thread]["h_TDC_top_bar"][i_hit_cut][plane][bar] = new TH1F(
              Form("h_TDC_top_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                   hit_cuts[i_hit_cut].cut_name.c_str()),
              Form("TDC Top Plane %s Bar %d %s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_TDC_top_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_TDC_top_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_TDC_btm_bar"][i_hit_cut][plane][bar] = new TH1F(
              Form("h_TDC_btm_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                   hit_cuts[i_hit_cut].cut_name.c_str()),
              Form("TDC Btm Plane %s Bar %d %s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_TDC_btm_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_TDC_btm_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ADC_top_bar"][i_hit_cut][plane][bar] = new TH1F(
              Form("h_ADC_top_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                   hit_cuts[i_hit_cut].cut_name.c_str()),
              Form("ADC Top Plane %s Bar %d %s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_ADC_top_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ADC Time (ns)");
          hist_map_vec[thread]["h_ADC_top_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ADC_btm_bar"][i_hit_cut][plane][bar] = new TH1F(
              Form("h_ADC_btm_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                   hit_cuts[i_hit_cut].cut_name.c_str()),
              Form("ADC Btm Plane %s Bar %d %s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_ADC_btm_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ADC Time (ns)");
          hist_map_vec[thread]["h_ADC_btm_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_Time_Avg_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_Time_Avg_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time Avg Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_Time_Avg_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("Time Avg (ns)");
          hist_map_vec[thread]["h_Time_Avg_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_bar"][i_hit_cut][plane][bar] = new TH1F(
              Form("h_ToF_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              Form("ToF Plane %s Bar %d %s", plane_names[plane].c_str(), bar, hit_cuts[i_hit_cut].cut_name.c_str()),
              tof_params.NBINS, tof_params.MIN, tof_params.MAX);
          hist_map_vec[thread]["h_ToF_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF (ns)");
          hist_map_vec[thread]["h_ToF_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_per_m_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_per_m_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF per m Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_per_m_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF per m (ns/m)");
          hist_map_vec[thread]["h_ToF_per_m_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // ToF Vertex per m
          hist_map_vec[thread]["h_ToF_vertex_per_m_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_per_m_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex per m Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_per_m_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_per_m_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // ToF Vertex RFcorr per m
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_per_m_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr per m Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_trigger_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_trigger_bar_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Trigger Bar Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_trigger_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF Trigger (ns)");
          hist_map_vec[thread]["h_ToF_trigger_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // ToF relative to vertex (raw)
          hist_map_vec[thread]["h_ToF_vertex_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_bar_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF Vertex (ns)");
          hist_map_vec[thread]["h_ToF_vertex_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // ToF relative to vertex (RF-corrected)
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_bar_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // Versions with edepCut before _bar
          hist_map_vec[thread]["h_Time_Avg_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_Time_Avg_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time Avg edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_Time_Avg_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("Time Avg (ns)");
          hist_map_vec[thread]["h_Time_Avg_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params.NBINS, tof_params.MIN, tof_params.MAX);
          hist_map_vec[thread]["h_ToF_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF (ns)");
          hist_map_vec[thread]["h_ToF_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_trigger_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_trigger_edepCut_bar_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Trigger edepCut Bar Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_trigger_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Trigger (ns)");
          hist_map_vec[thread]["h_ToF_trigger_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // edepCut versions for ToF Vertex
          hist_map_vec[thread]["h_ToF_vertex_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex (ns)");
          hist_map_vec[thread]["h_ToF_vertex_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          // edepCut versions for ToF Vertex per m
          hist_map_vec[thread]["h_ToF_vertex_per_m_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_per_m_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex per m edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_per_m_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_per_m_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_per_m_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr per m edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_edepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_edepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          hist_map_vec[thread]["h_Time_Avg_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_Time_Avg_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time Avg antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_Time_Avg_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "Time Avg (ns)");
          hist_map_vec[thread]["h_Time_Avg_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params.NBINS, tof_params.MIN, tof_params.MAX);
          hist_map_vec[thread]["h_ToF_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle("ToF (ns)");
          hist_map_vec[thread]["h_ToF_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_trigger_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_trigger_antiEdepCut_bar_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Trigger antiEdepCut Bar Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_trigger_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Trigger (ns)");
          hist_map_vec[thread]["h_ToF_trigger_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          // antiEdepCut versions for ToF Vertex per m
          hist_map_vec[thread]["h_ToF_vertex_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_per_m_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex per m antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_per_m_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr per m antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_per_m_params.NBINS, tof_per_m_params.MIN, tof_per_m_params.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar]
              ->GetXaxis()
              ->SetTitle("ToF Vertex RFcorr per m (ns/m)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_per_m_antiEdepCut_bar"][i_hit_cut][plane][bar]
              ->GetYaxis()
              ->SetTitle("Counts");

          // antiEdepCut versions for ToF Vertex
          hist_map_vec[thread]["h_ToF_vertex_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex (ns)");
          hist_map_vec[thread]["h_ToF_vertex_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_antiEdepCut_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          // photonCorr versions for ToF Vertex RFcorr
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_photonCorr_plane_%s_bar_%d_%s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr photonCorr Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_bar"][i_hit_cut][plane][bar]->GetXaxis()->SetTitle(
              "ToF Vertex RFcorr photonCorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_bar"][i_hit_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_edepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_photonCorr_edepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(),
                            bar, hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr photonCorr edepCut Plane %s Bar %d %s", plane_names[plane].c_str(), bar,
                            hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_edepCut_bar"][i_hit_cut][plane][bar]
              ->GetXaxis()
              ->SetTitle("ToF Vertex RFcorr photonCorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_edepCut_bar"][i_hit_cut][plane][bar]
              ->GetYaxis()
              ->SetTitle("Counts");

          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_bar"][i_hit_cut][plane][bar] =
              new TH1F(Form("h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_plane_%s_bar_%d_%s", plane_names[plane].c_str(),
                            bar, hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF Vertex RFcorr photonCorr antiEdepCut Plane %s Bar %d %s", plane_names[plane].c_str(),
                            bar, hit_cuts[i_hit_cut].cut_name.c_str()),
                       tof_params_full.NBINS, tof_params_full.MIN, tof_params_full.MAX);
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_bar"][i_hit_cut][plane][bar]
              ->GetXaxis()
              ->SetTitle("ToF Vertex RFcorr photonCorr (ns)");
          hist_map_vec[thread]["h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_bar"][i_hit_cut][plane][bar]
              ->GetYaxis()
              ->SetTitle("Counts");
        }
      }
    }
  }

  // Create an output ROOT file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return -1;
  }

  // Helper lambda to format numbers with commas
  auto format_with_commas = [](int64_t value) {
    std::string num    = std::to_string(value);
    int insertPosition = num.length() - 3;
    while (insertPosition > 0) {
      num.insert(insertPosition, ",");
      insertPosition -= 3;
    }
    return num;
  };

  // Process all events in one thread
  cout << "Starting single-threaded processing for " << format_with_commas(nEntries) << " events." << endl;

  // Hodo-level data
  Double_t tdc_time[N_SIDES][N_PLANES][MAX_DATA];
  Double_t tdc_counter[N_SIDES][N_PLANES][MAX_DATA];
  Double_t tdc_time_UnCorr[N_SIDES][N_PLANES][MAX_DATA], tdc_time_TWCorr[N_SIDES][N_PLANES][MAX_DATA],
      tdc_time_Corr[N_SIDES][N_PLANES][MAX_DATA];
  Double_t fullhit_time_avg[N_PLANES][MAX_DATA], fullhit_adc_avg[N_PLANES][MAX_DATA],
      fullhit_tof_avg[N_PLANES][MAX_DATA], fullhit_paddle[N_PLANES][MAX_DATA], fullhit_y_pos[N_PLANES][MAX_DATA];
  Int_t fullhit_n[N_PLANES];

  Double_t adc_time[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_counter[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_time_good[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_amp[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_amp_good[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_int[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_int_good[N_SIDES][N_PLANES][MAX_DATA];

  Int_t nData_adc[N_SIDES][N_PLANES];
  Int_t nData_tdc[N_SIDES][N_PLANES];
  Double_t hodo_start_time, pTRIG1, pTRIG2, pTRIG3, pTRIG4, evtyp, vertex_time, vertex_time_RFcorr, etottracknorm;

  T->SetCacheSize(CACHE_SIZE);
  T->SetBranchStatus("*", 0); // Disable all branches initially

  add_branch(T, Form("%c.hod.starttime", spec_prefix), &hodo_start_time);
  add_branch(T, "T.hms.pTRIG1_tdcTime", &pTRIG1);
  add_branch(T, "T.hms.pTRIG2_tdcTime", &pTRIG2);
  add_branch(T, "T.hms.pTRIG3_tdcTime", &pTRIG3);
  add_branch(T, "T.hms.pTRIG4_tdcTime", &pTRIG4);
  add_branch(T, "g.evtyp", &evtyp);
  add_branch(T, Form("%c.ladkin.t_vertex", spec_prefix), &vertex_time);
  add_branch(T, Form("%c.ladkin.t_vertex_RFcorr", spec_prefix), &vertex_time_RFcorr);
  // Normalized total track energy for electron selection
  add_branch(T, Form("%c.cal.etottracknorm", spec_prefix), &etottracknorm);

  for (int side = 0; side < N_SIDES; ++side) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      add_branch(T, Form("%c.ladhod.%s.%sTdcTime", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &tdc_time[side][plane]);
      add_branch(
          T,
          Form("%c.ladhod.%s.Good%sTdcTimeUnCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_UnCorr[side][plane]);
      add_branch(
          T,
          Form("%c.ladhod.%s.Good%sTdcTimeWalkCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_TWCorr[side][plane]);
      add_branch(
          T, Form("%c.ladhod.%s.Good%sTdcTimeCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_Corr[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcPulseTime", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_time[side][plane]);
      add_branch(
          T, Form("%c.ladhod.%s.Good%sAdcPulseTime", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_time_good[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcPulseAmp", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_amp[side][plane]);
      add_branch(
          T, Form("%c.ladhod.%s.Good%sAdcPulseAmp", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_amp_good[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcPulseInt", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_int[side][plane]);
      add_branch(
          T, Form("%c.ladhod.%s.Good%sAdcPulseInt", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_int_good[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &tdc_counter[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_counter[side][plane]);
      add_branch(
          T, Form("Ndata.%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_adc[side][plane]);
      add_branch(
          T, Form("Ndata.%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_tdc[side][plane]);
    }
  }
  for (int plane = 0; plane < N_PLANES; ++plane) {
    add_branch(T, Form("%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()), &fullhit_time_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitTOF", spec_prefix, plane_names[plane].c_str()), &fullhit_tof_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitEdepAmp", spec_prefix, plane_names[plane].c_str()),
               &fullhit_adc_avg[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitPaddleNum", spec_prefix, plane_names[plane].c_str()),
               &fullhit_paddle[plane]);
    add_branch(T, Form("Ndata.%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()), &fullhit_n[plane]);
    add_branch(T, Form("%c.ladhod.%s.HodoHitPos", spec_prefix, plane_names[plane].c_str()), &fullhit_y_pos[plane]);
  }

  Bool_t has_hodo_hit_avg[N_PLANES][N_PADDLES];
  Bool_t has_hodo_hit_tdc[N_SIDES][N_PLANES][N_PADDLES];
  Bool_t has_hodo_hit_adc[N_SIDES][N_PLANES][N_PADDLES];
  ////////////////////////////////////////////////////
  // Start Event Loop

  int start = 0;
  int end   = nEntries;
  if (end < 1) {
    end = T->GetEntries(); // If end is less than 1, set it to the total number of entries
  }

  for (int i = start; i < end; ++i) {
    T->GetEntry(i);
    if ((spec_prefix == 'H' && int(evtyp) != 2) || (spec_prefix == 'P' && int(evtyp) != 1))
      continue; // Only process events with the correct evtyp based on spec_prefix

    // Require electron-like event based on normalized calorimeter energy
    if (etottracknorm <= electron_cut)
      continue;

    // Loop through tdc_time, and set all entries to themselves times TDC2NS
    for (int side = 0; side < N_SIDES; ++side) {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int i_hit = 0; i_hit < nData_tdc[side][plane]; ++i_hit) {
          tdc_time[side][plane][i_hit] = (tdc_time[side][plane][i_hit] + 31700) * TDC2NS;
        }
      }
    }

    // clear has_hodo_hit array
    for (int side = 0; side < N_SIDES; ++side) {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int bar = 0; bar < N_PADDLES; ++bar) {
          if (side == 0)
            has_hodo_hit_avg[plane][bar] = false;
          has_hodo_hit_tdc[side][plane][bar] = false;
          has_hodo_hit_adc[side][plane][bar] = false;
        }
      }
    }
    // Loop over the hits and fill the arrays
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int i_top = 0; i_top < fullhit_n[plane]; ++i_top) {
        int bar = int(fullhit_paddle[plane][i_top]) - 1;
        if (bar < 0 || bar >= N_PADDLES)
          continue; // Skip invalid bars
        has_hodo_hit_avg[plane][bar] = true;
      }
    }
    for (int side = 0; side < N_SIDES; ++side) {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int i_hit = 0; i_hit < nData_tdc[side][plane]; ++i_hit) {
          int bar = int(tdc_counter[side][plane][i_hit]) - 1;
          if (bar < 0 || bar >= N_PADDLES)
            continue; // Skip invalid bars
          if (tdc_time[side][plane][i_hit] < MIN_TDC_TIME || tdc_time[side][plane][i_hit] > MAX_TDC_TIME)
            continue; // Skip invalid TDC times
          has_hodo_hit_tdc[side][plane][bar] = true;
        }
        for (int i_hit = 0; i_hit < nData_adc[side][plane]; ++i_hit) {
          int bar = int(adc_counter[side][plane][i_hit]) - 1;
          if (bar < 0 || bar >= N_PADDLES)
            continue; // Skip invalid bars
          has_hodo_hit_adc[side][plane][bar] = true;
        }
      }
    }

    // Calculate fullhit_tof_avg by subtracting hodo_start_time from fullhit_time_avg
    Double_t fullhit_tof_trigger_avg[N_PLANES][MAX_DATA] = {-999};
    Double_t fullhit_tof_avg_per_m[N_PLANES][MAX_DATA]   = {-999};
    // New: ToF relative to vertex time (raw and RF-corrected) and per-meter versions
    Double_t fullhit_tof_vertex_avg[N_PLANES][MAX_DATA]               = {-999};
    Double_t fullhit_tof_vertex_RFcorr_avg[N_PLANES][MAX_DATA]        = {-999};
    Double_t fullhit_tof_vertex_per_m[N_PLANES][MAX_DATA]             = {-999};
    Double_t fullhit_tof_vertex_RFcorr_per_m[N_PLANES][MAX_DATA]      = {-999};
    Double_t fullhit_tof_vertex_RFcorr_photonCorr[N_PLANES][MAX_DATA] = {-999};
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int i_hit = 0; i_hit < fullhit_n[plane]; ++i_hit) {
        double path_length =
            sqrt(pow(fullhit_y_pos[plane][i_hit], 2) + pow(22 * (fullhit_paddle[plane][i_hit] - 6), 2));
        path_length = sqrt(path_length * path_length + hodo_radii[plane] * hodo_radii[plane]) / 100;
        fullhit_tof_trigger_avg[plane][i_hit] =
            (fullhit_time_avg[plane][i_hit] - hodo_start_time + 60 - 1750 + 19 - 2 * 6) / path_length;
        // ToF relative to vertex (simple subtraction as requested)
        fullhit_tof_vertex_avg[plane][i_hit] = (fullhit_time_avg[plane][i_hit] - vertex_time - 10 - 1750 + 19 - 2 * 6);
        fullhit_tof_vertex_RFcorr_avg[plane][i_hit] =
            (fullhit_time_avg[plane][i_hit] - vertex_time_RFcorr - 10 - 1750 + 19 - 2 * 6);
        // (fullhit_time_avg[plane][i_hit] - (2 * vertex_time - vertex_time_RFcorr) - 10 - 1750 + 19 - 2 * 6);
        // Per-meter versions (divide by path length)
        if (path_length != 0) {
          fullhit_tof_vertex_per_m[plane][i_hit]        = fullhit_tof_vertex_avg[plane][i_hit] / path_length;
          fullhit_tof_vertex_RFcorr_per_m[plane][i_hit] = fullhit_tof_vertex_RFcorr_avg[plane][i_hit] / path_length;
          fullhit_tof_avg_per_m[plane][i_hit]           = (fullhit_tof_avg[plane][i_hit] - 40) / path_length;
          // Photon correction: subtract (path_length - hodo_radii[plane])/3 from RFcorr time
          fullhit_tof_vertex_RFcorr_photonCorr[plane][i_hit] =
              fullhit_tof_vertex_RFcorr_avg[plane][i_hit] - (path_length / 100.0) / 0.3;
        } else {
          fullhit_tof_vertex_per_m[plane][i_hit]             = -999;
          fullhit_tof_vertex_RFcorr_per_m[plane][i_hit]      = -999;
          fullhit_tof_avg_per_m[plane][i_hit]                = -999;
          fullhit_tof_vertex_RFcorr_photonCorr[plane][i_hit] = -999;
        }
      }
    }

    Double_t fullhit_time_avg_edepCut[N_PLANES][MAX_DATA]                     = {-999};
    Double_t fullhit_tof_avg_edepCut[N_PLANES][MAX_DATA]                      = {-999};
    Double_t fullhit_tof_trigger_avg_edepCut[N_PLANES][MAX_DATA]              = {-999};
    Double_t fullhit_tof_vertex_avg_edepCut[N_PLANES][MAX_DATA]               = {-999};
    Double_t fullhit_tof_vertex_RFcorr_avg_edepCut[N_PLANES][MAX_DATA]        = {-999};
    Double_t fullhit_tof_vertex_per_m_edepCut[N_PLANES][MAX_DATA]             = {-999};
    Double_t fullhit_tof_vertex_RFcorr_per_m_edepCut[N_PLANES][MAX_DATA]      = {-999};
    Double_t fullhit_tof_vertex_RFcorr_photonCorr_edepCut[N_PLANES][MAX_DATA] = {-999};
    Int_t fullhit_n_edepCut[N_PLANES]                                         = {0};

    Double_t fullhit_time_avg_antiEdepCut[N_PLANES][MAX_DATA]                     = {-999};
    Double_t fullhit_tof_avg_antiEdepCut[N_PLANES][MAX_DATA]                      = {-999};
    Double_t fullhit_tof_trigger_avg_antiEdepCut[N_PLANES][MAX_DATA]              = {-999};
    Double_t fullhit_tof_vertex_avg_antiEdepCut[N_PLANES][MAX_DATA]               = {-999};
    Double_t fullhit_tof_vertex_RFcorr_avg_antiEdepCut[N_PLANES][MAX_DATA]        = {-999};
    Double_t fullhit_tof_vertex_per_m_antiEdepCut[N_PLANES][MAX_DATA]             = {-999};
    Double_t fullhit_tof_vertex_RFcorr_per_m_antiEdepCut[N_PLANES][MAX_DATA]      = {-999};
    Double_t fullhit_tof_vertex_RFcorr_photonCorr_antiEdepCut[N_PLANES][MAX_DATA] = {-999};
    Int_t fullhit_n_antiEdepCut[N_PLANES]                                         = {0};

    for (int plane = 0; plane < N_PLANES; ++plane) {
      fullhit_n_edepCut[plane]     = 0;
      fullhit_n_antiEdepCut[plane] = 0;
      for (int i_hit = 0; i_hit < fullhit_n[plane]; ++i_hit) {
        if (fullhit_adc_avg[plane][i_hit] > edep_cut) {
          int idx                                                  = fullhit_n_edepCut[plane];
          fullhit_time_avg_edepCut[plane][idx]                     = fullhit_time_avg[plane][i_hit];
          fullhit_tof_avg_edepCut[plane][idx]                      = fullhit_tof_avg[plane][i_hit];
          fullhit_tof_trigger_avg_edepCut[plane][idx]              = fullhit_tof_trigger_avg[plane][i_hit];
          fullhit_tof_vertex_avg_edepCut[plane][idx]               = fullhit_tof_vertex_avg[plane][i_hit];
          fullhit_tof_vertex_RFcorr_avg_edepCut[plane][idx]        = fullhit_tof_vertex_RFcorr_avg[plane][i_hit];
          fullhit_tof_vertex_per_m_edepCut[plane][idx]             = fullhit_tof_vertex_per_m[plane][i_hit];
          fullhit_tof_vertex_RFcorr_per_m_edepCut[plane][idx]      = fullhit_tof_vertex_RFcorr_per_m[plane][i_hit];
          fullhit_tof_vertex_RFcorr_photonCorr_edepCut[plane][idx] = fullhit_tof_vertex_RFcorr_photonCorr[plane][i_hit];
          fullhit_n_edepCut[plane]++;
        } else {
          int idx                                                 = fullhit_n_antiEdepCut[plane];
          fullhit_time_avg_antiEdepCut[plane][idx]                = fullhit_time_avg[plane][i_hit];
          fullhit_tof_avg_antiEdepCut[plane][idx]                 = fullhit_tof_avg[plane][i_hit];
          fullhit_tof_trigger_avg_antiEdepCut[plane][idx]         = fullhit_tof_trigger_avg[plane][i_hit];
          fullhit_tof_vertex_avg_antiEdepCut[plane][idx]          = fullhit_tof_vertex_avg[plane][i_hit];
          fullhit_tof_vertex_RFcorr_avg_antiEdepCut[plane][idx]   = fullhit_tof_vertex_RFcorr_avg[plane][i_hit];
          fullhit_tof_vertex_per_m_antiEdepCut[plane][idx]        = fullhit_tof_vertex_per_m[plane][i_hit];
          fullhit_tof_vertex_RFcorr_per_m_antiEdepCut[plane][idx] = fullhit_tof_vertex_RFcorr_per_m[plane][i_hit];
          fullhit_tof_vertex_RFcorr_photonCorr_antiEdepCut[plane][idx] =
              fullhit_tof_vertex_RFcorr_photonCorr[plane][i_hit];
          fullhit_n_antiEdepCut[plane]++;
        }
      }
    }

    // Fill the histograms
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_avg, hist_map_vec[0]["h_ToF_bar"],
                            -200, 500);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_avg_edepCut,
                            hist_map_vec[0]["h_ToF_edepCut_bar"], -200, 500);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_avg_antiEdepCut,
                            hist_map_vec[0]["h_ToF_antiEdepCut_bar"], -200, 500);
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_avg_per_m,
                            hist_map_vec[0]["h_ToF_per_m_bar"], -10, 50);
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_trigger_avg,
                            hist_map_vec[0]["h_ToF_trigger_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_trigger_avg_edepCut,
                            hist_map_vec[0]["h_ToF_trigger_edepCut_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_trigger_avg_antiEdepCut, hist_map_vec[0]["h_ToF_trigger_antiEdepCut_bar"], -100,
                            100);
    // New: ToF relative to vertex and RF-corrected vertex
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_avg,
                            hist_map_vec[0]["h_ToF_vertex_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_avg_edepCut,
                            hist_map_vec[0]["h_ToF_vertex_edepCut_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_avg_antiEdepCut,
                            hist_map_vec[0]["h_ToF_vertex_antiEdepCut_bar"], -100, 100);

    // ToF vertex per-meter versions
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_per_m,
                            hist_map_vec[0]["h_ToF_vertex_per_m_bar"], -10, 10);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_per_m_edepCut,
                            hist_map_vec[0]["h_ToF_vertex_per_m_edepCut_bar"], -10, 10);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_per_m_antiEdepCut, hist_map_vec[0]["h_ToF_vertex_per_m_antiEdepCut_bar"],
                            -10, 10);

    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_RFcorr_avg,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_RFcorr_avg_edepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_edepCut_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_RFcorr_avg_antiEdepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_antiEdepCut_bar"], -100, 100);
    // ToF vertex RFcorr per-meter versions
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_RFcorr_per_m,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_bar"], -10, 10);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_RFcorr_per_m_edepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_edepCut_bar"], -10, 10);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_RFcorr_per_m_antiEdepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_antiEdepCut_bar"], -10, 10);
    // ToF vertex RFcorr photon-corrected versions
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_tof_vertex_RFcorr_photonCorr,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_RFcorr_photonCorr_edepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_edepCut_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg,
                            fullhit_tof_vertex_RFcorr_photonCorr_antiEdepCut,
                            hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_bar"], -100, 100);
    make_matching_hit_histo(fullhit_n, fullhit_paddle, has_hodo_hit_avg, fullhit_time_avg,
                            hist_map_vec[0]["h_Time_Avg_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(fullhit_n_edepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_time_avg_edepCut,
                            hist_map_vec[0]["h_Time_Avg_edepCut_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(fullhit_n_antiEdepCut, fullhit_paddle, has_hodo_hit_avg, fullhit_time_avg_antiEdepCut,
                            hist_map_vec[0]["h_Time_Avg_antiEdepCut_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(nData_tdc[0], tdc_counter[0], has_hodo_hit_tdc[0], tdc_time[0],
                            hist_map_vec[0]["h_TDC_top_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(nData_tdc[1], tdc_counter[1], has_hodo_hit_tdc[1], tdc_time[1],
                            hist_map_vec[0]["h_TDC_btm_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(nData_adc[0], adc_counter[0], has_hodo_hit_adc[0], adc_time[0],
                            hist_map_vec[0]["h_ADC_top_bar"], MIN_TDC_TIME, MAX_TDC_TIME);
    make_matching_hit_histo(nData_adc[1], adc_counter[1], has_hodo_hit_adc[1], adc_time[1],
                            hist_map_vec[0]["h_ADC_btm_bar"], MIN_TDC_TIME, MAX_TDC_TIME);

    // Loop over the hits and fill the histograms
    // Print the status as a percentage
    if (i % ((end - start) / 100) == 0) {
      std::cout << "\rProcessing: " << int((i - start) * 100.0 / (end - start)) << "% completed." << std::flush;
    }
  } // End Event Loop

  std::cout << "\rProcessing: 100% completed. \nWriting histograms" << std::endl;

  cout << "Writing histograms to file..." << endl;
  // Write histograms to the output file
  outputFile->cd();

  // Prepare progress tracking for writing canvases
  int totalCallsPerCut = 29; // number of write calls per hit-cut. Don't love hardcoding this, but is what it is for now
  int totalCalls       = totalCallsPerCut * nHitCuts;
  int globalWriteCount = 0;
  auto write_and_progress = [&](TH1 *histArr[N_PLANES][N_PADDLES], TFile *file, TString dir, TString var_name,
                                bool include_comp_plt) {
    write_to_canvas_plane(histArr, file, dir, var_name, include_comp_plt);
    ++globalWriteCount;
    int pct = int(100.0 * globalWriteCount / totalCalls);
    std::cout << "\rWriting histograms: " << pct << "% completed." << std::flush;
  };

  for (int i_hit_cut = 0; i_hit_cut < nHitCuts; ++i_hit_cut) {
    write_and_progress(hist_map_vec[0]["h_TDC_top_bar"][i_hit_cut], outputFile,
                       Form("KIN/TDC_TOP/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("TDC_TOP_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_TDC_btm_bar"][i_hit_cut], outputFile,
                       Form("KIN/TDC_BTM/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("TDC_BTM_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ADC_top_bar"][i_hit_cut], outputFile,
                       Form("KIN/ADC_TOP/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ADC_TOP_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ADC_btm_bar"][i_hit_cut], outputFile,
                       Form("KIN/ADC_BTM/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ADC_BTM_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_Time_Avg_bar"][i_hit_cut], outputFile,
                       Form("KIN/Time_Avg/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time_Avg_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_%s", hit_cuts[i_hit_cut].cut_name.c_str()), hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_per_m_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_per_m/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_per_m_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_per_m_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_per_m/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_per_m_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_per_m/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_per_m_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_trigger_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Trigger/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Trigger_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_photonCorr/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_photonCorr_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_Time_Avg_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/Time_Avg_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time_Avg_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_trigger_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Trigger_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Trigger_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_photonCorr_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_photonCorr_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_per_m_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_per_m_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_per_m_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_edepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_per_m_edepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_per_m_edepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_Time_Avg_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/Time_Avg_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("Time_Avg_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_trigger_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Trigger_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Trigger_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_photonCorr_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_photonCorr_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_photonCorr_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_per_m_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_per_m_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_per_m_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
    write_and_progress(hist_map_vec[0]["h_ToF_vertex_RFcorr_per_m_antiEdepCut_bar"][i_hit_cut], outputFile,
                       Form("KIN/ToF_Vertex_RFcorr_per_m_antiEdepCut/%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       Form("ToF_Vertex_RFcorr_per_m_antiEdepCut_%s", hit_cuts[i_hit_cut].cut_name.c_str()),
                       hit_cuts[i_hit_cut].include_comp_plt);
  }
  // Ensure progress reaches 100% and end the line
  std::cout << "\rWriting histograms: 100% completed." << std::endl;

  outputFile->Close();
  delete outputFile;

  std::cout << "Done! All histograms written to " << outputFileName << endl << "Cleaning up..." << std::endl;

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