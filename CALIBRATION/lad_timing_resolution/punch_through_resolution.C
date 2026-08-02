#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <vector>

// Function to calculate FWHM
double CalculateFWHM(TH1D *hist) {
  int maxBin     = hist->GetMaximumBin();
  double maxVal  = hist->GetBinContent(maxBin);
  double halfMax = maxVal / 2.0;

  // Find left edge
  int leftBin = maxBin;
  for (int i = maxBin; i >= 1; i--) {
    if (hist->GetBinContent(i) < halfMax) {
      leftBin = i;
      break;
    }
  }

  // Find right edge
  int rightBin = maxBin;
  for (int i = maxBin; i <= hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i) < halfMax) {
      rightBin = i;
      break;
    }
  }

  double leftEdge  = hist->GetBinCenter(leftBin);
  double rightEdge = hist->GetBinCenter(rightBin);
  return rightEdge - leftEdge;
}

void punch_through_resolution() {
  // Set batch mode to prevent X11 forwarding
  gROOT->SetBatch(kTRUE);

  // Create TChain and add all ROOT files
  TChain *T            = new TChain("T");
  std::string basePath = "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/";

  T->Add((basePath + "LAD_COIN_23744_6_6_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23747_0_0_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23740_0_0_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23745_5_5_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23746_7_7_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23745_4_4_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23744_1_1_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23745_6_6_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23730_16_16_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23744_8_8_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23743_0_0_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23744_0_0_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23746_12_12_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23745_10_10_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23746_18_18_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23731_0_0_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23744_2_2_-1.root").c_str());
  T->Add((basePath + "LAD_COIN_23732_0_0_-1.root").c_str());

  if (T->GetNtrees() == 0) {
    std::cerr << "Error: Could not load any ROOT files into TChain" << std::endl;
    return;
  }

  const int NBARS   = 11;
  const int NPLANES = 2;

  // Create canvases for each plane
  TCanvas *canvases[NPLANES];
  TH1D *histograms[NPLANES][NBARS];

  // Initialize canvases and histograms
  for (int plane = 0; plane < NPLANES; plane++) {
    canvases[plane] = new TCanvas(Form("c_plane_%d", plane), Form("Time Difference - Plane %d", plane), 1200, 800);
    canvases[plane]->Divide(4, 3);

    for (int bar = 0; bar < NBARS; bar++) {
      histograms[plane][bar] =
          new TH1D(Form("h_plane_%d_bar_%d", plane, bar), Form("Plane %d, Bar %d", plane, bar), 100, 0, 15);
      histograms[plane][bar]->SetXTitle("HitTime_1 - HitTime_0 (ns)");
      histograms[plane][bar]->SetYTitle("Counts");
    }
  }

  // Set branch addresses
  double *hittime_0 = new double[100];
  double *hittime_1 = new double[100];
  double *paddle    = new double[100];
  double *plane     = new double[100];
  int *ndata        = new int;

  T->SetBranchAddress("H.ladhod.goodhit_hittime_0", hittime_0);
  T->SetBranchAddress("H.ladhod.goodhit_hittime_1", hittime_1);
  T->SetBranchAddress("H.ladhod.goodhit_paddle_0", paddle);
  T->SetBranchAddress("H.ladhod.goodhit_plane_0", plane);
  T->SetBranchAddress("Ndata.H.ladhod.goodhit_hittime_0", ndata);

  // Loop through all entries
  Long64_t nentries = T->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    // Display progress
    if (i % 1000 == 0 && i > 0) {
      int progress = (int)((double)i / nentries * 100.0);
      std::cout << "\rProgress: " << progress << "%" << std::flush;
    }

    T->GetEntry(i);

    int n = *ndata;
    for (int j = 0; j < n && j < 100; j++) {
      double diff = hittime_1[j] - hittime_0[j];

      // Skip if difference is exactly zero
      if (diff == 0.0)
        continue;

      int bar_idx   = (int)paddle[j];
      int plane_idx = (int)plane[j] / 2;

      // Make sure indices are within range
      if (bar_idx >= 0 && bar_idx < NBARS && plane_idx >= 0 && plane_idx < NPLANES) {
        histograms[plane_idx][bar_idx]->Fill(diff);
      }
    }
  }

  // Draw histograms and calculate FWHM
  for (int p = 0; p < NPLANES; p++) {
    canvases[p]->cd();
    for (int bar = 0; bar < NBARS; bar++) {
      canvases[p]->cd(bar + 1);
      histograms[p][bar]->Draw();

      // Calculate and display FWHM
      double fwhm = CalculateFWHM(histograms[p][bar]);

      // Add text with FWHM and SIGMA
      TLatex *tex1 = new TLatex(0.6, 0.7, Form("FWHM: %.2f ns", fwhm));
      tex1->SetNDC();
      tex1->SetTextSize(0.05);
      tex1->Draw();

      TLatex *tex2 = new TLatex(0.6, 0.6, Form("SIGMA: %.2f ns", fwhm / 2.355));
      tex2->SetNDC();
      tex2->SetTextSize(0.05);
      tex2->Draw();
    }
  }

  // Write canvases to output file
  TFile *outputFile = new TFile("punch_through_resolution_output.root", "RECREATE");
  for (int p = 0; p < NPLANES; p++) {
    canvases[p]->Write();
  }
  outputFile->Close();

  // Clean up
  delete[] hittime_0;
  delete[] hittime_1;
  delete[] paddle;
  delete[] plane;
  delete ndata;
  delete T;

  std::cout << "\nAnalysis complete!" << std::endl;
}
