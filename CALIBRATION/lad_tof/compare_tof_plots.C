#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TObjArray.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Each entry pairs a canvas (path inside the ROOT file) with the legend label
// to use when drawing it.
struct PlotEntry {
  std::string canvasName;
  std::string label;
};

// Helper function to create a comparison canvas AND a ratio canvas for an
// arbitrary number of input canvases. Ratios are taken as entries[i] / entries[0]
// for i = 1..N-1, after normalizing each histogram so that the integral of its
// first `ratioNormBins` bins equals 1. Returns {comparisonCanvas, ratioCanvas}.
std::pair<TCanvas *, TCanvas *> createComparisonCanvases(TFile *f, const std::vector<PlotEntry> &entries,
                                                         const char *canvasTitle, const char *ratioTitle,
                                                         int ratioNormBins = 20) {

  const int N = (int)entries.size();
  if (N < 1) {
    std::cerr << "Error: createComparisonCanvases needs at least 1 entry" << std::endl;
    return {nullptr, nullptr};
  }

  // Retrieve the source canvases
  std::vector<TCanvas *> sourceCanvases;
  for (int i = 0; i < N; i++) {
    TCanvas *c = (TCanvas *)f->Get(entries[i].canvasName.c_str());
    if (!c) {
      std::cerr << "Error: Could not find canvas " << entries[i].canvasName << std::endl;
      return {nullptr, nullptr};
    }
    sourceCanvases.push_back(c);
  }

  // Extract histograms from each source canvas (assume up to 6 histograms per canvas)
  std::vector<std::vector<TH1F *>> histograms(N);
  int minHistsAcrossCanvases = 6;
  for (int canvasIdx = 0; canvasIdx < N; canvasIdx++) {
    TList *primitives = sourceCanvases[canvasIdx]->GetListOfPrimitives();
    int histCount     = 0;

    for (int i = 0; i < primitives->GetSize() && histCount < 6; i++) {
      TObject *obj = primitives->At(i);
      TPad *pad    = dynamic_cast<TPad *>(obj);
      if (!pad)
        continue; // Skip if not a pad
      TH1F *hist = dynamic_cast<TH1F *>(pad->GetListOfPrimitives()->At(0));
      if (hist) {
        histograms[canvasIdx].push_back(hist);
        histCount++;
      }
    }

    if (histCount < 6) {
      std::cerr << "Warning: Canvas " << canvasIdx << " (" << entries[canvasIdx].canvasName << ") has only "
                << histCount << " histograms, expected 6" << std::endl;
    }
    if (histCount < minHistsAcrossCanvases)
      minHistsAcrossCanvases = histCount;
  }

  // Create the comparison canvas with 3 rows and 2 columns (6 pads total)
  TCanvas *compareCanvas = new TCanvas(canvasTitle, canvasTitle, 1200, 1000);
  compareCanvas->Divide(2, 3); // 2 columns, 3 rows

  // Create the ratio canvas with the same layout
  TCanvas *ratioCanvas = new TCanvas(ratioTitle, ratioTitle, 1200, 1000);
  ratioCanvas->Divide(2, 3);

  // Default color cycle; falls back to ROOT palette indices for N > 8
  int defaultColors[]      = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kViolet};
  const int nDefaultColors = (int)(sizeof(defaultColors) / sizeof(defaultColors[0]));
  auto colorFor            = [&](int i) -> int {
    if (i < nDefaultColors)
      return defaultColors[i];
    return TColor::GetColorPalette(i % TColor::GetNumberOfColors());
  };

  // For each of the (up to 6) pads
  for (int padIdx = 0; padIdx < 6; padIdx++) {
    compareCanvas->cd(padIdx + 1);

    float maxY = 0;
    std::vector<TH1F *> clonedHists;

    // Find the maximum y value across all input histograms for this pad
    for (int canvasIdx = 0; canvasIdx < N; canvasIdx++) {
      if (padIdx < (int)histograms[canvasIdx].size()) {
        TH1F *hist = histograms[canvasIdx][padIdx];
        float max  = hist->GetMaximum();
        if (max > maxY) {
          maxY = max;
        }
      }
    }

    // Normalize and plot the histograms on the comparison canvas
    for (int canvasIdx = 0; canvasIdx < N; canvasIdx++) {
      if (padIdx < (int)histograms[canvasIdx].size()) {
        TH1F *hist = histograms[canvasIdx][padIdx];

        TH1F *clonedHist = (TH1F *)hist->Clone(Form("%s_clone_pad%d_c%d", canvasTitle, padIdx, canvasIdx));
        clonedHists.push_back(clonedHist);

        // Normalize so that each histogram's max equals maxY
        if (clonedHist->GetMaximum() > 0) {
          clonedHist->Scale(maxY / clonedHist->GetMaximum());
        }

        clonedHist->SetLineColor(colorFor(canvasIdx));
        clonedHist->SetLineWidth(2);

        if (canvasIdx == 0) {
          clonedHist->Draw("HIST");
        } else {
          clonedHist->Draw("HIST SAME");
        }
      }
    }

    // Legend on the comparison pad
    TLegend *leg = new TLegend(0.7, 0.7, 0.95, 0.95);
    leg->SetBorderSize(0);
    for (size_t canvasIdx = 0; canvasIdx < clonedHists.size(); canvasIdx++) {
      leg->AddEntry(clonedHists[canvasIdx], entries[canvasIdx].label.c_str(), "l");
    }
    leg->Draw();

    gPad->SetGrid();

    // -----------------------------
    // Ratio canvas: plot[i]/plot[0] for i = 1..N-1
    // -----------------------------
    ratioCanvas->cd(padIdx + 1);

    if (clonedHists.size() < 2) {
      // Need at least the reference + 1 other to form a ratio
      continue;
    }

    // Build separate copies of the originals for the ratio, normalized by the
    // integral of the first `ratioNormBins` bins (so that integral == 1 there).
    std::vector<TH1F *> ratioNormHists;
    for (int canvasIdx = 0; canvasIdx < N; canvasIdx++) {
      if (padIdx >= (int)histograms[canvasIdx].size()) {
        ratioNormHists.push_back(nullptr);
        continue;
      }
      TH1F *src           = histograms[canvasIdx][padIdx];
      TH1F *nh            = (TH1F *)src->Clone(Form("%s_rationorm_pad%d_c%d", ratioTitle, padIdx, canvasIdx));
      int lastBin         = std::min(ratioNormBins, nh->GetNbinsX());
      double headIntegral = nh->Integral(1, lastBin);
      if (headIntegral > 0) {
        nh->Scale(1.0 / headIntegral);
      }
      ratioNormHists.push_back(nh);
    }

    if (!ratioNormHists[0]) {
      // No reference histogram on this pad
      continue;
    }

    TLegend *legR = new TLegend(0.7, 0.7, 0.95, 0.95);
    legR->SetBorderSize(0);

    bool firstRatio = true;
    double xmin = 0, xmax = 0;

    for (size_t canvasIdx = 1; canvasIdx < ratioNormHists.size(); canvasIdx++) {
      if (!ratioNormHists[canvasIdx])
        continue;

      TH1F *ratio = (TH1F *)ratioNormHists[canvasIdx]->Clone(Form("%s_ratio%zu_pad%d", ratioTitle, canvasIdx, padIdx));
      ratio->SetTitle(Form("%s;%s;Ratio", ratioNormHists[0]->GetTitle(), ratioNormHists[0]->GetXaxis()->GetTitle()));
      ratio->Divide(ratioNormHists[0]);
      ratio->SetLineColor(colorFor(canvasIdx));
      ratio->SetLineWidth(2);
      ratio->SetStats(0);
      ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
      ratio->GetYaxis()->SetTitle("Ratio");

      if (firstRatio) {
        ratio->Draw("HIST");
        xmin       = ratio->GetXaxis()->GetXmin();
        xmax       = ratio->GetXaxis()->GetXmax();
        firstRatio = false;
      } else {
        ratio->Draw("HIST SAME");
      }

      legR->AddEntry(ratio, Form("%s / %s", entries[canvasIdx].label.c_str(), entries[0].label.c_str()), "l");
    }

    // Reference line at y = 1
    TLine *one = new TLine(xmin, 1.0, xmax, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 2);
    one->Draw("SAME");

    legR->Draw();
    gPad->SetGrid();
  }

  return {compareCanvas, ratioCanvas};
}

void compare_tof_plots(const char *inputFile = "/w/hallc-scshelf2102/c-lad/ehingerl/software/lad_replay/CALIBRATION/"
                                               "lad_tof/files/root/updated_timing_plots_C3_22745-23590_P.root") {

  // Open the input ROOT file
  TFile *f = new TFile(inputFile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Could not open file " << inputFile << std::endl;
    return;
  }

  // Define the three sets of plots. Each set can have any number of entries;
  // ratios are taken with respect to the first entry.
  std::vector<PlotEntry> entries_AllHits = {
      {"KIN/ToF_Vertex_RFcorr_photonCorr/All_Hits/c_ToF_Vertex_RFcorr_photonCorr_All_Hits_all_planes", "All Hits"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_edepCut/All_Hits/c_ToF_Vertex_RFcorr_photonCorr_edepCut_All_Hits_all_planes",
       "edep > const"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_antiEdepCut/All_Hits/"
       "c_ToF_Vertex_RFcorr_photonCorr_antiEdepCut_All_Hits_all_planes",
       "edep < const"}};

  std::vector<PlotEntry> entries_MatchingHit = {
      {"KIN/ToF_Vertex_RFcorr_photonCorr/Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_Matching_Hit_Tol_0_all_planes",
       "All Hits"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_edepCut/Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_edepCut_Matching_Hit_Tol_0_all_planes",
       "edep > const"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_antiEdepCut/Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_antiEdepCut_Matching_Hit_Tol_0_all_planes",
       "edep < const"}};

  std::vector<PlotEntry> entries_AntiMatchingHit = {
      {"KIN/ToF_Vertex_RFcorr_photonCorr/Anti-Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_Anti-Matching_Hit_Tol_0_all_planes",
       "All Hits"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_edepCut/Anti-Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_edepCut_Anti-Matching_Hit_Tol_0_all_planes",
       "edep > const"},
      {"KIN/ToF_Vertex_RFcorr_photonCorr_antiEdepCut/Anti-Matching_Hit_Tol_0/"
       "c_ToF_Vertex_RFcorr_photonCorr_antiEdepCut_Anti-Matching_Hit_Tol_0_all_planes",
       "edep < const"}};

  // Create comparison + ratio canvases for each set
  auto pair_AllHits = createComparisonCanvases(f, entries_AllHits, "comparison_AllHits", "ratio_AllHits");
  if (!pair_AllHits.first) {
    std::cerr << "Failed to create AllHits comparison canvas" << std::endl;
    return;
  }

  auto pair_MatchingHit =
      createComparisonCanvases(f, entries_MatchingHit, "comparison_MatchingHit", "ratio_MatchingHit");
  if (!pair_MatchingHit.first) {
    std::cerr << "Failed to create MatchingHit comparison canvas" << std::endl;
    return;
  }

  auto pair_AntiMatchingHit =
      createComparisonCanvases(f, entries_AntiMatchingHit, "comparison_AntiMatchingHit", "ratio_AntiMatchingHit");
  if (!pair_AntiMatchingHit.first) {
    std::cerr << "Failed to create AntiMatchingHit comparison canvas" << std::endl;
    return;
  }

  // Save canvases (comparison + ratio) to a ROOT file, with one subdirectory per set
  TFile *outputFile = new TFile("tof_comparison.root", "RECREATE");

  auto writeSet = [&](const char *dirName, const std::pair<TCanvas *, TCanvas *> &canvases) {
    TDirectory *dir = outputFile->mkdir(dirName);
    if (!dir)
      return;
    dir->cd();
    if (canvases.first)
      canvases.first->Write();
    if (canvases.second)
      canvases.second->Write();
    outputFile->cd();
  };

  writeSet("AllHits", pair_AllHits);
  writeSet("MatchingHit", pair_MatchingHit);
  writeSet("AntiMatchingHit", pair_AntiMatchingHit);

  outputFile->Close();

  std::cout << "Comparison + ratio canvases saved to tof_comparison.root" << std::endl;

  f->Close();
}