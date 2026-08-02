#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <vector>

// Function to calculate FWHM from a histogram
struct FWHMResults {
  double fwhm;
  double sigma;
  double mean;
  bool success;
};

FWHMResults calculateFWHM(TH1D *hist) {
  FWHMResults result = {0.0, 0.0, 0.0, false};

  if (!hist || hist->GetEntries() == 0) {
    return result;
  }

  // Fit with Gaussian
  TF1 *gaus = new TF1("gaus", "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  hist->Fit(gaus, "Q0");

  if (gaus->GetParameter(0) > 0) {
    result.mean    = gaus->GetParameter(1);
    result.sigma   = gaus->GetParameter(2);
    result.fwhm    = 2.355 * result.sigma; // FWHM = 2.355 * sigma for Gaussian
    result.success = true;
  }

  delete gaus;
  return result;
}

void laser_resolution() {
  // Set batch mode to prevent X11 forwarding
  gROOT->SetBatch(kTRUE);

  // Open the ROOT file and get the tree
  TFile *file =
      new TFile("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/COSMICS/laser/LAD_wGEM_cosmic_hall_23798_-1.root");
  TTree *T = (TTree *)file->Get("T");

  if (!T) {
    std::cerr << "Error: Could not find tree T in file" << std::endl;
    return;
  }

  // Define plane IDs and number of planes
  std::vector<std::string> planes = {"000", "001", "100", "101", "200"};
  const int NPLANES               = 11;

  // Create canvases for HodoHitTime and time difference
  std::vector<TCanvas *> canvases_time;
  std::vector<TCanvas *> canvases_diff;
  std::vector<std::vector<TH1D *>> histograms_time;
  std::vector<std::vector<TH1D *>> histograms_diff;

  // Initialize canvases and histograms
  for (size_t p = 0; p < planes.size(); p++) {
    TCanvas *c_time =
        new TCanvas(Form("c_time_%s", planes[p].c_str()), Form("HodoHitTime - Plane %s", planes[p].c_str()), 1200, 800);
    TCanvas *c_diff = new TCanvas(Form("c_diff_%s", planes[p].c_str()),
                                  Form("HodoHitTime - Plane 000 Paddle 0 - Plane %s", planes[p].c_str()), 1200, 800);

    c_time->Divide(4, 3);
    c_diff->Divide(4, 3);

    canvases_time.push_back(c_time);
    canvases_diff.push_back(c_diff);

    std::vector<TH1D *> hists_time;
    std::vector<TH1D *> hists_diff;

    for (int i = 0; i < NPLANES; i++) {
      TH1D *h_time = new TH1D(Form("h_time_%s_plane%d", planes[p].c_str(), i), Form("Paddle %d", i), 200, -20, 20);
      TH1D *h_diff = new TH1D(Form("h_diff_%s_plane%d", planes[p].c_str(), i), Form("Paddle %d", i), 200, -1, 1);

      h_time->SetXTitle("HodoHitTime (ns)");
      h_time->SetYTitle("Counts");
      h_diff->SetXTitle("HodoHitTime - Plane 000 Paddle 0 (ns)");
      h_diff->SetYTitle("Counts");

      hists_time.push_back(h_time);
      hists_diff.push_back(h_diff);
    }

    histograms_time.push_back(hists_time);
    histograms_diff.push_back(hists_diff);
  }

  // Set branch addresses
  // Arrays for each plane
  std::vector<double *> hodohittime_arrays;
  std::vector<double *> paddlenum_arrays;
  std::vector<int *> ndata_arrays;

  for (const auto &plane : planes) {
    double *hittime = new double[100];
    double *paddle  = new double[100];
    int *ndata      = new int;

    T->SetBranchAddress(Form("L.ladhod.%s.HodoHitTime", plane.c_str()), hittime);
    T->SetBranchAddress(Form("L.ladhod.%s.HodoHitPaddleNum", plane.c_str()), paddle);
    T->SetBranchAddress(Form("Ndata.L.ladhod.%s.HodoHitTime", plane.c_str()), ndata);

    hodohittime_arrays.push_back(hittime);
    paddlenum_arrays.push_back(paddle);
    ndata_arrays.push_back(ndata);
  }

  // Loop through all entries
  Long64_t nentries = T->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    // Display progress
    if (i % 1000 == 0 && i > 0) {
      int progress = (int)((double)i / nentries * 100.0);
      std::cout << "\rProgress: " << progress << "%" << std::flush;
    }

    T->GetEntry(i);

    // Get reference time from plane 000, paddle 0
    double reference_time = 0.0;
    if (*ndata_arrays[0] > 0) {
      reference_time = hodohittime_arrays[0][0];
    }

    // Loop through each plane
    for (size_t p = 0; p < planes.size(); p++) {
      int n = *ndata_arrays[p];

      // Loop through all hits in this plane for this event
      for (int j = 0; j < n && j < 100; j++) {
        int paddle     = paddlenum_arrays[p][j] - 1;
        double hittime = hodohittime_arrays[p][j];

        // Make sure paddle number is within range
        if (paddle >= 0 && paddle < NPLANES) {
          histograms_time[p][paddle]->Fill(hittime);
          histograms_diff[p][paddle]->Fill(hittime - reference_time);
        }
      }
    }
  }

  // Draw histograms
  for (size_t p = 0; p < planes.size(); p++) {
    canvases_time[p]->cd();
    for (int i = 0; i < NPLANES; i++) {
      canvases_time[p]->cd(i + 1);
      histograms_time[p][i]->Draw();

      // Calculate and display FWHM and sigma
      FWHMResults results = calculateFWHM(histograms_time[p][i]);
      if (results.success) {
        TPaveText *pt = new TPaveText(0.6, 0.7, 0.95, 0.95, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.08);
        pt->AddText(Form("FWHM: %.3f ns", results.fwhm));
        pt->AddText(Form("#sigma: %.3f ns", results.sigma));
        pt->Draw("same");
      }
    }

    canvases_diff[p]->cd();
    for (int i = 0; i < NPLANES; i++) {
      canvases_diff[p]->cd(i + 1);
      histograms_diff[p][i]->Draw();

      // Calculate and display FWHM and sigma
      FWHMResults results = calculateFWHM(histograms_diff[p][i]);
      if (results.success) {
        TPaveText *pt = new TPaveText(0.6, 0.7, 0.95, 0.95, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.08);
        pt->AddText(Form("FWHM: %.3f ns", results.fwhm));
        pt->AddText(Form("#sigma: %.3f ns", results.sigma));
        pt->Draw("same");
      }
    }
  }

  // Write canvases to output file
  TFile *outputFile = new TFile("laser_resolution_output.root", "RECREATE");
  for (size_t p = 0; p < planes.size(); p++) {
    canvases_time[p]->Write();
    canvases_diff[p]->Write();
  }
  outputFile->Close();

  // Clean up
  for (size_t p = 0; p < planes.size(); p++) {
    delete hodohittime_arrays[p];
    delete paddlenum_arrays[p];
    delete ndata_arrays[p];
  }

  file->Close();
}
