#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>
#include <TF1.h>
#include <TPaveText.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

// Helper function to fit two Gaussians and calculate average FWHM
void FitTwoGaussians(TH1D *hist, double &fwhm_avg, double &sigma) {
  fwhm_avg = 0.0;
  sigma = 0.0;

  // Create a fit function: sum of two Gaussians
  // [0] = amplitude of first Gaussian
  // [1] = center of first Gaussian
  // [2] = sigma of first Gaussian
  // [3] = amplitude of second Gaussian
  // [4] = center of second Gaussian
  // [5] = sigma of second Gaussian
  TF1 *fit = new TF1("twoGaus", 
    "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2)",
    hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  // Set initial parameters
  double max_bin_content = hist->GetMaximum();
  fit->SetParameter(0, max_bin_content);      // Amplitude 1
  fit->SetParameter(1, -1.0);                  // Center 1 (start at -1)
  fit->SetParameter(2, 0.3);                   // Sigma 1
  fit->SetParameter(3, max_bin_content);      // Amplitude 2
  fit->SetParameter(4, 1.0);                   // Center 2 (start at +1)
  fit->SetParameter(5, 0.3);                   // Sigma 2

  // Perform the fit
  hist->Fit(fit, "Q0");  // Q = quiet, 0 = don't draw

  // Extract sigma values from fit
  double sigma1 = fit->GetParameter(2);
  double sigma2 = fit->GetParameter(5);

  // Calculate average FWHM from the two sigmas
  // FWHM = 2.355 * sigma
  double fwhm1 = 2.355 * sigma1;
  double fwhm2 = 2.355 * sigma2;
  fwhm_avg = (fwhm1 + fwhm2) / 2.0;
  sigma = (sigma1 + sigma2) / 2.0;

  delete fit;
}

void cosmic_resolution() {
  // Set batch mode to prevent X11 forwarding
  gROOT->SetBatch(kTRUE);

  // Adjustable parameter for minimum adjacent hits
  const int MIN_ADJACENT_HITS = 5;

  // Open the ROOT file and get the tree
  TFile *file =
      new TFile("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_22485_-1.root");
  TTree *T = (TTree *)file->Get("T");

  if (!T) {
    std::cerr << "Error: Could not find tree T in file" << std::endl;
    return;
  }

  // Define plane IDs and number of planes
  std::vector<std::string> planes = {"000", "001", "100", "101", "200"};
  const int NPLANES               = 11;
  const int NPLANE_IDS            = planes.size();

  // Structure to hold one hit per plane + paddle
  struct HitData {
    double time;
    int paddle;
    bool has_hit;
  };

  // Create canvases and histograms for time differences between adjacent bars
  // For each plane, we'll have 10 histograms (NPLANES - 1)
  std::vector<TCanvas *> canvases;
  std::vector<std::vector<TH1D *>> histograms_diff;

  for (size_t p = 0; p < planes.size(); p++) {
    // Create canvas for this plane
    TCanvas *c = new TCanvas(Form("c_plane_%s", planes[p].c_str()), 
                             Form("Plane %s: Time Differences Between Adjacent Bars", planes[p].c_str()), 
                             1400, 900);
    c->Divide(5, 2);  // 5 columns x 2 rows for 10 histograms
    canvases.push_back(c);

    std::vector<TH1D *> hists_diff;

    for (int i = 0; i < NPLANES - 1; i++) {
      TH1D *h_diff = new TH1D(Form("h_diff_%s_bars%d_%d", planes[p].c_str(), i, i + 1),
                              Form("Plane %s: Bar %d - Bar %d", planes[p].c_str(), i, i + 1), 200, -3, 3);

      h_diff->SetXTitle("Time Difference (ns)");
      h_diff->SetYTitle("Counts");

      hists_diff.push_back(h_diff);
    }

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

    // Create hit data array for this event: one hit per plane + paddle
    // hits[plane_id][paddle_num]
    HitData hits[NPLANE_IDS][NPLANES];

    // Initialize all hits as empty
    for (int p = 0; p < NPLANE_IDS; p++) {
      for (int pad = 0; pad < NPLANES; pad++) {
        hits[p][pad].has_hit = false;
        hits[p][pad].time    = 0.0;
        hits[p][pad].paddle  = -1;
      }
    }

    // Fill the hit array with one hit per plane + paddle
    for (size_t p = 0; p < planes.size(); p++) {
      int n = *ndata_arrays[p];

      // For each plane, find the first hit in each paddle
      for (int j = 0; j < n && j < 100; j++) {
        int paddle     = paddlenum_arrays[p][j] - 1; // Convert to 0-indexed
        double hittime = hodohittime_arrays[p][j];

        // Only store if paddle is valid and we haven't already stored a hit for this paddle
        if (paddle >= 0 && paddle < NPLANES && !hits[p][paddle].has_hit) {
          hits[p][paddle].has_hit = true;
          hits[p][paddle].time    = hittime;
          hits[p][paddle].paddle  = paddle;
        }
      }
    }

    // Check for >=5 adjacent hits and fill histograms
    for (size_t p = 0; p < planes.size(); p++) {
      // Look for sequences of adjacent hits with >=MIN_ADJACENT_HITS hits
      for (int start_pad = 0; start_pad < NPLANES; start_pad++) {
        // Count adjacent hits starting from start_pad
        int adjacent_count = 0;
        for (int pad = start_pad; pad < NPLANES; pad++) {
          if (hits[p][pad].has_hit) {
            adjacent_count++;
          } else {
            break; // Stop at first gap
          }
        }

        // If we have >=MIN_ADJACENT_HITS adjacent hits, fill time difference histograms
        if (adjacent_count >= MIN_ADJACENT_HITS) {
          for (int j = 0; j < adjacent_count - 1; j++) {
            int bar1     = start_pad + j;
            int bar2     = start_pad + j + 1;
            double time1 = hits[p][bar1].time;
            double time2 = hits[p][bar2].time;
            double diff  = time2 - time1;

            // Fill the histogram for this bar pair
            // The histogram index is (bar1), since we have histograms for 0-1, 1-2, ..., 9-10
            if (bar1 < NPLANES - 1) {
              histograms_diff[p][bar1]->Fill(diff);
            }
          }
        }
      }
    }
  }

  // Draw histograms on canvases and fit with two Gaussians
  for (size_t p = 0; p < planes.size(); p++) {
    for (int i = 0; i < NPLANES - 1; i++) {
      canvases[p]->cd(i + 1);
      histograms_diff[p][i]->Draw();

      // Fit with two Gaussians
      double fwhm_avg = 0.0;
      double sigma = 0.0;
      FitTwoGaussians(histograms_diff[p][i], fwhm_avg, sigma);

      // Add text box with FWHM and sigma
      TPaveText *pt = new TPaveText(0.6, 0.7, 0.95, 0.9, "NDC");
      pt->SetBorderSize(1);
      pt->SetFillColor(0);
      pt->AddText(Form("FWHM = %.3f ns", fwhm_avg));
      pt->AddText(Form("#sigma = %.3f ns", sigma));
      pt->Draw();
    }
  }

  // Write canvases and histograms to output file
  TFile *outputFile = new TFile("cosmic_resolution_output.root", "RECREATE");
  for (size_t p = 0; p < planes.size(); p++) {
    canvases[p]->Write();
    // for (int i = 0; i < NPLANES - 1; i++) {
    //   histograms_diff[p][i]->Write();
    // }
  }
  outputFile->Close();

  std::cout << "\nAnalysis complete. Output written to cosmic_resolution_output.root" << std::endl;

  // Clean up
  for (size_t p = 0; p < planes.size(); p++) {
    delete hodohittime_arrays[p];
    delete paddlenum_arrays[p];
    delete ndata_arrays[p];
  }

  file->Close();
}
