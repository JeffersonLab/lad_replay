#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TTree.h>

struct bcm_calib_param {
  double A;
  double B;
};

bcm_calib_param param_bcm1  = {2489.35, 0.0001753829};
bcm_calib_param param_bcm2  = {250361.6, 0.0001811995};
bcm_calib_param param_unser = {327953, 0.000250};
bcm_calib_param param_bcm4a = {4445.02, 0.00001222493};
bcm_calib_param param_bcm4b = {1782.0, 0.0000146627};
bcm_calib_param param_bcm4c = {2358.5, 0.000010040};

double frq2I(double frq, const bcm_calib_param &param) {
  // Convert frequency to current using calibration parameters
  return (frq - param.A) * param.B;
}

double I2frq(double I, const bcm_calib_param &param) {
  // Convert current to frequency using calibration parameters
  return I / param.B + param.A;
}

void epics2bcm() {
  // Open the input ROOT file (replace with your actual filename)
  TFile *infile = TFile::Open("/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/SCALAR/"
                              "SCALAR_LAD_COIN_22337_0_0_-1.root",
                              "READ");
  if (!infile || infile->IsZombie()) {
    printf("Error: Cannot open input.root\n");
    return;
  }

  // Get the trees
  TTree *E   = (TTree *)infile->Get("E");
  TTree *TSP = (TTree *)infile->Get("TSP");

  // Set up variables to read from tree E
  Double_t ibcm3H04A, ibcm3H04B, ibcm3H04C = 0;
  E->SetBranchAddress("ibcm3H04B", &ibcm3H04B);

  Long64_t nentries = E->GetEntries();
  std::vector<double> event_E, ibcm3H04B_vec, ibcm3H04B_vec_corrected;
  event_E.reserve(nentries);
  ibcm3H04B_vec.reserve(nentries);
  ibcm3H04B_vec_corrected.reserve(nentries);

  // Loop over entries and fill vectors
  for (Long64_t i = 0; i < nentries; ++i) {
    E->GetEntry(i);
    event_E.push_back(i);
    ibcm3H04B_vec.push_back(ibcm3H04B);
    // Apply calibration to ibcm3H04B
    double I = I2frq(ibcm3H04B, param_bcm4b); // Use the appropriate calibration parameters
    ibcm3H04B_vec_corrected.push_back(I);
  }

  Double_t bcm4b_scalerRate = 0;
  TSP->SetBranchAddress("P.BCM4B.scalerRate", &bcm4b_scalerRate);

  Long64_t nentries_TSP = TSP->GetEntries();
  std::vector<double> event_TSP, bcm4b_scalerRate_vec;
  event_TSP.reserve(nentries_TSP);
  bcm4b_scalerRate_vec.reserve(nentries_TSP);

  for (Long64_t i = 0; i < nentries_TSP; ++i) {
    TSP->GetEntry(i);
    event_TSP.push_back(i);
    bcm4b_scalerRate_vec.push_back(bcm4b_scalerRate);
  }

  // Save to output ROOT file
  TFile *outfile = TFile::Open("output.root", "RECREATE");

  TGraph *gr_bcm4b = new TGraph(nentries_TSP, &event_TSP[0], &bcm4b_scalerRate_vec[0]);
  gr_bcm4b->SetTitle("P.BCM4B.scalerRate vs Event Number;Event Number;P.BCM4B.scalerRate");
  gr_bcm4b->Write("P_BCM4B_scalerRate_vs_event");

  // Create a graph
  TGraph *gr = new TGraph(nentries, &event_E[0], &ibcm3H04B_vec[0]);
  gr->SetTitle("ibcm3H04B vs Event Number;Event Number;ibcm3H04B");

  // Apply calibration to ibcm3H04B
  TGraph *gr_corrected = new TGraph(nentries, &event_E[0], &ibcm3H04B_vec_corrected[0]);
  gr_corrected->SetTitle("ibcm3H04B (corrected) vs Event Number;Event Number;ibcm3H04B (corrected)");
  gr_corrected->Write("ibcm3H04B_corrected_vs_event");

  gr->Write("ibcm3H04B_vs_event");
  outfile->Close();

  // Clean up
  infile->Close();
}