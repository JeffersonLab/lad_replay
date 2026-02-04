/*HMS Hodo Calibration script to determine scin. prop velocity, cable time offsets per paddle,
  and time offsets between paddles in different planes
  Author: Carlos Yero
  Dated: June 5, 2018
*/

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
#include <TDecompSVD.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TPolyLine.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

const static double TdcTimeTWCorr_MAX = 100000.0; // Was originally 100.
const static double AdcTimeTWCorr_MAX = 100000.0;
const static int minADC_PulseAmp      = 40; // Minimum ADC Pulse Amplitude to consider for Time Walk correction

// void fitHodoCalib(TString filename, Int_t runNUM, Bool_t cosmic_flag = kFALSE) {
void fitHodoCalib(Int_t runNUM) {

  // Flag to use average velocity for all paddles instead of individual values
  Bool_t use_vp_avg = kTRUE;

  gROOT->SetBatch(kTRUE);
  // TString filename = Form("../../ROOTfiles/COSMICS/LAD_wREF_cosmic_hall_%d_-1.root", runNUM);
  // TString filename = Form("../../ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_%d_-1.root", runNUM);

  gStyle->SetOptFit();
  // gROOT->SetBatch(kTRUE); // do not display plots

  Int_t evtNUM = 30000;

  // TFile *data_file = new TFile(filename, "READ");
  // TTree *T         = (TTree *)data_file->Get("T");

  // Create a TChain named "T" and add the listed files
  TChain *T = new TChain("T");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23063_6_7_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23038_9_10_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23064_0_0_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23065_0_0_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23041_0_1_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23071_0_1_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23046_3_4_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23073_9_9_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23047_9_9_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23053_6_6_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23076_6_7_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23054_0_0_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23077_0_1_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23055_9_10_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23095_3_3_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23096_9_9_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23057_9_10_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23097_3_3_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23058_9_9_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23098_9_9_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23061_6_7_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23062_6_6_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23100_9_10_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23038_6_8_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23038_0_2_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23038_3_5_-1.root");

  T->Add("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_22483_-1.root");
  T->Add("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_22487_-1.root");
  // Report what was added
  std::cout << "TChain T created: " << T->GetNtrees() << " files, total entries = " << T->GetEntries() << std::endl;

  /******Define Fixed Quantities********/
  static const Int_t NPLANES                     = 6;
  static const Int_t MAX_PADDLES                 = 11;
  static const Int_t SIDES                       = 2;
  TString spec                                   = "L"; // Usually "L"
  TString det                                    = "ladhod";
  TString pl_names[NPLANES]                      = {"000", "001", "100", "101", "200", "REFBAR"};
  TString side_names[2]                          = {"GoodTop", "GoodBtm"};
  TString nsign[2]                               = {"+", "-"};
  Int_t maxPMT[NPLANES]                          = {11, 11, 11, 11, 11, 1};
  Int_t refBar                                   = 0;
  Int_t refPlane                                 = 0;
  Double_t lladhodo_velArr[NPLANES][MAX_PADDLES] = {0.0}; // store hhodo velocity parameters (1/slope of the line fit)
  Double_t lladhodo_velArr_err[NPLANES][MAX_PADDLES][2] = {
      0.0};                                                 // store hhodo velocity parameters (1/slope of the line fit)
  Double_t lladhodo_cableArr[NPLANES][MAX_PADDLES] = {0.0}; // store hhodo cableLength differences (y-int of line fit)
  Double_t lladhodo_LCoeff[NPLANES][MAX_PADDLES]   = {0.0}; // Variables to write out LCoeff. parameter file
  Double_t lladhodo_sigArr[NPLANES][MAX_PADDLES]   = {0.0}; // store hhodo sigma parameters
  Double_t lladhodo_sigArr_LCoeff[NPLANES][MAX_PADDLES] = {0.0}; // store hhodo sigma parameters for LCoeff.
  Double_t lladhodo_velArr_FADC[NPLANES][MAX_PADDLES]   = {
      0.0}; // store hhodo velocity parameters (1/slope of the line fit) for FADC
  Double_t lladhodo_cableArr_FADC[NPLANES][MAX_PADDLES] = {
      0.0}; // store hhodo cableLength differences (y-int of line fit) for FADC
  Double_t lladhodo_LCoeff_FADC[NPLANES][MAX_PADDLES] = {0.0}; // Variables to write out LCoeff. parameter file for FADC
  Double_t lladhodo_sigArr_FADC[NPLANES][MAX_PADDLES] = {0.0}; // store hhodo sigma parameters for FADC
  Double_t lladhodo_sigArr_LCoeff_FADC[NPLANES][MAX_PADDLES] = {0.0};        // store hhodo sigma parameters for FADC
  Double_t vp                                                = 30.0;         // speed of light [cm/ns]
  Double_t paddle_vel_fit_max_threshold                      = 0.10;         // 0-1 range, 0.05 = 5% threshold
  Double_t paddle_vel_fit_max_threshold_err[2]               = {0.05, 0.15}; // 0-1 range, 0.05 = 5% threshold
  Double_t barlength[NPLANES][MAX_PADDLES]                   = {
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1}};

  // Plotting parameters
  static const Int_t TDC_T_NBINS         = 1000;
  static const Double_t TDC_T_MIN        = -60.0; // 60
  static const Double_t TDC_T_MAX        = 60.0;  //-60
  static const Int_t ADC_TDC_diff_NBINS  = 100;
  static const Double_t ADC_TDC_diff_MIN = -1000.0; // 360.0
  static const Double_t ADC_TDC_diff_MAX = -330.0;  //-330.0
  static const Int_t DiffTime_NBINS      = 500;
  static const Double_t DiffTime_MIN     = -40.0; // -40.0
  static const Double_t DiffTime_MAX     = 40.0;  // 40.0
  static const Int_t FADC_T_NBINS        = 1000;
  static const Double_t FADC_T_MIN       = -60.0;
  static const Double_t FADC_T_MAX       = 60.0;
  /******Define Leafs to be read from TTree******/

  //---Names---
  TString base;
  TString nTdcTimeUnCorr;
  TString nTdcTimeTWCorr;
  TString nAdcPulseTime;
  TString nAdcPulseAmp;

  TString nTrackXPos;
  TString nTrackYPos;
  TString nDiffTWCorr;

  //---Variables---
  Double_t TdcTimeUnCorr[NPLANES][SIDES][MAX_PADDLES];
  Double_t TdcTimeTWCorr[NPLANES][SIDES][MAX_PADDLES];
  Double_t AdcPulseTime[NPLANES][SIDES][MAX_PADDLES];
  Double_t AdcPulseAmp[NPLANES][SIDES][MAX_PADDLES];

  Double_t hcal_etrkNorm;
  Double_t hcer_npeSum;
  Double_t hdc_ntrack;
  Double_t hod_nhits[NPLANES];
  Double_t beta;

  /******Define Matrices/Vectors and related *******/
  // Int_t good_pad[PLANES]; // keep track of good paddle hits
  // Int_t ngood = 0;

  // static const Int_t npar = 55 - 1; // reference paddle ??? is fixed (so we actually have 54)
  // static const Int_t ROW  = 6 * evtNUM;
  // static const Int_t COL  = npar;

  // Int_t row1, row2, row3, row4, row5, row6; // keep track of ith row element
  // Int_t cnt;                                // keep track of good plane hits

  // /*******Define Canvas and Histograms*******/

  // //----Canvas----
  TCanvas *TWAvg_canv[NPLANES];
  TCanvas *TWAvg_canv_FADC[NPLANES];
  TCanvas *TWAvg_canv_2D[NPLANES];
  TCanvas *TWDiff_canv[NPLANES];
  TCanvas *TWDiff_canv_FADC[NPLANES];

  TCanvas *TWUnCorr_canv[NPLANES][SIDES];
  TCanvas *TWCorr_canv[NPLANES][SIDES];
  TCanvas *TWCorr_canv_FADC[NPLANES][SIDES];

  TCanvas *TWCorr_BarLCoef[NPLANES];
  TCanvas *TWCorr_BarLCoef_FADC[NPLANES];
  // //----Histograms----
  TH1F *h1Hist_TWAvg[NPLANES][MAX_PADDLES];     // (TWCorr_Top + TWCorr_Btm) / 2       <-------
  TH1F *h1Hist_TWAvg_CUT[NPLANES][MAX_PADDLES]; //<------
  TH1F *h1Hist_TWDiff[NPLANES][MAX_PADDLES];    // (TWCorr_Top - TWCorr_Btm) / 2       <-------

  TH1F *h1Hist_TWAvg_FADC[NPLANES][MAX_PADDLES];     // (TWCorr_Top + TWCorr_Btm) / 2 for FADC      <-------
  TH1F *h1Hist_TWAvg_CUT_FADC[NPLANES][MAX_PADDLES]; //<------
  TH1F *h1Hist_TWDiff_FADC[NPLANES][MAX_PADDLES];    // (TWCorr_Top - TWCorr_Btm) / 2 for FADC      <-------

  TH2F *h2Hist_TW_UnCorr[NPLANES][SIDES][MAX_PADDLES]; // Time-Walk Uncorrected vs. ADC Pulse Amp Hist
  TH2F *h2Hist_TW_Corr[NPLANES][SIDES][MAX_PADDLES];   // Time-Walk Corrected vs. ADC Pulse Amp Hist

  TH1F *h1Hist_TWCorr_BarLCoef[NPLANES][MAX_PADDLES];
  TH1F *h1Hist_TWCorr_BarLCoef_FADC[NPLANES][MAX_PADDLES];
  // /*******Define Fit Functions and Related Variables*******/

  Double_t Mean, mean_FADC;     // variable to get Mean to make a 3sig cut
  Double_t StdDev, StdDev_FADC; // variable to get satndar deviation to make a 3sig cut
  Double_t nSig, nSig_FADC;     // multiple of Sigma used for sigmaCut

  /********Initialize HISTOS and GET TTREE VARIABLES*********/

  // Loop over hodo planes
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Loop over hodo side
    for (Int_t side = 0; side < SIDES; side++) {

      // Loop over hodo PMTs
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        // Initialize Histograms
        h2Hist_TW_UnCorr[npl][side][ipmt] =
            new TH2F(Form("TW_UnCorr PMT %s%d%s", pl_names[npl].Data(), ipmt + 1, nsign[side].Data()),
                     Form("PMT %s%d%s: UnCorr. (TDC - ADC) Pulse Time vs. ADC Pulse Amplitude ", pl_names[npl].Data(),
                          ipmt + 1, nsign[side].Data()),
                     600, 0, 420, ADC_TDC_diff_NBINS, ADC_TDC_diff_MIN, ADC_TDC_diff_MAX);
        h2Hist_TW_Corr[npl][side][ipmt] =
            new TH2F(Form("TW_Corr PMT %s%d%s", pl_names[npl].Data(), ipmt + 1, nsign[side].Data()),
                     Form("PMT %s%d%s: Corr. (TDC - ADC) Pulse Time vs. ADC Pulse Amplitude ", pl_names[npl].Data(),
                          ipmt + 1, nsign[side].Data()),
                     600, 0, 420, ADC_TDC_diff_NBINS, ADC_TDC_diff_MIN, ADC_TDC_diff_MAX);

        h2Hist_TW_UnCorr[npl][side][ipmt]->GetYaxis()->SetTitle("Time Walk UnCorr.(TDC - ADC) Pulse Time (ns)");
        h2Hist_TW_UnCorr[npl][side][ipmt]->GetXaxis()->SetTitle("ADC Pulse Amplitude (mV)");
        h2Hist_TW_UnCorr[npl][side][ipmt]->GetXaxis()->CenterTitle();

        h2Hist_TW_Corr[npl][side][ipmt]->GetYaxis()->SetTitle("Time Walk Corr.(TDC - ADC) Pulse Time (ns)");

        h2Hist_TW_Corr[npl][side][ipmt]->GetXaxis()->SetTitle("ADC Pulse Amplitude (mV)");
        h2Hist_TW_Corr[npl][side][ipmt]->GetXaxis()->CenterTitle();

        //----Define TTree Leaf Names-----
        base = spec + "." + det + "." + pl_names[npl];

        nTdcTimeUnCorr = base + "." + side_names[side] + "TdcTimeUnCorr";
        nTdcTimeTWCorr = base + "." + side_names[side] + "TdcTimeWalkCorr";
        nAdcPulseTime  = base + "." + side_names[side] + "AdcPulseTime";
        nAdcPulseAmp   = base + "." + side_names[side] + "AdcPulseAmp";

        //------Set Branch Address-------
        T->SetBranchAddress(nTdcTimeUnCorr, &TdcTimeUnCorr[npl][side]);
        T->SetBranchAddress(nTdcTimeTWCorr, &TdcTimeTWCorr[npl][side]);
        T->SetBranchAddress(nAdcPulseTime, &AdcPulseTime[npl][side]);
        T->SetBranchAddress(nAdcPulseAmp, &AdcPulseAmp[npl][side]);

      } // end loop over hodo PMTs

    } // end loop over hodo side

  } // end loop over hodo planes

  // Get the average value of the branch with the name nTdcTimeTWCorr (first 100 entries only)
  double tdc_avg           = 0.0;
  double sum_TdcTimeTWCorr = 0.0;
  int count_TdcTimeTWCorr  = 0;
  double sum_AdcPulseTime  = 0.0;
  int count_AdcPulseTime   = 0;
  Long64_t nAvgEntries     = std::min((Long64_t)100, T->GetEntries());
  for (Long64_t i = 0; i < nAvgEntries; i++) {
    T->GetEntry(i);
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      for (Int_t side = 0; side < SIDES; side++) {
        for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
          if (TdcTimeTWCorr[npl][side][ipmt] < TdcTimeTWCorr_MAX) {
            sum_TdcTimeTWCorr += TdcTimeTWCorr[npl][side][ipmt];
            ++count_TdcTimeTWCorr;
          }
          if (AdcPulseTime[npl][side][ipmt] < AdcTimeTWCorr_MAX) {
            sum_AdcPulseTime += AdcPulseTime[npl][side][ipmt];
            ++count_AdcPulseTime;
          }
        }
      }
    }
  }
  double avg_TdcTimeTWCorr = (count_TdcTimeTWCorr > 0) ? (sum_TdcTimeTWCorr / count_TdcTimeTWCorr) : 0.0;
  double avg_AdcPulseTime  = (count_AdcPulseTime > 0) ? (sum_AdcPulseTime / count_AdcPulseTime) : 0.0;
  cout << "Average value of nTdcTimeTWCorr (first 100 entries): " << avg_TdcTimeTWCorr << endl;
  cout << "Average value of AdcPulseTime (first 100 entries): " << avg_AdcPulseTime << endl;

  // Create 1D histograms for TWAvg, TWAvg_CUT, and TWDiff for each plane and paddle (side == 0)
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      h1Hist_TWAvg[npl][ipmt] =
          new TH1F(Form("Avg. Time: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("Paddle %s%d: Time-Walk Corrected Average Time", pl_names[npl].Data(), ipmt + 1), TDC_T_NBINS,
                   TDC_T_MIN + avg_TdcTimeTWCorr, TDC_T_MAX + avg_TdcTimeTWCorr);

      h1Hist_TWAvg_CUT[npl][ipmt] =
          new TH1F(Form("Avg. Time CUT: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("Paddle %s%d: Time-Walk Corrected Average (CUT)", pl_names[npl].Data(), ipmt + 1), TDC_T_NBINS,
                   TDC_T_MIN + avg_TdcTimeTWCorr, TDC_T_MAX + avg_TdcTimeTWCorr);

      h1Hist_TWDiff[npl][ipmt] =
          new TH1F(Form("Diff. Time: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("Paddle %s%d: Time-Walk Corrected Difference Time (Top - Btm)", pl_names[npl].Data(), ipmt + 1),
                   DiffTime_NBINS, DiffTime_MIN, DiffTime_MAX);

      // Set Axis Titles
      h1Hist_TWAvg[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. TDC Average Paddle Time (ns)");
      h1Hist_TWAvg_CUT[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. TDC Average Paddle Time (ns)");

      h1Hist_TWAvg[npl][ipmt]->GetXaxis()->CenterTitle();
      h1Hist_TWAvg_CUT[npl][ipmt]->GetXaxis()->CenterTitle();

      // Initialize FADC 1D histograms analogous to the TDC ones, using FADC average time instead of TDC avg
      h1Hist_TWAvg_FADC[npl][ipmt] =
          new TH1F(Form("Avg. Time FADC: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("Paddle %s%d: FADC Time-Walk Corrected Average Time", pl_names[npl].Data(), ipmt + 1),
                   FADC_T_NBINS, FADC_T_MIN + avg_AdcPulseTime, FADC_T_MAX + avg_AdcPulseTime);

      h1Hist_TWAvg_CUT_FADC[npl][ipmt] =
          new TH1F(Form("Avg. Time CUT FADC: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("Paddle %s%d: FADC Time-Walk Corrected Average (CUT)", pl_names[npl].Data(), ipmt + 1),
                   FADC_T_NBINS, FADC_T_MIN + avg_AdcPulseTime, FADC_T_MAX + avg_AdcPulseTime);

      h1Hist_TWDiff_FADC[npl][ipmt] = new TH1F(
          Form("Diff. Time FADC: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
          Form("Paddle %s%d: FADC Time-Walk Corrected Difference Time (Top - Btm)", pl_names[npl].Data(), ipmt + 1),
          DiffTime_NBINS, DiffTime_MIN, DiffTime_MAX);

      // Set axis titles and center them
      h1Hist_TWAvg_FADC[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. FADC Average Paddle Time (ns)");
      h1Hist_TWAvg_CUT_FADC[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. FADC Average Paddle Time (ns)");
      h1Hist_TWDiff_FADC[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. FADC Difference Time (Top - Btm)");

      h1Hist_TWAvg_FADC[npl][ipmt]->GetXaxis()->CenterTitle();
      h1Hist_TWAvg_CUT_FADC[npl][ipmt]->GetXaxis()->CenterTitle();
      h1Hist_TWDiff_FADC[npl][ipmt]->GetXaxis()->CenterTitle();
    }
  }

  //**************************************************************//
  // FIRST PASS OF EVENT LOOP (Get the StdDev) of (TDC+ + TDC-)/2 //
  //**************************************************************//
  cout << "Initializing 1st Pass of Event Loop: " << endl;

  Long64_t nentries = T->GetEntries();
  // nentries          = 20000; // Hard coded for testing

  // Loop over all entries
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);

    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over pmt
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        if (AdcPulseAmp[npl][0][ipmt] < minADC_PulseAmp || AdcPulseAmp[npl][1][ipmt] < minADC_PulseAmp) {
          continue; // skip if ADC Pulse Amplitude is below threshold
        }
        if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {
          // Fill Average TW Corr TDC Time
          h1Hist_TWAvg[npl][ipmt]->Fill((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.);

        } // end time cut

        if (AdcPulseTime[npl][0][ipmt] < AdcTimeTWCorr_MAX && AdcPulseTime[npl][1][ipmt] < AdcTimeTWCorr_MAX) {
          // Fill Average TW Corr FADC Time
          h1Hist_TWAvg_FADC[npl][ipmt]->Fill((AdcPulseTime[npl][0][ipmt] + AdcPulseTime[npl][1][ipmt]) / 2.);
        }

      } // end pmt loop

    } // end plane loop

    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";

  } // end loop over entries

  // Set cut on Sigma,
  nSig = nSig_FADC = 1;

  //************************************//
  //    SECOND PASS OF EVENT LOOP       //
  //************************************//

  cout << "Initializing 2nd Pass of Event Loop: " << endl;

  // Loop over all entries
  for (Long64_t i = 0; i < nentries; i++) {

    T->GetEntry(i);

    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over plane side
      for (Int_t side = 0; side < SIDES; side++) {

        // Loop over pmt
        for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

          // Get Standard deviation from initial entry fill
          StdDev = h1Hist_TWAvg[npl][ipmt]->GetStdDev();
          Mean   = h1Hist_TWAvg[npl][ipmt]->GetMean();
          // Get Standard deviation from initial entry fill for FADC
          StdDev_FADC = h1Hist_TWAvg_FADC[npl][ipmt]->GetStdDev();
          mean_FADC   = h1Hist_TWAvg_FADC[npl][ipmt]->GetMean();

          // FIll Uncorrected/Corrected Time Walk Histos
          h2Hist_TW_UnCorr[npl][side][ipmt]->Fill(AdcPulseAmp[npl][side][ipmt],
                                                  TdcTimeUnCorr[npl][side][ipmt] - AdcPulseTime[npl][side][ipmt]);
          h2Hist_TW_Corr[npl][side][ipmt]->Fill(AdcPulseAmp[npl][side][ipmt],
                                                TdcTimeTWCorr[npl][side][ipmt] - AdcPulseTime[npl][side][ipmt]);

          // Add Time Cuts to get rid of kBig - kBig values, which yielded high evt density at zero
          if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {

            if (side == 0) // require only one side, as a time diff between the two ends of a paddle is take
            {
              if ((((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.) > (Mean - nSig * StdDev)) &&
                  (((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.) < (Mean + nSig * StdDev))) {

                h1Hist_TWAvg_CUT[npl][ipmt]->Fill((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.);
                h1Hist_TWDiff[npl][ipmt]->Fill((TdcTimeTWCorr[npl][1][ipmt] - TdcTimeTWCorr[npl][0][ipmt]) / 2.);

              } // end 3SIGMA CUT of TW Corr Time

            } // end require ONLY single side

          } // end time cuts

          // Repeat analogous cuts/fills for FADC
          if (AdcPulseTime[npl][0][ipmt] < AdcTimeTWCorr_MAX && AdcPulseTime[npl][1][ipmt] < AdcTimeTWCorr_MAX) {
            if (side == 0) {
              Double_t avgFADC = (AdcPulseTime[npl][0][ipmt] + AdcPulseTime[npl][1][ipmt]) / 2.0;
              if (avgFADC > (mean_FADC - nSig_FADC * StdDev_FADC) && avgFADC < (mean_FADC + nSig_FADC * StdDev_FADC)) {
                h1Hist_TWAvg_CUT_FADC[npl][ipmt]->Fill(avgFADC);
                h1Hist_TWDiff_FADC[npl][ipmt]->Fill((AdcPulseTime[npl][1][ipmt] - AdcPulseTime[npl][0][ipmt]) / 2.0);
              }
            }
          }

        } // end pmt loop

      } // end side loop

    } // end plane loop

    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";

  } // end entry loop

  /***********DRAW HISTOGRAMS TO CANVAS***************/
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Create Canvas to store TW-Corr Time/Dist vs. trk position
    TWAvg_canv[npl] = new TCanvas(Form("TWAvg_%d", npl), Form("TWAvg, plane %s", pl_names[npl].Data()), 1000, 700);
    TWAvg_canv_FADC[npl] =
        new TCanvas(Form("TWAvgFADC_%d", npl), Form("TWAvg FADC, plane %s", pl_names[npl].Data()), 1000, 700);
    TWAvg_canv_2D[npl] =
        new TCanvas(Form("TWAvg2D_%d", npl), Form("TWAvg2D, plane %s", pl_names[npl].Data()), 1000, 700);
    TWDiff_canv[npl] = new TCanvas(Form("TWDiff_%d", npl), Form("TWDiff, plane %s", pl_names[npl].Data()), 1000, 700);
    TWDiff_canv_FADC[npl] =
        new TCanvas(Form("TWDiffFADC_%d", npl), Form("TWDiff FADC, plane %s", pl_names[npl].Data()), 1000, 700);

    if (npl < 5) {
      TWAvg_canv[npl]->Divide(4, 3);
      TWAvg_canv_FADC[npl]->Divide(4, 3);
      TWDiff_canv[npl]->Divide(4, 3);
      TWDiff_canv_FADC[npl]->Divide(4, 3);
      TWAvg_canv_2D[npl]->Divide(4, 3);
    }

    // Loop over plane side
    for (Int_t side = 0; side < SIDES; side++) {
      // Create Canvas
      TWUnCorr_canv[npl][side] =
          new TCanvas(Form("TWUnCorrCanv%d%d", npl, side),
                      Form("plane %s_%s", pl_names[npl].Data(), side_names[side].Data()), 1000, 700);
      TWCorr_canv[npl][side] =
          new TCanvas(Form("TWCorrCanv%d%d", npl, side),
                      Form("plane %s_%s", pl_names[npl].Data(), side_names[side].Data()), 1000, 700);

      // TWCorr_canv_FADC[npl][side] =
      //     new TCanvas(Form("TWCorrCanvFADC%d%d", npl, side),
      //                 Form("FADC plane %s_%s", pl_names[npl].Data(), side_names[side].Data()), 1000, 700);
      // Divide Canvas
      if (npl < 5) {
        TWUnCorr_canv[npl][side]->Divide(4, 3);
        TWCorr_canv[npl][side]->Divide(4, 3);
        // TWCorr_canv_FADC[npl][side]->Divide(4, 3);
      }

      // Loop over pmt
      Int_t fit_status;

      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        TWUnCorr_canv[npl][side]->cd(ipmt + 1);
        h2Hist_TW_UnCorr[npl][side][ipmt]->Draw("colz");

        TWCorr_canv[npl][side]->cd(ipmt + 1);
        h2Hist_TW_Corr[npl][side][ipmt]->Draw("colz");
        fit_status = -1;

        // Require ONLY one side
        if (side == 0) {

          lladhodo_cableArr[npl][ipmt] = h1Hist_TWDiff[npl][ipmt]->GetMean();
          lladhodo_sigArr[npl][ipmt]   = h1Hist_TWDiff[npl][ipmt]->GetStdDev();

          // Calculate the value at which h1Hist_TWDiff drops below 0.1 times its peak for both ends
          double peak      = h1Hist_TWDiff[npl][ipmt]->GetMaximum();
          double peak_bin  = h1Hist_TWDiff[npl][ipmt]->GetMaximumBin();
          double threshold = paddle_vel_fit_max_threshold * peak;
          int bin_left     = peak_bin;
          int bin_right    = peak_bin;

          // Find the bin to the left where the value drops below the threshold (main threshold)
          while (bin_left > 1 && h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_left) > threshold) {
            bin_left--;
          }
          double value_below_threshold_btm = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_left);

          // Find the bin to the right where the value drops below the threshold (main threshold)
          while (bin_right < h1Hist_TWDiff[npl][ipmt]->GetNbinsX() &&
                 h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_right) > threshold) {
            bin_right++;
          }
          double value_below_threshold_top = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_right);

          // Also compute the same pair of edge-values for the two error-bar threshold points
          // Build fractional thresholds, clamp into [0,1]
          double low_frac = paddle_vel_fit_max_threshold_err[0];
          if (low_frac < 0.0)
            low_frac = 0.0;
          double high_frac = paddle_vel_fit_max_threshold_err[1];
          if (high_frac > 1.0)
            high_frac = 1.0;

          // Convert fractional thresholds to counts (multiplying by peak)
          double threshold_low  = low_frac * peak;
          double threshold_high = high_frac * peak;

          // Find left/right bins for low threshold
          int bin_left_low  = peak_bin;
          int bin_right_low = peak_bin;
          while (bin_left_low > 1 && h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_left_low) > threshold_low) {
            --bin_left_low;
          }
          while (bin_right_low < h1Hist_TWDiff[npl][ipmt]->GetNbinsX() &&
                 h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_right_low) > threshold_low) {
            ++bin_right_low;
          }
          double value_below_threshold_btm_low = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_left_low);
          double value_below_threshold_top_low = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_right_low);

          // Find left/right bins for high threshold
          int bin_left_high  = peak_bin;
          int bin_right_high = peak_bin;
          while (bin_left_high > 1 && h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_left_high) > threshold_high) {
            --bin_left_high;
          }
          while (bin_right_high < h1Hist_TWDiff[npl][ipmt]->GetNbinsX() &&
                 h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_right_high) > threshold_high) {
            ++bin_right_high;
          }
          double value_below_threshold_btm_high = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_left_high);
          double value_below_threshold_top_high = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_right_high);

          // Compute velocities for main, low and high thresholds (protect against zero/very small denom)
          double delta_main = value_below_threshold_top - value_below_threshold_btm;
          double vel_main   = (fabs(delta_main) > 1e-6) ? (barlength[npl][ipmt] / delta_main) : 0.0;

          double delta_low = value_below_threshold_top_low - value_below_threshold_btm_low;
          double vel_low   = (fabs(delta_low) > 1e-6) ? (barlength[npl][ipmt] / delta_low) : 0.0;

          double delta_high = value_below_threshold_top_high - value_below_threshold_btm_high;
          double vel_high   = (fabs(delta_high) > 1e-6) ? (barlength[npl][ipmt] / delta_high) : 0.0;

          // Store main velocity as before; vel_low and vel_high are available for error-band use
          lladhodo_velArr[npl][ipmt] = vel_main;

          lladhodo_velArr_err[npl][ipmt][0] = vel_low;
          lladhodo_velArr_err[npl][ipmt][1] = vel_high;

          // Draw 1D TWAvg Histos
          TWAvg_canv[npl]->cd(ipmt + 1);
          h1Hist_TWAvg[npl][ipmt]->SetLineColor(kBlack);
          h1Hist_TWAvg_CUT[npl][ipmt]->SetLineColor(kRed);
          h1Hist_TWAvg[npl][ipmt]->Draw();
          h1Hist_TWAvg_CUT[npl][ipmt]->Draw("same");

          // Draw 1D TWDiff Histos
          TWDiff_canv[npl]->cd(ipmt + 1);
          h1Hist_TWDiff[npl][ipmt]->Draw();

          // Draw a line at lladhodo_cableArr[npl][ipmt]
          TLine *line_cable = new TLine(lladhodo_cableArr[npl][ipmt], 0, lladhodo_cableArr[npl][ipmt],
                                        h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_cable->SetLineColor(kBlack);
          line_cable->SetLineStyle(2);
          line_cable->Draw("same");

          // Draw lines at bin_left and bin_right
          TLine *line_left = new TLine(value_below_threshold_btm, 0, value_below_threshold_btm,
                                       h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_left->SetLineColor(kRed);
          line_left->SetLineStyle(2);
          line_left->Draw("same");

          TLine *line_right = new TLine(value_below_threshold_top, 0, value_below_threshold_top,
                                        h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_right->SetLineColor(kRed);
          line_right->SetLineStyle(2);
          line_right->Draw("same");

          // Repeate above for FADC
          // FADC: analogous processing to the TDC block above
          lladhodo_cableArr_FADC[npl][ipmt] = h1Hist_TWDiff_FADC[npl][ipmt]->GetMean();
          lladhodo_sigArr_FADC[npl][ipmt]   = h1Hist_TWDiff_FADC[npl][ipmt]->GetStdDev();

          // Determine velocity bounds from FADC diff histogram peak width at threshold
          double peak_fadc      = h1Hist_TWDiff_FADC[npl][ipmt]->GetMaximum();
          int peak_bin_fadc     = h1Hist_TWDiff_FADC[npl][ipmt]->GetMaximumBin();
          double threshold_fadc = paddle_vel_fit_max_threshold * peak_fadc;
          int bin_left_fadc     = peak_bin_fadc;
          int bin_right_fadc    = peak_bin_fadc;

          // Move left until below threshold (or hit first bin)
          while (bin_left_fadc > 1 && h1Hist_TWDiff_FADC[npl][ipmt]->GetBinContent(bin_left_fadc) > threshold_fadc) {
            --bin_left_fadc;
          }
          double value_below_threshold_btm_fadc = h1Hist_TWDiff_FADC[npl][ipmt]->GetBinCenter(bin_left_fadc);

          // Move right until below threshold (or hit last bin)
          while (bin_right_fadc < h1Hist_TWDiff_FADC[npl][ipmt]->GetNbinsX() &&
                 h1Hist_TWDiff_FADC[npl][ipmt]->GetBinContent(bin_right_fadc) > threshold_fadc) {
            ++bin_right_fadc;
          }
          double value_below_threshold_top_fadc = h1Hist_TWDiff_FADC[npl][ipmt]->GetBinCenter(bin_right_fadc);

          // Protect against zero division
          double delta_fadc = value_below_threshold_top_fadc - value_below_threshold_btm_fadc;
          if (fabs(delta_fadc) > 1e-6) {
            lladhodo_velArr_FADC[npl][ipmt] = barlength[npl][ipmt] / delta_fadc;
          } else {
            lladhodo_velArr_FADC[npl][ipmt] = 0.0;
          }

          // Draw FADC average histos on the 2D/alternate canvas
          TWAvg_canv_FADC[npl]->cd(ipmt + 1);
          h1Hist_TWAvg_FADC[npl][ipmt]->SetLineColor(kBlue);
          h1Hist_TWAvg_CUT_FADC[npl][ipmt]->SetLineColor(kMagenta);
          h1Hist_TWAvg_FADC[npl][ipmt]->Draw();
          h1Hist_TWAvg_CUT_FADC[npl][ipmt]->Draw("same");

          // Draw FADC diff histogram and vertical lines on the TWDiff canvas
          TWDiff_canv_FADC[npl]->cd(ipmt + 1);
          h1Hist_TWDiff_FADC[npl][ipmt]->Draw();

          // Vertical line at cable offset (FADC)
          TLine *line_cable_fadc = new TLine(lladhodo_cableArr_FADC[npl][ipmt], 0, lladhodo_cableArr_FADC[npl][ipmt],
                                             h1Hist_TWDiff_FADC[npl][ipmt]->GetMaximum());
          line_cable_fadc->SetLineColor(kBlue);
          line_cable_fadc->SetLineStyle(2);
          line_cable_fadc->Draw("same");

          // Lines at threshold-determined edges
          TLine *line_left_fadc = new TLine(value_below_threshold_btm_fadc, 0, value_below_threshold_btm_fadc,
                                            h1Hist_TWDiff_FADC[npl][ipmt]->GetMaximum());
          line_left_fadc->SetLineColor(kGreen + 2);
          line_left_fadc->SetLineStyle(2);
          line_left_fadc->Draw("same");

          TLine *line_right_fadc = new TLine(value_below_threshold_top_fadc, 0, value_below_threshold_top_fadc,
                                             h1Hist_TWDiff_FADC[npl][ipmt]->GetMaximum());
          line_right_fadc->SetLineColor(kGreen + 2);
          line_right_fadc->SetLineStyle(2);
          line_right_fadc->Draw("same");

        } // end single SIDE requirement

      } // end pmt loop

    } // end side loop

  } // end plane loop

  /***********DRAW CANVAS FOR lladhodo_cableArr***************/
  TCanvas *cableArr_canv = new TCanvas("cableArr_canv", "Cable Time Offsets per Paddle", 1000, 700);
  cableArr_canv->Divide(3, 2);

  TCanvas *cableArr_FADC_canv = new TCanvas("cableArr_FADC_canv", "FADC Cable Time Offsets per Paddle", 1000, 700);
  cableArr_FADC_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    // Build TGraphErrors for TDC cable offsets (analogous to FADC below)
    Int_t npts = maxPMT[npl];
    double *x  = new double[npts];
    double *y  = new double[npts];
    double *ex = new double[npts];
    double *ey = new double[npts];

    for (Int_t ipmt = 0; ipmt < npts; ++ipmt) {
      x[ipmt]  = ipmt + 1;
      y[ipmt]  = lladhodo_cableArr[npl][ipmt];
      ex[ipmt] = 0.0;
      ey[ipmt] = lladhodo_sigArr[npl][ipmt];
    }

    TGraphErrors *g_CableArr = new TGraphErrors(npts, x, y, ex, ey);
    g_CableArr->SetTitle(Form("Cable Time Offsets per Paddle for Plane %s", pl_names[npl].Data()));
    g_CableArr->GetXaxis()->SetTitle("Paddle Number");
    g_CableArr->GetYaxis()->SetTitle("Cable Time Offset (ns)");
    g_CableArr->SetMarkerStyle(20);    // Set marker style to points
    g_CableArr->SetMarkerColor(kBlue); // Set marker color

    double minY = *min_element(lladhodo_cableArr[npl], lladhodo_cableArr[npl] + maxPMT[npl]);
    double maxY = *max_element(lladhodo_cableArr[npl], lladhodo_cableArr[npl] + maxPMT[npl]);
    g_CableArr->GetYaxis()->SetRangeUser(minY - 1, maxY + 1); // Adjust y range

    cableArr_canv->cd(npl + 1);
    g_CableArr->Draw("AP"); // draw axes and points with error bars

    delete[] x;
    delete[] y;
    delete[] ex;
    delete[] ey;

    // Draw a dashed line at y = 0
    TLine *line = new TLine(0, 0, maxPMT[npl] + 1, 0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2); // Dashed line
    line->Draw("same");

    // Build arrays for TGraphErrors (x, y, ex, ey) so points have y error bars = sigArr_FADC
    npts            = maxPMT[npl];
    double *x_fadc  = new double[npts];
    double *y_fadc  = new double[npts];
    double *ex_fadc = new double[npts];
    double *ey_fadc = new double[npts];

    for (Int_t ipmt = 0; ipmt < npts; ipmt++) {
      x_fadc[ipmt]  = ipmt + 1;
      y_fadc[ipmt]  = lladhodo_cableArr_FADC[npl][ipmt];
      ex_fadc[ipmt] = 0.0;                             // no x error
      ey_fadc[ipmt] = lladhodo_sigArr_FADC[npl][ipmt]; // y error from sigArr_FADC
    }

    TGraphErrors *g_CableArr_FADC = new TGraphErrors(npts, x_fadc, y_fadc, ex_fadc, ey_fadc);
    g_CableArr_FADC->SetTitle(Form("FADC Cable Time Offsets per Paddle for Plane %s", pl_names[npl].Data()));
    g_CableArr_FADC->GetXaxis()->SetTitle("Paddle Number");
    g_CableArr_FADC->GetYaxis()->SetTitle("Cable Time Offset (ns)");
    g_CableArr_FADC->SetMarkerStyle(21);
    g_CableArr_FADC->SetMarkerColor(kMagenta);

    double minY_fadc = *min_element(lladhodo_cableArr_FADC[npl], lladhodo_cableArr_FADC[npl] + maxPMT[npl]);
    double maxY_fadc = *max_element(lladhodo_cableArr_FADC[npl], lladhodo_cableArr_FADC[npl] + maxPMT[npl]);
    g_CableArr_FADC->GetYaxis()->SetRangeUser(minY_fadc - 1, maxY_fadc + 1);

    // cleanup temporary arrays
    delete[] x_fadc;
    delete[] y_fadc;
    delete[] ex_fadc;
    delete[] ey_fadc;
    cableArr_canv->cd(npl + 1);
    g_CableArr_FADC->Draw("Psame");

    TLine *line_fadc = new TLine(0, 0, maxPMT[npl] + 1, 0);
    line_fadc->SetLineColor(kBlack);
    line_fadc->SetLineStyle(2);
    line_fadc->Draw("same");
  }

  /***********DRAW CANVAS FOR lladhodo_velArr***************/
  TCanvas *velArr_canv = new TCanvas("velArr_canv", "Propagation Velocities per Paddle", 1000, 700);
  velArr_canv->Divide(3, 2);

  TCanvas *velArr_FADC_canv = new TCanvas("velArr_FADC_canv", "FADC Propagation Velocities per Paddle", 1000, 700);
  velArr_FADC_canv->Divide(3, 2);

  // Calculate global average velocity (TDC) excluding values >18 or <15
  double sum_vel = 0.0;
  int count_vel = 0;
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      double val = lladhodo_velArr[npl][ipmt];
      if (val >= 15.0 && val <= 18.0) {
        sum_vel += val;
        count_vel++;
      }
    }
  }
  double avg_vel = (count_vel > 0) ? (sum_vel / count_vel) : 0.0;
  cout << "Average TDC velocity (14-18 cm/ns): " << avg_vel << " cm/ns (" << count_vel << " values)" << endl;

  // Calculate global average velocity (FADC) excluding values >18 or <14
  double sum_vel_fadc = 0.0;
  int count_vel_fadc = 0;
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      double val = lladhodo_velArr_FADC[npl][ipmt];
      if (val >= 14.0 && val <= 18.0) {
        sum_vel_fadc += val;
        count_vel_fadc++;
      }
    }
  }
  double avg_vel_fadc = (count_vel_fadc > 0) ? (sum_vel_fadc / count_vel_fadc) : 0.0;
  cout << "Average FADC velocity (14-18 cm/ns): " << avg_vel_fadc << " cm/ns (" << count_vel_fadc << " values)" << endl;

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    Int_t npts = maxPMT[npl];
    double *x = new double[npts];
    double *y = new double[npts];
    double *ex = new double[npts];
    double *eylow = new double[npts];
    double *eyhigh = new double[npts];

    for (Int_t ipmt = 0; ipmt < npts; ipmt++) {
      x[ipmt] = ipmt + 1;
      y[ipmt] = lladhodo_velArr[npl][ipmt];
      ex[ipmt] = 0.0;
      eylow[ipmt] = lladhodo_velArr[npl][ipmt] - lladhodo_velArr_err[npl][ipmt][0];
      eyhigh[ipmt] = lladhodo_velArr_err[npl][ipmt][1] - lladhodo_velArr[npl][ipmt];
    }

    TGraphAsymmErrors *g_VelArr = new TGraphAsymmErrors(npts, x, y, ex, ex, eylow, eyhigh);
    g_VelArr->SetTitle(Form("Propagation Velocities per Paddle for Plane %s", pl_names[npl].Data()));
    g_VelArr->GetXaxis()->SetTitle("Paddle Number");
    g_VelArr->GetYaxis()->SetTitle("Propagation Velocity (cm/ns)");
    g_VelArr->SetMarkerStyle(20);
    g_VelArr->SetMarkerColor(kRed);

    double minY = *min_element(lladhodo_velArr[npl], lladhodo_velArr[npl] + maxPMT[npl]);
    double maxY = *max_element(lladhodo_velArr[npl], lladhodo_velArr[npl] + maxPMT[npl]);
    g_VelArr->GetYaxis()->SetRangeUser(minY - 1, maxY + 1);

    velArr_canv->cd(npl + 1);
    g_VelArr->Draw("AP");

    // Draw horizontal line at global average velocity
    TLine *line_avg_vel = new TLine(0, avg_vel, maxPMT[npl] + 1, avg_vel);
    line_avg_vel->SetLineColor(kGreen + 2);
    line_avg_vel->SetLineStyle(2);
    line_avg_vel->SetLineWidth(2);
    line_avg_vel->Draw("same");

    delete[] x;
    delete[] y;
    delete[] ex;
    delete[] eylow;
    delete[] eyhigh;

    TGraph *g_VelArr_FADC = new TGraph(maxPMT[npl]);
    g_VelArr_FADC->SetTitle(Form("FADC Propagation Velocities per Paddle for Plane %s", pl_names[npl].Data()));
    g_VelArr_FADC->GetXaxis()->SetTitle("Paddle Number");
    g_VelArr_FADC->GetYaxis()->SetTitle("Propagation Velocity (cm/ns)");
    g_VelArr_FADC->SetMarkerStyle(21);
    g_VelArr_FADC->SetMarkerColor(kBlue);

    double minY_fadc = *min_element(lladhodo_velArr_FADC[npl], lladhodo_velArr_FADC[npl] + maxPMT[npl]);
    double maxY_fadc = *max_element(lladhodo_velArr_FADC[npl], lladhodo_velArr_FADC[npl] + maxPMT[npl]);
    g_VelArr_FADC->GetYaxis()->SetRangeUser(minY_fadc - 1, maxY_fadc + 1);

    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      g_VelArr_FADC->SetPoint(ipmt, ipmt + 1, lladhodo_velArr_FADC[npl][ipmt]);
    }
    velArr_FADC_canv->cd(npl + 1);
    g_VelArr_FADC->Draw("AP");

    // Draw horizontal line at global average velocity (FADC)
    TLine *line_avg_vel_fadc = new TLine(0, avg_vel_fadc, maxPMT[npl] + 1, avg_vel_fadc);
    line_avg_vel_fadc->SetLineColor(kGreen + 2);
    line_avg_vel_fadc->SetLineStyle(2);
    line_avg_vel_fadc->SetLineWidth(2);
    line_avg_vel_fadc->Draw("same");

    // Draw a dashed line at y = 0
    TLine *line_fadc_zero = new TLine(0, 0, maxPMT[npl] + 1, 0);
    line_fadc_zero->SetLineColor(kBlack);
    line_fadc_zero->SetLineStyle(2);
    line_fadc_zero->Draw("same");
  }

  /************WRITE FIT RESULTS TO PARAMETER FILE***************/

  ofstream outPARAM;
  outPARAM.open(Form("../../PARAM/LAD/HODO/ladhodo_Vpcalib_%d.param", runNUM));

  outPARAM << "; LAD Hodoscope Parameter File Containing propagation velocities per paddle " << endl;
  outPARAM << "; and signal cable time diff. offsets per paddle " << endl;
  outPARAM << Form("; Run %d ", runNUM) << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";Propagation Velocity Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_velFit = ";

  //--------Write Velocity Parameters to Param File-----
  outPARAM << fixed << setprecision(4);
  if (use_vp_avg) {
    // Write average velocity for all paddles
    for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
      for (Int_t npl = 0; npl < NPLANES; npl++) {
        if (npl == 0) {
          if (ipmt == 0) {
            outPARAM << avg_vel;
          } else {
            outPARAM << setw(26) << avg_vel;
          }
        } else {
          outPARAM << ", " << setw(15) << avg_vel;
        }
      }
      outPARAM << endl;
    }
  } else {
    // Write individual velocity values
    for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
      for (Int_t npl = 0; npl < NPLANES; npl++) {
        if (npl == 0) {
          if (ipmt == 0) {
            outPARAM << lladhodo_velArr[npl][ipmt];
          } else {
            outPARAM << setw(26) << lladhodo_velArr[npl][ipmt];
          }
        } else {
          outPARAM << ", " << setw(15) << lladhodo_velArr[npl][ipmt];
        }
      }
      outPARAM << endl;
    }
  }

  // Also write FADC velocities
  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << ";FADC Propagation Velocity Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_velFit_FADC = ";
  outPARAM << fixed << setprecision(4);
  if (use_vp_avg) {
    // Write average FADC velocity for all paddles
    for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
      for (Int_t npl = 0; npl < NPLANES; npl++) {
        if (npl == 0) {
          if (ipmt == 0) {
            outPARAM << avg_vel_fadc;
          } else {
            outPARAM << setw(26) << avg_vel_fadc;
          }
        } else {
          outPARAM << ", " << setw(15) << avg_vel_fadc;
        }
      }
      outPARAM << endl;
    }
  } else {
    // Write individual FADC velocity values
    for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
      for (Int_t npl = 0; npl < NPLANES; npl++) {
        if (npl == 0) {
          if (ipmt == 0) {
            outPARAM << lladhodo_velArr_FADC[npl][ipmt];
          } else {
            outPARAM << setw(26) << lladhodo_velArr_FADC[npl][ipmt];
          }
        } else {
          outPARAM << ", " << setw(15) << lladhodo_velArr_FADC[npl][ipmt];
        }
      }
      outPARAM << endl;
    }
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";PMTs Signal Cable Time Diff. Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_cableFit = ";

  // Write Cable Length Time Diff. Parameters to Param File
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_cableArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_cableArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_cableArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  // Also write FADC cable offsets
  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << ";FADC PMTs Signal Cable Time Diff. Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_cableFit_FADC = ";
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_cableArr_FADC[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_cableArr_FADC[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_cableArr_FADC[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";PMTs Time Diff. Sigma Parameters" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_TopSigma = ";

  // Write Sigma Parameters to file (TDC)
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_sigArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_sigArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_sigArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << "lladhodo_BtmSigma = ";

  // Write Btm Sigma Parameters to file (use same array if separate not available)
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_sigArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_sigArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_sigArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  // Also write FADC sigma parameters (Top & Btm)
  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << ";FADC PMTs Time Diff. Sigma Parameters" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_TopSigma_FADC = ";
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_sigArr_FADC[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_sigArr_FADC[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_sigArr_FADC[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << "lladhodo_BtmSigma_FADC = ";
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_sigArr_FADC[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_sigArr_FADC[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_sigArr_FADC[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  cout << "FINISHED Getting Vp and Cable Fits . . . " << endl;
  cout << "Starting the code to fit Hodo Matrix . . . " << endl;

  /**********BEGIN CODE TO FIT HODO MATRIX**************/

  // Initialize h1Hist_TWCorr_BarLCoef histograms
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      h1Hist_TWCorr_BarLCoef[npl][ipmt] =
          new TH1F(Form("TWCorr_BarLCoef_%s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("TWCorr_BarLCoef, plane %s, paddle %d", pl_names[npl].Data(), ipmt + 1), DiffTime_NBINS,
                   DiffTime_MIN, DiffTime_MAX);
      // Define FADC version of the BarLCoef histogram (analogous to TDC one)
      h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt] =
          new TH1F(Form("TWCorr_BarLCoef_FADC_%s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("TWCorr_BarLCoef FADC, plane %s, paddle %d", pl_names[npl].Data(), ipmt + 1), DiffTime_NBINS,
                   DiffTime_MIN, DiffTime_MAX);
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetXaxis()->SetTitle("FADC-corrected Δt (ns)");
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetYaxis()->SetTitle("Entries");
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetXaxis()->CenterTitle();
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetYaxis()->CenterTitle();
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->SetLineColor(kBlue);
      // h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->SetFillColor(0);
    }
  }
  Double_t good_TW_top, good_TW_btm, good_TW_top_ref, good_TW_btm_ref;
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over pmt
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
        if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {

          good_TW_top = TdcTimeTWCorr[npl][0][ipmt];
          good_TW_btm = TdcTimeTWCorr[npl][1][ipmt] - 2 * lladhodo_cableArr[npl][ipmt]; // IMPORTANT: Apply cable time
          // correction obtained from fits
          good_TW_top_ref = TdcTimeTWCorr[refPlane][0][refBar];
          good_TW_btm_ref = TdcTimeTWCorr[refPlane][1][refBar] - 2 * lladhodo_cableArr[refPlane][refBar];
          if (good_TW_top_ref == 0) {
            continue;
          }
          h1Hist_TWCorr_BarLCoef[npl][ipmt]->Fill(((good_TW_top_ref + good_TW_btm_ref) - (good_TW_top + good_TW_btm)) /
                                                  2);

          // Fill FADC version of the BarLCoef histogram (analogous to TDC one)
          if (AdcPulseTime[npl][0][ipmt] < AdcTimeTWCorr_MAX && AdcPulseTime[npl][1][ipmt] < AdcTimeTWCorr_MAX) {
            Double_t good_FADC_top = AdcPulseTime[npl][0][ipmt];
            Double_t good_FADC_btm =
                AdcPulseTime[npl][1][ipmt] - 2 * lladhodo_cableArr_FADC[npl][ipmt]; // Apply cable time
            // correction obtained from fits
            Double_t good_FADC_top_ref = AdcPulseTime[refPlane][0][refBar];
            Double_t good_FADC_btm_ref =
                AdcPulseTime[refPlane][1][refBar] - 2 * lladhodo_cableArr_FADC[refPlane][refBar];
            if (good_FADC_top_ref == 0) {
              continue;
            }
            h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->Fill(
                ((good_FADC_top_ref + good_FADC_btm_ref) - (good_FADC_top + good_FADC_btm)) / 2);

          } // end time cut
        }
      } // end pmt loop

    } // end plane loop
    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";
  } // end event loop

  // Loop through h1Hist_TWCorr_BarLCoef and fill the mean into lladhodo_LCoeff
  double refBarTime     = h1Hist_TWCorr_BarLCoef[refPlane][refBar]->GetMean();
  double refBarTimeFADC = h1Hist_TWCorr_BarLCoef_FADC[refPlane][refBar]->GetMean();
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      // if (npl == refPlane && ipmt == refBar) {
      //   continue;
      // }
      lladhodo_LCoeff[npl][ipmt]             = (h1Hist_TWCorr_BarLCoef[npl][ipmt]->GetMean());
      lladhodo_sigArr_LCoeff[npl][ipmt]      = (h1Hist_TWCorr_BarLCoef[npl][ipmt]->GetStdDev());
      lladhodo_LCoeff_FADC[npl][ipmt]        = (h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetMean());
      lladhodo_sigArr_LCoeff_FADC[npl][ipmt] = (h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetStdDev());
    }
  }
  // Create Canvas to store TW-Corr Time/Dist vs. trk position
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Create Canvas for TDC BarLCoef
    TWCorr_BarLCoef[npl] = new TCanvas(Form("TWCorr_BarLCoef_%d", npl),
                                       Form("TWCorr_BarLCoef, plane %s", pl_names[npl].Data()), 1000, 700);

    // Divide Canvas
    if (npl < 5) {
      TWCorr_BarLCoef[npl]->Divide(4, 3);
    }

    // Create Canvas for FADC BarLCoef
    TWCorr_BarLCoef_FADC[npl] = new TCanvas(Form("TWCorr_BarLCoef_FADC_%d", npl),
                                            Form("TWCorr_BarLCoef FADC, plane %s", pl_names[npl].Data()), 1000, 700);
    if (npl < 5) {
      TWCorr_BarLCoef_FADC[npl]->Divide(4, 3);
    }

    // Loop over pmt
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

      // Draw TDC BarLCoef histogram
      TWCorr_BarLCoef[npl]->cd(ipmt + 1);
      h1Hist_TWCorr_BarLCoef[npl][ipmt]->Draw();

      // Compute lladhodo_LCoeff relative to reference (TDC)
      lladhodo_LCoeff[npl][ipmt] =
          h1Hist_TWCorr_BarLCoef[npl][ipmt]->GetMean() - h1Hist_TWCorr_BarLCoef[refPlane][refBar]->GetMean();

      // Draw FADC BarLCoef histogram
      TWCorr_BarLCoef_FADC[npl]->cd(ipmt + 1);
      h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->Draw();

      // Compute lladhodo_LCoeff_FADC relative to reference (FADC)
      lladhodo_LCoeff_FADC[npl][ipmt] =
          h1Hist_TWCorr_BarLCoef_FADC[npl][ipmt]->GetMean() - h1Hist_TWCorr_BarLCoef_FADC[refPlane][refBar]->GetMean();

    } // end pmt loop

  } // end plane loop

  /***********DRAW CANVAS FOR lladhodo_LCoeff (TDC)***************/
  TCanvas *LCoeff_canv = new TCanvas("LCoeff_canv", "Timing Corrections per Paddle (TDC)", 1000, 700);
  LCoeff_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    Int_t npts = maxPMT[npl];
    double *x  = new double[npts];
    double *y  = new double[npts];
    double *ex = new double[npts];
    double *ey = new double[npts];

    for (Int_t ipmt = 0; ipmt < npts; ++ipmt) {
      x[ipmt]  = ipmt + 1;
      y[ipmt]  = lladhodo_LCoeff[npl][ipmt];
      ex[ipmt] = 0.0;
      ey[ipmt] = lladhodo_sigArr_LCoeff[npl][ipmt];
    }

    TGraphErrors *g_LCoeff = new TGraphErrors(npts, x, y, ex, ey);
    g_LCoeff->SetTitle(Form("Timing Corrections per Paddle for Plane %s (TDC)", pl_names[npl].Data()));
    g_LCoeff->GetXaxis()->SetTitle("Paddle Number");
    g_LCoeff->GetYaxis()->SetTitle("Timing Correction (ns)");
    g_LCoeff->SetMarkerStyle(20);
    g_LCoeff->SetMarkerColor(kRed);

    double minY = *min_element(lladhodo_LCoeff[npl], lladhodo_LCoeff[npl] + npts);
    double maxY = *max_element(lladhodo_LCoeff[npl], lladhodo_LCoeff[npl] + npts);
    g_LCoeff->GetYaxis()->SetRangeUser(minY - 1, maxY + 1);

    LCoeff_canv->cd(npl + 1);
    g_LCoeff->Draw("AP"); // axes, points and error bars

    // dashed zero line
    TLine *zero_line = new TLine(0, 0, npts + 1, 0);
    zero_line->SetLineColor(kBlack);
    zero_line->SetLineStyle(2);
    zero_line->Draw("same");

    delete[] x;
    delete[] y;
    delete[] ex;
    delete[] ey;
  }

  /***********DRAW CANVAS FOR lladhodo_LCoeff_FADC (FADC)***************/
  TCanvas *LCoeff_FADC_canv = new TCanvas("LCoeff_FADC_canv", "Timing Corrections per Paddle (FADC)", 1000, 700);
  LCoeff_FADC_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    Int_t npts      = maxPMT[npl];
    double *x_fadc  = new double[npts];
    double *y_fadc  = new double[npts];
    double *ex_fadc = new double[npts];
    double *ey_fadc = new double[npts];

    for (Int_t ipmt = 0; ipmt < npts; ++ipmt) {
      x_fadc[ipmt]  = ipmt + 1;
      y_fadc[ipmt]  = lladhodo_LCoeff_FADC[npl][ipmt];
      ex_fadc[ipmt] = 0.0;
      ey_fadc[ipmt] = lladhodo_sigArr_LCoeff_FADC[npl][ipmt];
    }

    TGraphErrors *g_LCoeff_FADC = new TGraphErrors(npts, x_fadc, y_fadc, ex_fadc, ey_fadc);
    g_LCoeff_FADC->SetTitle(Form("Timing Corrections per Paddle for Plane %s (FADC)", pl_names[npl].Data()));
    g_LCoeff_FADC->GetXaxis()->SetTitle("Paddle Number");
    g_LCoeff_FADC->GetYaxis()->SetTitle("Timing Correction (ns)");
    g_LCoeff_FADC->SetMarkerStyle(21);
    g_LCoeff_FADC->SetMarkerColor(kBlue);

    double minY_fadc = *min_element(lladhodo_LCoeff_FADC[npl], lladhodo_LCoeff_FADC[npl] + npts);
    double maxY_fadc = *max_element(lladhodo_LCoeff_FADC[npl], lladhodo_LCoeff_FADC[npl] + npts);
    g_LCoeff_FADC->GetYaxis()->SetRangeUser(minY_fadc - 1, maxY_fadc + 1);

    LCoeff_canv->cd(npl + 1);
    g_LCoeff_FADC->Draw("Psame"); // axes, points and error bars

    TLine *zero_line_fadc = new TLine(0, 0, npts + 1, 0);
    zero_line_fadc->SetLineColor(kBlack);
    zero_line_fadc->SetLineStyle(2);
    zero_line_fadc->Draw("same");

    delete[] x_fadc;
    delete[] y_fadc;
    delete[] ex_fadc;
    delete[] ey_fadc;
  }

  // Write Fit Results to Parameter File

  outPARAM << "" << endl;
  outPARAM << "" << endl;
  outPARAM << "" << endl;
  outPARAM << ";Timing Corrections Per Paddle, where 000 Paddle 00 has been set as the reference paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_LCoeff = ";

  // Write Lambda Time Coeff. Parameters to Param File (TDC)
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_LCoeff[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_LCoeff[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_LCoeff[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  // Also write the FADC timing corrections
  outPARAM << "" << endl;
  outPARAM << "" << endl;
  outPARAM << ";Timing Corrections Per Paddle from FADC, where 000 Paddle 00 has been set as the reference paddle"
           << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_LCoeff_FADC = ";

  // Write Lambda Time Coeff. Parameters to Param File (FADC)
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_LCoeff_FADC[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_LCoeff_FADC[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_LCoeff_FADC[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  cout << "FINISHED Fitting Hodo Matrix . . . " << endl;
  cout << "Parameter File Created: " << Form("../../PARAM/LAD/HODO/ladhodo_Vpcalib_%d.param", runNUM) << endl;

  // Write Histograms to ROOT file
  // Create output root file where histograms will be stored
  TFile *outROOT = new TFile(Form("HodoCalibPlots_%d.root", runNUM), "recreate");
  // Create directories for each type of canvas
  outROOT->mkdir("TWAvg");
  outROOT->mkdir("TWAvg_FADC");
  outROOT->mkdir("TWAvg2D");
  outROOT->mkdir("TWDiff");
  outROOT->mkdir("TWDiff_FADC");
  outROOT->mkdir("TWCorr_BarLCoef");
  outROOT->mkdir("TWUnCorr");
  outROOT->mkdir("TWCorr");
  outROOT->mkdir("Param");

  // Loop through all canvases and save each canvas in the appropriate directory
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outROOT->cd("TWAvg");
    TWAvg_canv[npl]->Write();
    outROOT->cd("TWAvg_FADC");
    TWAvg_canv_FADC[npl]->Write();
    outROOT->cd("TWAvg2D");
    TWAvg_canv_2D[npl]->Write();
    outROOT->cd("TWDiff");
    TWDiff_canv[npl]->Write();
    outROOT->cd("TWDiff_FADC");
    TWDiff_canv_FADC[npl]->Write();
    outROOT->cd("TWCorr_BarLCoef");
    TWCorr_BarLCoef[npl]->Write();
    for (Int_t side = 0; side < SIDES; side++) {
      outROOT->cd("TWUnCorr");
      TWUnCorr_canv[npl][side]->Write();
      outROOT->cd("TWCorr");
      TWCorr_canv[npl][side]->Write();
    }
  }
  outROOT->cd("Param");
  cableArr_canv->Write();
  // cableArr_FADC_canv->Write();
  velArr_canv->Write();
  velArr_FADC_canv->Write();
  LCoeff_canv->Write();
  outROOT->Write();
  outROOT->Close();
}
