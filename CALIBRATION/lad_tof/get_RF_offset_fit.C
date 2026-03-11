#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
double mod_value = 4.00801; // RF period in ns

void processFitting(int run_num, const vector<int> &run_nums, vector<double> &rf_offsets_P,
                    vector<double> &rf_offset_errs_P, vector<double> &rf_widths_P, vector<double> &rf_width_errs_P,
                    vector<double> &rf_offsets_H, vector<double> &rf_offset_errs_H, vector<double> &rf_widths_H,
                    vector<double> &rf_width_errs_H, TFile *inFile, TFile *out) {
  gROOT->SetBatch(kTRUE);

  vector<double> fit_ranges    = {0.2, 0.25, 0.3}; // Multiple fit ranges to test
  double fit_range_for_display = fit_ranges[fit_ranges.size() / 2];

  // Read histograms from input file
  TH1F *h_mod_P = (TH1F *)inFile->Get(Form("h_mod_P/h_mod_P_%d", run_num));
  TH1F *h_mod_H = (TH1F *)inFile->Get(Form("h_mod_H/h_mod_H_%d", run_num));
  TH1F *h_P     = (TH1F *)inFile->Get(Form("h_P/h_P_%d", run_num));
  TH1F *h_H     = (TH1F *)inFile->Get(Form("h_H/h_H_%d", run_num));

  if (!h_mod_P || !h_mod_H || !h_P || !h_H) {
    cout << "Warning: Could not read all histograms for run " << run_num << endl;
    return;
  }

  // Extended histograms for fitting: range [-N, 2N]
  TH1F *h_mod_P_ext = new TH1F(Form("h_mod_P_ext_%d", run_num), "RF Offset Mod P Ext", 1500, -mod_value, 2 * mod_value);
  TH1F *h_mod_H_ext = new TH1F(Form("h_mod_H_ext_%d", run_num), "RF Offset Mod H Ext", 1500, -mod_value, 2 * mod_value);

  // Fill extended histograms by tiling the modulo histogram content
  for (int i = 1; i <= h_mod_P->GetNbinsX(); i++) {
    double x = h_mod_P->GetBinCenter(i);
    double c = h_mod_P->GetBinContent(i);
    if (c == 0) {
      continue;
    }
    h_mod_P_ext->Fill(x, c);
    h_mod_P_ext->Fill(x - mod_value, c);
    h_mod_P_ext->Fill(x + mod_value, c);
  }

  for (int i = 1; i <= h_mod_H->GetNbinsX(); i++) {
    double x = h_mod_H->GetBinCenter(i);
    double c = h_mod_H->GetBinContent(i);
    if (c == 0) {
      continue;
    }
    h_mod_H_ext->Fill(x, c);
    h_mod_H_ext->Fill(x - mod_value, c);
    h_mod_H_ext->Fill(x + mod_value, c);
  }

  // Create canvas for P arm modulo histogram
  TCanvas *c_mod_P = new TCanvas(Form("c_mod_P_%d", run_num), "RF Offset Mod P", 800, 600);
  h_mod_P->Draw();

  // Find the peak position
  int max_bin    = h_mod_P->GetMaximumBin();
  double max_pos = h_mod_P->GetBinCenter(max_bin);

  // Perform multiple fits with different fit ranges
  vector<double> P_means, P_sigmas, P_norms;
  vector<double> P_fit_lows, P_fit_highs;

  for (double fit_range : fit_ranges) {
    // Calculate fit range ±fit_range around peak
    double fit_low_P  = max_pos - fit_range;
    double fit_high_P = max_pos + fit_range;
    // Standard Gaussian fit on extended histogram (no wrapping)
    TF1 *gaus_P_fit = new TF1(Form("gaus_P_fit_%d_%.2f", run_num, fit_range), "gaus", -mod_value, 2 * mod_value);

    // Set initial parameters
    gaus_P_fit->SetParameter(0, h_mod_P->GetMaximum());
    gaus_P_fit->SetParameter(1, max_pos);
    gaus_P_fit->SetParameter(2, 0.2);

    // Fit over the appropriate range(s)
    h_mod_P_ext->Fit(gaus_P_fit, "QN", "", fit_low_P, fit_high_P);

    double fit_mean  = gaus_P_fit->GetParameter(1);
    double fit_sigma = gaus_P_fit->GetParameter(2);

    // Wrap mean if needed
    if (fit_mean < 0) {
      fit_mean += mod_value;
    } else if (fit_mean > mod_value) {
      fit_mean -= mod_value;
    }

    P_means.push_back(fit_mean);
    P_sigmas.push_back(fit_sigma);
    P_norms.push_back(gaus_P_fit->GetParameter(0));
    P_fit_lows.push_back(fit_low_P);
    P_fit_highs.push_back(fit_high_P);
  }

  // Calculate average and std of means
  double rf_offset_P = 0, sigma_P = 0;
  for (double m : P_means) {
    rf_offset_P += m;
  }
  rf_offset_P /= P_means.size();

  // Calculate std of means as uncertainty
  double rf_offset_P_err = 0;
  for (double m : P_means) {
    rf_offset_P_err += (m - rf_offset_P) * (m - rf_offset_P);
  }
  rf_offset_P_err = sqrt(rf_offset_P_err / P_means.size());

  // Average sigma
  for (double s : P_sigmas) {
    sigma_P += s;
  }
  sigma_P /= P_sigmas.size();

  // Set histogram range to ±fit_range_for_display around peak (for display only)
  double display_low  = max_pos - fit_range_for_display;
  double display_high = max_pos + fit_range_for_display;

  // Handle wrapping for display range
  if (display_low < 0) {
    h_mod_P->GetXaxis()->SetRangeUser(0, display_high);
  } else if (display_high > mod_value) {
    h_mod_P->GetXaxis()->SetRangeUser(display_low, mod_value);
  } else {
    h_mod_P->GetXaxis()->SetRangeUser(display_low, display_high);
  }

  // Draw all fit functions on canvas, only in their fitted ranges
  int color_index = 0;
  int colors[]    = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2};
  for (size_t i = 0; i < P_means.size(); i++) {
    double fit_low  = P_fit_lows[i];
    double fit_high = P_fit_highs[i];
    double mean     = P_means[i];
    double sigma    = P_sigmas[i];
    double norm     = P_norms[i];

    int color = colors[color_index % 5];

    if (fit_low < 0) {
      // Draw wrapped segment near upper edge
      TF1 *gaus_P_display1 =
          new TF1(Form("gaus_P_display1_%d_%zu", run_num, i), "gaus", fit_low + mod_value, mod_value);
      gaus_P_display1->SetParameters(norm, mean + mod_value, sigma);
      gaus_P_display1->SetLineColor(color);
      gaus_P_display1->SetLineWidth(2);
      gaus_P_display1->Draw("same");

      TF1 *gaus_P_display2 = new TF1(Form("gaus_P_display2_%d_%zu", run_num, i), "gaus", 0, fit_high);
      gaus_P_display2->SetParameters(norm, mean, sigma);
      gaus_P_display2->SetLineColor(color);
      gaus_P_display2->SetLineWidth(2);
      gaus_P_display2->Draw("same");
    } else if (fit_high > mod_value) {
      // Draw wrapped segment near lower edge
      TF1 *gaus_P_display1 = new TF1(Form("gaus_P_display1_%d_%zu", run_num, i), "gaus", fit_low, mod_value);
      gaus_P_display1->SetParameters(norm, mean, sigma);
      gaus_P_display1->SetLineColor(color);
      gaus_P_display1->SetLineWidth(2);
      gaus_P_display1->Draw("same");

      TF1 *gaus_P_display2 = new TF1(Form("gaus_P_display2_%d_%zu", run_num, i), "gaus", 0, fit_high - mod_value);
      gaus_P_display2->SetParameters(norm, mean - mod_value, sigma);
      gaus_P_display2->SetLineColor(color);
      gaus_P_display2->SetLineWidth(2);
      gaus_P_display2->Draw("same");
    } else {
      // Normal case: draw only in fitted range
      TF1 *gaus_P_display = new TF1(Form("gaus_P_display_%d_%zu", run_num, i), "gaus", fit_low, fit_high);
      gaus_P_display->SetParameters(norm, mean, sigma);
      gaus_P_display->SetLineColor(color);
      gaus_P_display->SetLineWidth(2);
      gaus_P_display->Draw("same");
    }

    color_index++;
  }

  // Annotate P canvas with fit results
  TLatex latex_P;
  latex_P.SetNDC();
  latex_P.SetTextSize(0.04);
  latex_P.DrawLatex(0.15, 0.85, Form("Center: %.4f #pm %.4f ns", rf_offset_P, rf_offset_P_err));
  latex_P.DrawLatex(0.15, 0.80, Form("Width: %.4f ns (avg)", sigma_P));

  // Create canvas for H arm modulo histogram
  TCanvas *c_mod_H = new TCanvas(Form("c_mod_H_%d", run_num), "RF Offset Mod H", 800, 600);
  h_mod_H->Draw();

  // Find the peak position for H arm
  max_bin = h_mod_H->GetMaximumBin();
  max_pos = h_mod_H->GetBinCenter(max_bin);

  // Perform multiple fits with different fit ranges
  vector<double> H_means, H_sigmas, H_norms;
  vector<double> H_fit_lows, H_fit_highs;

  for (double fit_range : fit_ranges) {
    // Calculate fit range ±fit_range around peak
    double fit_low_H  = max_pos - fit_range;
    double fit_high_H = max_pos + fit_range;
    // Standard Gaussian fit on extended histogram (no wrapping)
    TF1 *gaus_H_fit = new TF1(Form("gaus_H_fit_%d_%.2f", run_num, fit_range), "gaus", -mod_value, 2 * mod_value);

    // Set initial parameters
    gaus_H_fit->SetParameter(0, h_mod_H->GetMaximum());
    gaus_H_fit->SetParameter(1, max_pos);
    gaus_H_fit->SetParameter(2, 0.2);

    // Fit over the appropriate range(s)
    h_mod_H_ext->Fit(gaus_H_fit, "QN", "", fit_low_H, fit_high_H);

    double fit_mean  = gaus_H_fit->GetParameter(1);
    double fit_sigma = gaus_H_fit->GetParameter(2);

    // Wrap mean if needed
    if (fit_mean < 0) {
      fit_mean += mod_value;
    } else if (fit_mean > mod_value) {
      fit_mean -= mod_value;
    }

    H_means.push_back(fit_mean);
    H_sigmas.push_back(fit_sigma);
    H_norms.push_back(gaus_H_fit->GetParameter(0));
    H_fit_lows.push_back(fit_low_H);
    H_fit_highs.push_back(fit_high_H);
  }

  // Calculate average and std of means
  double rf_offset_H = 0, sigma_H = 0;
  for (double m : H_means) {
    rf_offset_H += m;
  }
  rf_offset_H /= H_means.size();

  // Calculate std of means as uncertainty
  double rf_offset_H_err = 0;
  for (double m : H_means) {
    rf_offset_H_err += (m - rf_offset_H) * (m - rf_offset_H);
  }
  rf_offset_H_err = sqrt(rf_offset_H_err / H_means.size());

  // Average sigma
  for (double s : H_sigmas) {
    sigma_H += s;
  }
  sigma_H /= H_sigmas.size();

  // Set histogram range to ±fit_range_for_display around peak (for display only)
  display_low  = max_pos - fit_range_for_display;
  display_high = max_pos + fit_range_for_display;

  // Handle wrapping for display range
  if (display_low < 0) {
    h_mod_H->GetXaxis()->SetRangeUser(0, display_high);
  } else if (display_high > mod_value) {
    h_mod_H->GetXaxis()->SetRangeUser(display_low, mod_value);
  } else {
    h_mod_H->GetXaxis()->SetRangeUser(display_low, display_high);
  }

  // Draw all fit functions on canvas, only in their fitted ranges
  color_index = 0;
  for (size_t i = 0; i < H_means.size(); i++) {
    double fit_low  = H_fit_lows[i];
    double fit_high = H_fit_highs[i];
    double mean     = H_means[i];
    double sigma    = H_sigmas[i];
    double norm     = H_norms[i];

    int color = colors[color_index % 5];

    if (fit_low < 0) {
      // Draw wrapped segment near upper edge
      TF1 *gaus_H_display1 =
          new TF1(Form("gaus_H_display1_%d_%zu", run_num, i), "gaus", fit_low + mod_value, mod_value);
      gaus_H_display1->SetParameters(norm, mean + mod_value, sigma);
      gaus_H_display1->SetLineColor(color);
      gaus_H_display1->SetLineWidth(2);
      gaus_H_display1->Draw("same");

      TF1 *gaus_H_display2 = new TF1(Form("gaus_H_display2_%d_%zu", run_num, i), "gaus", 0, fit_high);
      gaus_H_display2->SetParameters(norm, mean, sigma);
      gaus_H_display2->SetLineColor(color);
      gaus_H_display2->SetLineWidth(2);
      gaus_H_display2->Draw("same");
    } else if (fit_high > mod_value) {
      // Draw wrapped segment near lower edge
      TF1 *gaus_H_display1 = new TF1(Form("gaus_H_display1_%d_%zu", run_num, i), "gaus", fit_low, mod_value);
      gaus_H_display1->SetParameters(norm, mean, sigma);
      gaus_H_display1->SetLineColor(color);
      gaus_H_display1->SetLineWidth(2);
      gaus_H_display1->Draw("same");

      TF1 *gaus_H_display2 = new TF1(Form("gaus_H_display2_%d_%zu", run_num, i), "gaus", 0, fit_high - mod_value);
      gaus_H_display2->SetParameters(norm, mean - mod_value, sigma);
      gaus_H_display2->SetLineColor(color);
      gaus_H_display2->SetLineWidth(2);
      gaus_H_display2->Draw("same");
    } else {
      // Normal case: draw only in fitted range
      TF1 *gaus_H_display = new TF1(Form("gaus_H_display_%d_%zu", run_num, i), "gaus", fit_low, fit_high);
      gaus_H_display->SetParameters(norm, mean, sigma);
      gaus_H_display->SetLineColor(color);
      gaus_H_display->SetLineWidth(2);
      gaus_H_display->Draw("same");
    }

    color_index++;
  }

  // Annotate H canvas with fit results
  TLatex latex_H;
  latex_H.SetNDC();
  latex_H.SetTextSize(0.04);
  latex_H.DrawLatex(0.15, 0.85, Form("Center: %.4f #pm %.4f ns", rf_offset_H, rf_offset_H_err));
  latex_H.DrawLatex(0.15, 0.80, Form("Width: %.4f ns (avg)", sigma_H));

  // Calculate std of widths for error estimate
  double sigma_P_err = 0;
  for (double s : P_sigmas) {
    sigma_P_err += (s - sigma_P) * (s - sigma_P);
  }
  sigma_P_err = sqrt(sigma_P_err / P_sigmas.size());

  double sigma_H_err = 0;
  for (double s : H_sigmas) {
    sigma_H_err += (s - sigma_H) * (s - sigma_H);
  }
  sigma_H_err = sqrt(sigma_H_err / H_sigmas.size());

  // Store results
  rf_offsets_P.push_back(rf_offset_P);
  rf_offset_errs_P.push_back(rf_offset_P_err);
  rf_widths_P.push_back(sigma_P);
  rf_width_errs_P.push_back(sigma_P_err);

  rf_offsets_H.push_back(rf_offset_H);
  rf_offset_errs_H.push_back(rf_offset_H_err);
  rf_widths_H.push_back(sigma_H);
  rf_width_errs_H.push_back(sigma_H_err);

  // Write to ROOT file
  TDirectory *dir_h_P = out->GetDirectory("h_P");
  if (!dir_h_P) {
    dir_h_P = out->mkdir("h_P");
  }
  dir_h_P->cd();
  h_P->Write();
  out->cd();

  TDirectory *dir_h_mod_P = out->GetDirectory("h_mod_P");
  if (!dir_h_mod_P) {
    dir_h_mod_P = out->mkdir("h_mod_P");
  }
  dir_h_mod_P->cd();
  c_mod_P->Write();
  out->cd();

  TDirectory *dir_h_H = out->GetDirectory("h_H");
  if (!dir_h_H) {
    dir_h_H = out->mkdir("h_H");
  }
  dir_h_H->cd();
  h_H->Write();
  out->cd();

  TDirectory *dir_h_mod_H = out->GetDirectory("h_mod_H");
  if (!dir_h_mod_H) {
    dir_h_mod_H = out->mkdir("h_mod_H");
  }
  dir_h_mod_H->cd();
  c_mod_H->Write();
  out->cd();
}

void get_RF_offset_fit(const char *run_list_file = "run_numbers_all.txt") {

  // Suppress ROOT warnings
  gErrorIgnoreLevel = kError;

  // Filter parameters for stdev plots
  double hminusp_center  = 3.5; // Center value for HminusP range
  double hminusp_width   = 0.2; // +/- width around center
  double uncertainty_max = 0.1; // Maximum uncertainty on data point
  double width_max       = 0.4; // Maximum allowed width

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

  vector<double> rf_offsets_P, rf_offsets_H, rf_offset_errs_P, rf_offset_errs_H;
  vector<double> rf_widths_P, rf_widths_H, rf_width_errs_P, rf_width_errs_H;

  // Open input file with histograms
  TFile *inFile = new TFile("rf_offset_histograms_all.root", "READ");
  TFile *out    = new TFile("rf_offset_results_comp_all.root", "RECREATE");

  int total_runs = run_nums.size();
  int current    = 0;
  for (const auto &run_num : run_nums) {
    current++;
    double percent = (100.0 * current) / total_runs;
    cout << "\rProgress: " << current << "/" << total_runs << " (" << fixed << setprecision(1) << percent << "%) "
         << flush;
    processFitting(run_num, run_nums, rf_offsets_P, rf_offset_errs_P, rf_widths_P, rf_width_errs_P, rf_offsets_H,
                   rf_offset_errs_H, rf_widths_H, rf_width_errs_H, inFile, out);
  }
  cout << endl; // Final newline after progress completes

  inFile->Close();

  // Define vertical line positions and labels
  struct VLine {
    double value;
    const char *label;
  };
  std::vector<VLine> vlines = {{23491, "Beam off"}, {22568, "Beam off"}, {23585, "Beam off"}, 
                               {22229, "Beam off"}, {22359, "Beam off"}, {22489, "Q3 down"}, 
                               {22715, "Beam off"}, {22825, "Beam off"}, {22938, "Beam off"}, 
                               {23248, "Beam off"}, {23380, "Beam off"}, {23430, "Beam off"}, 
                               {22570, "Beam off"}, {23719, "Beam off"}};

  // Helper to draw vertical dashed lines at specific run numbers
  auto draw_vlines = [&](TCanvas *canvas, TH1F *frame, const std::vector<int> &runs) {
    canvas->cd();
    double ymin = frame->GetYaxis()->GetXmin();
    double ymax = frame->GetYaxis()->GetXmax();
    for (const auto &vl : vlines) {
      // Find if the vline value is between two consecutive runs
      int idx_left = -1, idx_right = -1;
      for (size_t i = 0; i < runs.size() - 1; i++) {
        if (vl.value >= runs[i] && vl.value <= runs[i + 1]) {
          idx_left  = i;
          idx_right = i + 1;
          break;
        }
      }
      if (idx_left == -1 || idx_right == -1)
        continue;
      
      // Interpolate x position between the two runs
      double x_pos = runs[idx_left] + (vl.value - runs[idx_left]) * 0.5;
      
      TLine *line = new TLine(x_pos, ymin, x_pos, ymax);
      line->SetLineColor(kBlack);
      line->SetLineStyle(2);
      line->SetLineWidth(1);
      line->Draw();
      
      TLatex latex;
      latex.SetTextAlign(22);
      latex.SetTextSize(0.018);
      latex.SetTextColor(kBlack);
      latex.SetTextAngle(90);
      latex.SetNDC(false);
      latex.DrawLatex(x_pos, ymin + 0.02 * (ymax - ymin), vl.label);
    }
    canvas->Update();
  };

  TCanvas *c_offsets_P      = new TCanvas("c_offsets_P", "RF Offsets P", 800, 600);
  TGraphErrors *g_offsets_P = new TGraphErrors();

  int point_idx = 0;
  for (size_t i = 0; i < rf_offsets_P.size(); i++) {
    if (fabs(rf_offsets_P[i]) > 0.01) {
      g_offsets_P->SetPoint(point_idx, run_nums[i], rf_offsets_P[i]);
      g_offsets_P->SetPointError(point_idx, 0, rf_offset_errs_P[i]);
      point_idx++;
    }
  }

  g_offsets_P->SetMarkerColor(kBlue);
  g_offsets_P->SetLineColor(kBlue);
  g_offsets_P->SetMarkerStyle(20);

  double x_margin = (run_nums.back() - run_nums.front()) * 0.02;
  TH1F *h_frame_offsets_P =
      new TH1F("h_frame_offsets_P", "RF Offsets P Arm", 100, run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_offsets_P->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_P->GetYaxis()->SetTitle("Offset (ns)");
  h_frame_offsets_P->SetStats(0);
  h_frame_offsets_P->Draw();
  g_offsets_P->Draw("P same");
  gPad->Modified();
  gPad->Update();
  double ymin = g_offsets_P->GetYaxis()->GetXmin();
  double ymax = g_offsets_P->GetYaxis()->GetXmax();
  h_frame_offsets_P->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  draw_vlines(c_offsets_P, h_frame_offsets_P, run_nums);
  c_offsets_P->Write();

  // Create plot for H arm RF offsets
  TCanvas *c_offsets_H      = new TCanvas("c_offsets_H", "RF Offsets H", 800, 600);
  TGraphErrors *g_offsets_H = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    if (fabs(rf_offsets_H[i]) > 0.01) {
      g_offsets_H->SetPoint(point_idx, run_nums[i], rf_offsets_H[i]);
      g_offsets_H->SetPointError(point_idx, 0, rf_offset_errs_H[i]);
      point_idx++;
    }
  }

  g_offsets_H->SetMarkerColor(kRed);
  g_offsets_H->SetLineColor(kRed);
  g_offsets_H->SetMarkerStyle(20);

  TH1F *h_frame_offsets_H =
      new TH1F("h_frame_offsets_H", "RF Offsets H Arm", 100, run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_offsets_H->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_H->GetYaxis()->SetTitle("Offset (ns)");
  h_frame_offsets_H->SetStats(0);
  h_frame_offsets_H->Draw();
  g_offsets_H->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_offsets_H->GetYaxis()->GetXmin();
  ymax = g_offsets_H->GetYaxis()->GetXmax();
  h_frame_offsets_H->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  draw_vlines(c_offsets_H, h_frame_offsets_H, run_nums);
  c_offsets_H->Write();

  // Create plot for P arm RF widths
  TCanvas *c_widths_P      = new TCanvas("c_widths_P", "RF Widths P", 800, 600);
  TGraphErrors *g_widths_P = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < rf_widths_P.size(); i++) {
    if (fabs(rf_widths_P[i]) > 0.01) {
      g_widths_P->SetPoint(point_idx, run_nums[i], rf_widths_P[i]);
      g_widths_P->SetPointError(point_idx, 0, rf_width_errs_P[i]);
      point_idx++;
    }
  }

  g_widths_P->SetMarkerColor(kBlue);
  g_widths_P->SetLineColor(kBlue);
  g_widths_P->SetMarkerStyle(20);

  TH1F *h_frame_widths_P =
      new TH1F("h_frame_widths_P", "RF Widths P Arm", 100, run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_widths_P->GetXaxis()->SetTitle("Run Number");
  h_frame_widths_P->GetYaxis()->SetTitle("Width (ns)");
  h_frame_widths_P->SetStats(0);
  h_frame_widths_P->Draw();
  g_widths_P->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_widths_P->GetYaxis()->GetXmin();
  ymax = g_widths_P->GetYaxis()->GetXmax();
  h_frame_widths_P->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  c_widths_P->Write();

  // Create plot for H arm RF widths
  TCanvas *c_widths_H      = new TCanvas("c_widths_H", "RF Widths H", 800, 600);
  TGraphErrors *g_widths_H = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < rf_widths_H.size(); i++) {
    if (fabs(rf_widths_H[i]) > 0.01) {
      g_widths_H->SetPoint(point_idx, run_nums[i], rf_widths_H[i]);
      g_widths_H->SetPointError(point_idx, 0, rf_width_errs_H[i]);
      point_idx++;
    }
  }

  g_widths_H->SetMarkerColor(kRed);
  g_widths_H->SetLineColor(kRed);
  g_widths_H->SetMarkerStyle(20);

  TH1F *h_frame_widths_H =
      new TH1F("h_frame_widths_H", "RF Widths H Arm", 100, run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_widths_H->GetXaxis()->SetTitle("Run Number");
  h_frame_widths_H->GetYaxis()->SetTitle("Width (ns)");
  h_frame_widths_H->SetStats(0);
  h_frame_widths_H->Draw();
  g_widths_H->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_widths_H->GetYaxis()->GetXmin();
  ymax = g_widths_H->GetYaxis()->GetXmax();
  h_frame_widths_H->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  c_widths_H->Write();

  // Create plot for HMS - SHMS offset (mod mod_value)
  TCanvas *c_offsets_diff      = new TCanvas("c_offsets_HminusP", "HMS - SHMS Offsets", 800, 600);
  TGraphErrors *g_offsets_diff = new TGraphErrors();
  std::vector<double> diffs;
  diffs.reserve(rf_offsets_H.size());

  point_idx = 0;
  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    double diff = rf_offsets_H[i] - rf_offsets_P[i];
    while (diff < 0) {
      diff += mod_value;
    }
    while (diff >= mod_value) {
      diff -= mod_value;
    }
    diffs.push_back(diff);
    if (fabs(rf_offsets_H[i]) > 0.01 && fabs(rf_offsets_P[i]) > 0.01) {
      double err = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + rf_offset_errs_P[i] * rf_offset_errs_P[i]);
      g_offsets_diff->SetPoint(point_idx, run_nums[i], diff);
      g_offsets_diff->SetPointError(point_idx, 0, err);
      point_idx++;
    }
  }

  g_offsets_diff->SetMarkerColor(kBlack);
  g_offsets_diff->SetLineColor(kBlack);
  g_offsets_diff->SetMarkerStyle(20);

  TH1F *h_frame_offsets_diff = new TH1F("h_frame_offsets_diff", "HMS - SHMS (mod N)", 100, run_nums.front() - x_margin,
                                        run_nums.back() + x_margin);
  h_frame_offsets_diff->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_diff->GetYaxis()->SetTitle("HMS - SHMS Offset (ns)");
  h_frame_offsets_diff->SetStats(0);
  h_frame_offsets_diff->Draw();
  g_offsets_diff->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_offsets_diff->GetYaxis()->GetXmin();
  ymax = g_offsets_diff->GetYaxis()->GetXmax();
  h_frame_offsets_diff->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  c_offsets_diff->Write();

  // Compute stdev of HMS - SHMS offsets for error inflation, using only filtered values
  std::vector<double> filtered_diffs;
  for (size_t i = 0; i < diffs.size(); i++) {
    double diff = diffs[i];
    bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
    bool width_ok_P = rf_widths_P[i] < width_max;
    bool width_ok_H = rf_widths_H[i] < width_max;
    
    if (within_hminusp_range && width_ok_P && width_ok_H) {
      filtered_diffs.push_back(diff);
    }
  }
  
  double diff_stdev = 0;
  if (!filtered_diffs.empty()) {
    double diff_mean = 0;
    for (double d : filtered_diffs) {
      diff_mean += d;
    }
    diff_mean /= filtered_diffs.size();

    double diff_var = 0;
    for (double d : filtered_diffs) {
      diff_var += (d - diff_mean) * (d - diff_mean);
    }
    diff_stdev = sqrt(diff_var / filtered_diffs.size());
  }
  diff_stdev = 0; // Set to 0 for now, can be adjusted later if needed

  // Create plot for HMS - SHMS offset (filtered, passing cuts)
  TCanvas *c_offsets_diff_filtered      = new TCanvas("c_offsets_HminusP_filtered", "HMS - SHMS Offsets (Filtered)", 800, 600);
  TGraphErrors *g_offsets_diff_filtered = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < diffs.size(); i++) {
    double diff = diffs[i];
    bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
    bool width_ok_P = rf_widths_P[i] < width_max;
    bool width_ok_H = rf_widths_H[i] < width_max;
    
    if (within_hminusp_range && width_ok_P && width_ok_H) {
      if (fabs(rf_offsets_H[i]) > 0.01 && fabs(rf_offsets_P[i]) > 0.01) {
        double err = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + rf_offset_errs_P[i] * rf_offset_errs_P[i]);
        g_offsets_diff_filtered->SetPoint(point_idx, run_nums[i], diff);
        g_offsets_diff_filtered->SetPointError(point_idx, 0, err);
        point_idx++;
      }
    }
  }

  g_offsets_diff_filtered->SetMarkerColor(kBlack);
  g_offsets_diff_filtered->SetLineColor(kBlack);
  g_offsets_diff_filtered->SetMarkerStyle(20);

  TH1F *h_frame_offsets_diff_filtered = new TH1F("h_frame_offsets_diff_filtered", "HMS - SHMS (Filtered)", 100, 
                                                  run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_offsets_diff_filtered->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_diff_filtered->GetYaxis()->SetTitle("HMS - SHMS Offset (ns)");
  h_frame_offsets_diff_filtered->SetStats(0);
  h_frame_offsets_diff_filtered->Draw();
  g_offsets_diff_filtered->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_offsets_diff_filtered->GetYaxis()->GetXmin();
  ymax = g_offsets_diff_filtered->GetYaxis()->GetXmax();
  h_frame_offsets_diff_filtered->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  c_offsets_diff_filtered->Write();

  // Create plot for HMS - SHMS offset (mod mod_value)
  TCanvas *c_offsets_P_stdev      = new TCanvas("c_offsets_P_stdev", "RF Offsets P (stdev)", 800, 600);
  TGraphErrors *g_offsets_P_stdev = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < rf_offsets_P.size(); i++) {
    if (fabs(rf_offsets_P[i]) > 0.01) {
      // Check filter conditions
      double diff = diffs[i];
      double err  = sqrt(rf_offset_errs_P[i] * rf_offset_errs_P[i] + diff_stdev * diff_stdev);
      bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
      bool uncertainty_ok       = err < uncertainty_max;
      bool width_ok             = rf_widths_P[i] < width_max;

      if (within_hminusp_range && uncertainty_ok && width_ok) {
        g_offsets_P_stdev->SetPoint(point_idx, run_nums[i], rf_offsets_P[i]);
        g_offsets_P_stdev->SetPointError(point_idx, 0, err);
        point_idx++;
      }
    }
  }

  g_offsets_P_stdev->SetMarkerColor(kBlue);
  g_offsets_P_stdev->SetLineColor(kBlue);
  g_offsets_P_stdev->SetMarkerStyle(20);

  TH1F *h_frame_offsets_P_stdev = new TH1F("h_frame_offsets_P_stdev", "RF Offsets P Arm (stdev)", 100,
                                           run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_offsets_P_stdev->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_P_stdev->GetYaxis()->SetTitle("Offset (ns)");
  h_frame_offsets_P_stdev->SetStats(0);
  h_frame_offsets_P_stdev->Draw();
  g_offsets_P_stdev->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_offsets_P_stdev->GetYaxis()->GetXmin();
  ymax = g_offsets_P_stdev->GetYaxis()->GetXmax();
  h_frame_offsets_P_stdev->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  draw_vlines(c_offsets_P_stdev, h_frame_offsets_P_stdev, run_nums);
  c_offsets_P_stdev->Write();

  // Create plot for H arm RF offsets with inflated errors
  TCanvas *c_offsets_H_stdev      = new TCanvas("c_offsets_H_stdev", "RF Offsets H (stdev)", 800, 600);
  TGraphErrors *g_offsets_H_stdev = new TGraphErrors();

  point_idx = 0;
  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    if (fabs(rf_offsets_H[i]) > 0.01) {
      // Check filter conditions
      double diff               = diffs[i];
      double err                = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + diff_stdev * diff_stdev);
      bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
      bool uncertainty_ok       = err < uncertainty_max;
      bool width_ok             = rf_widths_H[i] < width_max;

      if (within_hminusp_range && uncertainty_ok && width_ok) {
        g_offsets_H_stdev->SetPoint(point_idx, run_nums[i], rf_offsets_H[i]);
        g_offsets_H_stdev->SetPointError(point_idx, 0, err);
        point_idx++;
      }
    }
  }

  g_offsets_H_stdev->SetMarkerColor(kRed);
  g_offsets_H_stdev->SetLineColor(kRed);
  g_offsets_H_stdev->SetMarkerStyle(20);

  TH1F *h_frame_offsets_H_stdev = new TH1F("h_frame_offsets_H_stdev", "RF Offsets H Arm (stdev)", 100,
                                           run_nums.front() - x_margin, run_nums.back() + x_margin);
  h_frame_offsets_H_stdev->GetXaxis()->SetTitle("Run Number");
  h_frame_offsets_H_stdev->GetYaxis()->SetTitle("Offset (ns)");
  h_frame_offsets_H_stdev->SetStats(0);
  h_frame_offsets_H_stdev->Draw();
  g_offsets_H_stdev->Draw("P same");
  gPad->Modified();
  gPad->Update();
  ymin = g_offsets_H_stdev->GetYaxis()->GetXmin();
  ymax = g_offsets_H_stdev->GetYaxis()->GetXmax();
  h_frame_offsets_H_stdev->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
  draw_vlines(c_offsets_H_stdev, h_frame_offsets_H_stdev, run_nums);
  c_offsets_H_stdev->Write();

  // Compute statistics for specified run ranges
  struct RunRangeStats {
    int run_low, run_high;
    int count_P = 0, count_H = 0;
    double avg_P = 0, avg_H = 0;
    double stdev_P = 0, stdev_H = 0;
  };

  std::vector<RunRangeStats> ranges = {
    {22178, 22480},
    {22560, 22660},
    {22661, 22711},
    {22745, 23590},
    {23597, 23795}
  };

  for (auto &range : ranges) {
    std::vector<double> values_P, values_H;

    // Collect values for P arm
    for (size_t i = 0; i < rf_offsets_P.size(); i++) {
      if (run_nums[i] >= range.run_low && run_nums[i] <= range.run_high) {
        if (fabs(rf_offsets_P[i]) > 0.01) {
          double diff = diffs[i];
          double err_P = sqrt(rf_offset_errs_P[i] * rf_offset_errs_P[i] + diff_stdev * diff_stdev);
          bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
          bool uncertainty_ok = err_P < uncertainty_max;
          bool width_ok = rf_widths_P[i] < width_max;
          
          if (within_hminusp_range && uncertainty_ok && width_ok) {
            values_P.push_back(rf_offsets_P[i]);
            range.count_P++;
          }
        }
      }
    }

    // Collect values for H arm
    for (size_t i = 0; i < rf_offsets_H.size(); i++) {
      if (run_nums[i] >= range.run_low && run_nums[i] <= range.run_high) {
        if (fabs(rf_offsets_H[i]) > 0.01) {
          double diff = diffs[i];
          double err_H = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + diff_stdev * diff_stdev);
          bool within_hminusp_range = (diff >= hminusp_center - hminusp_width) && (diff <= hminusp_center + hminusp_width);
          bool uncertainty_ok = err_H < uncertainty_max;
          bool width_ok = rf_widths_H[i] < width_max;
          
          if (within_hminusp_range && uncertainty_ok && width_ok) {
            values_H.push_back(rf_offsets_H[i]);
            range.count_H++;
          }
        }
      }
    }

    // Calculate averages
    if (!values_P.empty()) {
      for (double v : values_P) {
        range.avg_P += v;
      }
      range.avg_P /= values_P.size();
      
      for (double v : values_P) {
        range.stdev_P += (v - range.avg_P) * (v - range.avg_P);
      }
      range.stdev_P = sqrt(range.stdev_P / values_P.size());
    }

    if (!values_H.empty()) {
      for (double v : values_H) {
        range.avg_H += v;
      }
      range.avg_H /= values_H.size();
      
      for (double v : values_H) {
        range.stdev_H += (v - range.avg_H) * (v - range.avg_H);
      }
      range.stdev_H = sqrt(range.stdev_H / values_H.size());
    }
  }

  // Create canvas with statistics table
  TCanvas *c_stats = new TCanvas("c_range_statistics", "Range Statistics", 1000, 600);
  c_stats->cd();

  double y_pos = 0.95;
  double line_height = 0.08;
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.018);
  latex.SetTextFont(62);
  
  // Header
  latex.DrawLatex(0.05, y_pos, "Run Range");
  latex.DrawLatex(0.18, y_pos, "N(P)");
  latex.DrawLatex(0.28, y_pos, "N(H)");
  latex.DrawLatex(0.38, y_pos, "Avg P (ns)");
  latex.DrawLatex(0.55, y_pos, "Avg H (ns)");
  latex.DrawLatex(0.72, y_pos, "Stdev P (ns)");
  latex.DrawLatex(0.88, y_pos, "Stdev H (ns)");
  
  latex.SetTextFont(42);
  y_pos -= line_height;
  
  for (const auto &range : ranges) {
    latex.DrawLatex(0.05, y_pos, Form("%d-%d", range.run_low, range.run_high));
    latex.DrawLatex(0.18, y_pos, Form("%d", range.count_P));
    latex.DrawLatex(0.28, y_pos, Form("%d", range.count_H));
    latex.DrawLatex(0.38, y_pos, Form("%.4f", range.avg_P));
    latex.DrawLatex(0.55, y_pos, Form("%.4f", range.avg_H));
    latex.DrawLatex(0.72, y_pos, Form("%.4f", range.stdev_P));
    latex.DrawLatex(0.88, y_pos, Form("%.4f", range.stdev_H));
    y_pos -= line_height;
  }
  
  c_stats->Write();

  out->Close();
}
