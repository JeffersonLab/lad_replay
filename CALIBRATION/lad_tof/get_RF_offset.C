#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <cmath>
#include <vector>

using namespace std;
double mod_value = 4.00801; // RF period in ns

void processRun(int run_num, vector<double> &rf_offsets_P, vector<double> &rf_offset_errs_P,
                vector<double> &rf_widths_P, vector<double> &rf_width_errs_P, vector<double> &rf_offsets_H,
                vector<double> &rf_offset_errs_H, vector<double> &rf_widths_H, vector<double> &rf_width_errs_H,
                TFile *out) {
  gROOT->SetBatch(kTRUE);
  TFile *f = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                            "LAD_COIN_%d_0_0_50000.root",
                            run_num));
  TTree *T = (TTree *)f->Get("T");

  vector<double> fit_ranges = {0.2, 0.25, 0.3}; // Multiple fit ranges to test

  // Helper lambda to get center 98th percentile range
  auto getPercentileRange = [](TH1F *h, double &low, double &high) {
    double sum    = h->Integral();
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

  // P arm - first pass to get raw histogram for percentile calculation
  TH1F *h_P_temp = new TH1F(Form("h_P_temp_%d", run_num), "RF Offset P Temp", 1000, 0, 4000);
  T->Draw(Form("P.ladkin.t_vertex-T.shms.pRF_tdcTime>>h_P_temp_%d", run_num), "g.evtyp==1");

  double p_low, p_high;
  getPercentileRange(h_P_temp, p_low, p_high);

  // Create final P histogram with 98th percentile range
  TH1F *h_P = new TH1F(Form("h_P_%d", run_num), "RF Offset P", 500, p_low, p_high);
  T->Draw(Form("P.ladkin.t_vertex-T.shms.pRF_tdcTime>>h_P_%d", run_num), "g.evtyp==1");

  // H arm - first pass to get raw histogram for percentile calculation
  TH1F *h_H_temp = new TH1F(Form("h_H_temp_%d", run_num), "RF Offset H Temp", 1000, 0, 4000);
  T->Draw(Form("H.ladkin.t_vertex-T.hms.hRF_tdcTime>>h_H_temp_%d", run_num), "g.evtyp==2");

  double h_low, h_high;
  getPercentileRange(h_H_temp, h_low, h_high);

  // Create final H histogram with 98th percentile range
  TH1F *h_H = new TH1F(Form("h_H_%d", run_num), "RF Offset H", 500, h_low, h_high);
  T->Draw(Form("H.ladkin.t_vertex-T.hms.hRF_tdcTime>>h_H_%d", run_num), "g.evtyp==2");
  TH1F *h_mod_P = new TH1F(Form("h_mod_P_%d", run_num), "RF Offset Mod P", 500, 0, mod_value);
  T->Draw(Form("fmod(P.ladkin.t_vertex-T.shms.pRF_tdcTime,%.5f)>>h_mod_P_%d", mod_value, run_num), "g.evtyp==1");

  TH1F *h_mod_H = new TH1F(Form("h_mod_H_%d", run_num), "RF Offset Mod H", 500, 0, mod_value);
  T->Draw(Form("fmod(H.ladkin.t_vertex-T.hms.hRF_tdcTime,%.5f)>>h_mod_H_%d", mod_value, run_num), "g.evtyp==2");

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
  vector<double> P_fit_lows, P_fit_highs;                           // Store fit ranges
  double fit_range_for_display = fit_ranges[fit_ranges.size() / 2]; // Use middle fit range for display

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

  // Annotate P canvas with fit results (center and width with uncertainties)
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
  vector<double> H_fit_lows, H_fit_highs; // Store fit ranges

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

  // Annotate H canvas with fit results (center and width with uncertainties)
  TLatex latex_H;
  latex_H.SetNDC();
  latex_H.SetTextSize(0.04);
  latex_H.DrawLatex(0.15, 0.85, Form("Center: %.4f #pm %.4f ns", rf_offset_H, rf_offset_H_err));
  latex_H.DrawLatex(0.15, 0.80, Form("Width: %.4f ns (avg)", sigma_H));

  std::cout << "RF Offset P: " << rf_offset_P << " +/- " << rf_offset_P_err << " ns" << std::endl;
  std::cout << "RF Offset H: " << rf_offset_H << " +/- " << rf_offset_H_err << " ns" << std::endl;

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

  // Write to ROOT file in separate directories
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

void get_RF_offset() {

  vector<int> run_nums = {22290, 22353, 22449, 22602, 22760, 22919, 23237, 23414, 23448, 23471, 23598, 23770,
                          22291, 22354, 22450, 22603, 22761, 22920, 23238, 23415, 23449, 23472, 23600, 23771};
  // removed 23018, since no SHMS
  std::sort(run_nums.begin(), run_nums.end());

  vector<double> rf_offsets_P, rf_offsets_H, rf_offset_errs_P, rf_offset_errs_H;
  vector<double> rf_widths_P, rf_widths_H, rf_width_errs_P, rf_width_errs_H;

  TFile *out = new TFile("rf_offset_results_comp.root", "RECREATE");

  for (const auto &run_num : run_nums) {
    cout << "Processing run " << run_num << "..." << endl;

    // Call the core function to process each run
    processRun(run_num, rf_offsets_P, rf_offset_errs_P, rf_widths_P, rf_width_errs_P, rf_offsets_H, rf_offset_errs_H,
               rf_widths_H, rf_width_errs_H, out);
  }

  // Define vertical line positions and labels
  struct VLine {
    double value;
    const char *label;
  };
  std::vector<VLine> vlines = {
      {22500, "Switchboard failure"},
      {23013.5, "HMS rotation"},
      {23461.5, "SHMS rotation"}};

  // Helper to draw vertical dashed lines between closest two points
  auto draw_vlines = [&](TCanvas *canvas, TH1F *frame, const std::vector<int> &runs) {
    canvas->cd();
    double ymin = frame->GetYaxis()->GetXmin();
    double ymax = frame->GetYaxis()->GetXmax();
    for (const auto &vl : vlines) {
      // Find closest two run indices
      int idx_left = -1, idx_right = -1;
      for (size_t i = 0; i < runs.size() - 1; i++) {
        if (vl.value > runs[i] && vl.value < runs[i + 1]) {
          idx_left = i;
          idx_right = i + 1;
          break;
        }
      }
      if (idx_left == -1 || idx_right == -1)
        continue;
      double xline = idx_left + 0.5; // halfway between
      TLine *line = new TLine(xline, ymin, xline, ymax);
      line->SetLineColor(kBlack);
      line->SetLineStyle(2); // dashed
      line->SetLineWidth(2);
      line->Draw();
      TLatex latex;
      latex.SetTextAlign(22);
      latex.SetTextSize(0.0175);
      latex.SetTextColor(kBlack);
      latex.SetNDC(false);
      latex.DrawLatex(xline, 0.5 * (ymin + ymax), vl.label);
    }
    canvas->Update();
  };
  TCanvas *c_offsets_P      = new TCanvas("c_offsets_P", "RF Offsets P", 800, 600);
  TGraphErrors *g_offsets_P = new TGraphErrors(rf_offsets_P.size());

  for (size_t i = 0; i < rf_offsets_P.size(); i++) {
    g_offsets_P->SetPoint(i, i, rf_offsets_P[i]);
    g_offsets_P->SetPointError(i, 0, rf_offset_errs_P[i]);
  }

  g_offsets_P->SetMarkerColor(kBlue);
  g_offsets_P->SetLineColor(kBlue);
  g_offsets_P->SetMarkerStyle(20);

  TH1F *h_frame_offsets_P =
      new TH1F("h_frame_offsets_P", "RF Offsets P Arm", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_offsets_P->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_offsets_P->GetXaxis()->SetLabelSize(0.04);
  h_frame_offsets_P->GetXaxis()->LabelsOption("v");
  h_frame_offsets_P->GetXaxis()->SetLabelOffset(0.01);
  h_frame_offsets_P->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);

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
  TGraphErrors *g_offsets_H = new TGraphErrors(rf_offsets_H.size());

  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    g_offsets_H->SetPoint(i, i, rf_offsets_H[i]);
    g_offsets_H->SetPointError(i, 0, rf_offset_errs_H[i]);
  }

  g_offsets_H->SetMarkerColor(kRed);
  g_offsets_H->SetLineColor(kRed);
  g_offsets_H->SetMarkerStyle(20);

  TH1F *h_frame_offsets_H =
      new TH1F("h_frame_offsets_H", "RF Offsets H Arm", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_offsets_H->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_offsets_H->GetXaxis()->SetLabelSize(0.04);
  h_frame_offsets_H->GetXaxis()->LabelsOption("v");
  h_frame_offsets_H->GetXaxis()->SetLabelOffset(0.01);
  h_frame_offsets_H->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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
  TGraphErrors *g_widths_P = new TGraphErrors(rf_widths_P.size());

  for (size_t i = 0; i < rf_widths_P.size(); i++) {
    g_widths_P->SetPoint(i, i, rf_widths_P[i]);
    g_widths_P->SetPointError(i, 0, rf_width_errs_P[i]);
  }

  g_widths_P->SetMarkerColor(kBlue);
  g_widths_P->SetLineColor(kBlue);
  g_widths_P->SetMarkerStyle(20);

  TH1F *h_frame_widths_P =
      new TH1F("h_frame_widths_P", "RF Widths P Arm", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_widths_P->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_widths_P->GetXaxis()->SetLabelSize(0.04);
  h_frame_widths_P->GetXaxis()->LabelsOption("v");
  h_frame_widths_P->GetXaxis()->SetLabelOffset(0.01);
  h_frame_widths_P->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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
  TGraphErrors *g_widths_H = new TGraphErrors(rf_widths_H.size());

  for (size_t i = 0; i < rf_widths_H.size(); i++) {
    g_widths_H->SetPoint(i, i, rf_widths_H[i]);
    g_widths_H->SetPointError(i, 0, rf_width_errs_H[i]);
  }

  g_widths_H->SetMarkerColor(kRed);
  g_widths_H->SetLineColor(kRed);
  g_widths_H->SetMarkerStyle(20);

  TH1F *h_frame_widths_H =
      new TH1F("h_frame_widths_H", "RF Widths H Arm", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_widths_H->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_widths_H->GetXaxis()->SetLabelSize(0.04);
  h_frame_widths_H->GetXaxis()->LabelsOption("v");
  h_frame_widths_H->GetXaxis()->SetLabelOffset(0.01);
  h_frame_widths_H->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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
  TGraphErrors *g_offsets_diff = new TGraphErrors(rf_offsets_H.size());
  std::vector<double> diffs;
  diffs.reserve(rf_offsets_H.size());

  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    double diff = rf_offsets_H[i] - rf_offsets_P[i];
    while (diff < 0) {
      diff += mod_value;
    }
    while (diff >= mod_value) {
      diff -= mod_value;
    }
    diffs.push_back(diff);
    double err = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + rf_offset_errs_P[i] * rf_offset_errs_P[i]);
    g_offsets_diff->SetPoint(i, i, diff);
    g_offsets_diff->SetPointError(i, 0, err);
  }

  g_offsets_diff->SetMarkerColor(kBlack);
  g_offsets_diff->SetLineColor(kBlack);
  g_offsets_diff->SetMarkerStyle(20);

  TH1F *h_frame_offsets_diff =
      new TH1F("h_frame_offsets_diff", "HMS - SHMS (mod N)", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_offsets_diff->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_offsets_diff->GetXaxis()->SetLabelSize(0.04);
  h_frame_offsets_diff->GetXaxis()->LabelsOption("v");
  h_frame_offsets_diff->GetXaxis()->SetLabelOffset(0.01);
  h_frame_offsets_diff->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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

  // Compute stdev of HMS - SHMS offsets for error inflation
  double diff_mean = 0;
  for (double d : diffs) {
    diff_mean += d;
  }
  diff_mean /= diffs.size();

  double diff_var = 0;
  for (double d : diffs) {
    diff_var += (d - diff_mean) * (d - diff_mean);
  }
  double diff_stdev = sqrt(diff_var / diffs.size());

  // Create plot for P arm RF offsets with inflated errors
  TCanvas *c_offsets_P_stdev      = new TCanvas("c_offsets_P_stdev", "RF Offsets P (stdev)", 800, 600);
  TGraphErrors *g_offsets_P_stdev = new TGraphErrors(rf_offsets_P.size());

  for (size_t i = 0; i < rf_offsets_P.size(); i++) {
    double err = sqrt(rf_offset_errs_P[i] * rf_offset_errs_P[i] + diff_stdev * diff_stdev);
    g_offsets_P_stdev->SetPoint(i, i, rf_offsets_P[i]);
    g_offsets_P_stdev->SetPointError(i, 0, err);
  }

  g_offsets_P_stdev->SetMarkerColor(kBlue);
  g_offsets_P_stdev->SetLineColor(kBlue);
  g_offsets_P_stdev->SetMarkerStyle(20);

  TH1F *h_frame_offsets_P_stdev =
      new TH1F("h_frame_offsets_P_stdev", "RF Offsets P Arm (stdev)", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_offsets_P_stdev->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_offsets_P_stdev->GetXaxis()->SetLabelSize(0.04);
  h_frame_offsets_P_stdev->GetXaxis()->LabelsOption("v");
  h_frame_offsets_P_stdev->GetXaxis()->SetLabelOffset(0.01);
  h_frame_offsets_P_stdev->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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
  TGraphErrors *g_offsets_H_stdev = new TGraphErrors(rf_offsets_H.size());

  for (size_t i = 0; i < rf_offsets_H.size(); i++) {
    double err = sqrt(rf_offset_errs_H[i] * rf_offset_errs_H[i] + diff_stdev * diff_stdev);
    g_offsets_H_stdev->SetPoint(i, i, rf_offsets_H[i]);
    g_offsets_H_stdev->SetPointError(i, 0, err);
  }

  g_offsets_H_stdev->SetMarkerColor(kRed);
  g_offsets_H_stdev->SetLineColor(kRed);
  g_offsets_H_stdev->SetMarkerStyle(20);

  TH1F *h_frame_offsets_H_stdev =
      new TH1F("h_frame_offsets_H_stdev", "RF Offsets H Arm (stdev)", run_nums.size(), -0.5, run_nums.size() - 0.5);
  for (size_t i = 0; i < run_nums.size(); i++) {
    h_frame_offsets_H_stdev->GetXaxis()->SetBinLabel(i + 1, Form("%d", run_nums[i]));
  }
  h_frame_offsets_H_stdev->GetXaxis()->SetLabelSize(0.04);
  h_frame_offsets_H_stdev->GetXaxis()->LabelsOption("v");
  h_frame_offsets_H_stdev->GetXaxis()->SetLabelOffset(0.01);
  h_frame_offsets_H_stdev->GetXaxis()->SetNdivisions(run_nums.size(), 0, 0, kTRUE);
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

  out->Close();
}
