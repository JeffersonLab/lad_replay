#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TObject.h>
#include <TPad.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

const char spec_prefix = 'P'; // Spectrometer to replay

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

const double fit_func_xmin_back  = 2.5; // 2.2; // Minimum x value for fit function
const double fit_func_xmax_back  = 7.0; // Maximum x value for fit function
const double fit_func_xmin_front = 2.5; // Minimum x value for fit function
const double fit_func_xmax_front = 7.0; // Maximum x value for fit function
const double fit_func_mid        = 3.2; // Midpoint for piecewise functions

// const double quad_fit_param[2] = {164, 42};
// const double quad_fit_param[2] = {5.31637, 12.9101};
const double quad_fit_param[2] = {0, 13};
double cut_line_pt1[2]         = {2, 120};  // Point 1 for cut line
const double cut_line_pt2[2]   = {6, 0};    // Point 2 for cut line
const double int_amp_ratio     = 1.0 / 5.0; // Ratio of integral to amplitude ADC values

const double true_adc_MeV = 80;   // True ADC value in MeV
const double line_pt_tol  = 50.0; // adjust tolerance as needed

const int rebinX = 2; // set >1 to rebin in X
const int rebinY = 2; // set >1 to rebin in Y

const double front_box_cut[2] = {4, 70}; // points >x = [0] and <y = [1] are killed

double max_edep[2][N_PLANES][N_PADDLES][2] = {0}; // [amp/int][x/y][plane][paddle]

double edep_fit_func_profile(double *x, double *par) {
  // x[0]: x, x[1]: y
  // par[0]: amplitude, par[1]: x-mean, par[2]: x-sigma, par[3]: y-mean, par[4]: y-sigma, par[5]: x0, par[6]: constant
  // 2D Gaussian for x < x0, else constant
  // Only one free parameter: par[0] is the y-axis scaling factor
  double y;
  if (x[0] < fit_func_mid) {
    // Line between (1.5,40) and (fit_func_mid,180)
    double x1 = 1.5, y1 = 40.0;
    double x2 = fit_func_mid, y2 = 180.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  } else {
    // Line between (fit_func_mid,180) and (6,100)
    double x1 = fit_func_mid, y1 = 180.0;
    double x2 = 6.0, y2 = 100.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  }
  return par[0] * y;
}

// double edep_fit_func_profile_quad_back(double *x, double *par) {
//   double dx = par[1];
//   double y;
//   // Fit two quadratics: one for x < fit_func_mid, one for x >= fit_func_mid
//   // par[0]: scaling
//   // par[1]: a_left, par[2]: b_left, par[3]: c_left
//   // par[4]: a_right, par[5]: b_right, par[6]: c_right
//   double x1, y1, x2, y2, a;
//   if (x[0] - dx < fit_func_mid) {
//     // Quadratic passing through (1,20) and (3.3,190)
//     // One free parameter
//     x1 = 2, y1 = 20.0;
//     // x1 = par[2], y1 = par[3];
//     x2 = fit_func_mid, y2 = 190.0;
//     a = quad_fit_param[0];
//   } else {
//     // Quadratic passing through (3.3,190) and (7,30)
//     x1 = fit_func_mid, y1 = 190.0;
//     x2 = 6.0, y2 = 30.0;
//     a = quad_fit_param[1];
//   }
//   // double a = par[0];
//   double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
//   y            = a * ((x[0] - dx) * (x[0] - dx) - x1 * x1) + b_tmp * ((x[0] - dx) - x1) + y1;

//   return y * par[0]; // Scale by par[0]
// }

double edep_fit_func_profile_quad_back(double *x, double *par) {
  double dx = par[1];
  double y;
  // Fit two quadratics: one for x < fit_func_mid, one for x >= fit_func_mid
  // par[0]: scaling
  // par[1]: a_left, par[2]: b_left, par[3]: c_left
  // par[4]: a_right, par[5]: b_right, par[6]: c_right
  double x1, y1, x2, y2, a;
  double fit_func_midtmp = fit_func_mid - dx;
  if (x[0] < fit_func_midtmp) {
    // Quadratic passing through (1,20) and (3.3,190)
    // One free parameter
    // x1 = 2, y1 = 20.0;
    // x1 = 2, y1 = 20.0;
    x1 = par[2], y1 = 20.0;
    x2 = fit_func_midtmp, y2 = 190.0;
    a = quad_fit_param[0];
  } else {
    // Quadratic passing through (3.3,190) and (7,50)
    x1 = fit_func_midtmp, y1 = 190.0;
    x2 = 6.0, y2 = 44.0;
    a = quad_fit_param[1];
  }
  // double a = par[0];
  dx           = 0;
  double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
  y            = a * ((x[0] - dx) * (x[0] - dx) - x1 * x1) + b_tmp * ((x[0] - dx) - x1) + y1;

  return y * par[0]; // Scale by par[0]
}

// double edep_fit_func_profile_quad_back(double *x, double *par) {
//   double y;
//   // Fit two quadratics: one for x < fit_func_mid, one for x >= fit_func_mid
//   // par[0]: scaling
//   // par[1]: a_left, par[2]: b_left, par[3]: c_left
//   // par[4]: a_right, par[5]: b_right, par[6]: c_right
//   double x1, y1, x2, y2, a;
//   if (x[0] < fit_func_mid) {
//     // Quadratic passing through (1,20) and (3.3,190)
//     // One free parameter 
//     x1 = 0.8, y1 = 116;
//     x2 = fit_func_mid, y2 = 210.0;
//     // a = par[1];
//     a = quad_fit_param[0];
//   } else {
//     // Quadratic passing through (3.3,190) and (7,20)
//     x1 = fit_func_mid, y1 = 210;
//     x2 = 3.0, y2 = 107;
//     // a = par[4];
//     a = quad_fit_param[1];
//   }
//   // double a = par[0];
//   double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
//   y            = a * (x[0] * x[0] - x1 * x1) + b_tmp * (x[0] - x1) + y1;

//   return y * par[0]; // Scale by par[0]
// }

double edep_fit_func_profile_quad_front(double *x, double *par) {
  // double dx = par[1];
  double dx = 0;
  double y;
  // Fit two quadratics: one for x < fit_func_mid, one for x >= fit_func_mid
  // par[0]: scaling
  // par[1]: a_left, par[2]: b_left, par[3]: c_left
  // par[4]: a_right, par[5]: b_right, par[6]: c_right

  // Quadratic passing through (fit_func_mid,180) and (6,20)
  // double x1 = 1 - dx, y1 = 80.0;
  // double x2 = 3 - dx, y2 = 170.0;
  // double a = quad_fit_param[1];
  double a  = par[0];
  double y1 = par[1];
  double x0 = par[2];
  double y0 = par[3];
  // double a = par[1];
  // double a = par[0];
  double x1 = fit_func_xmax_front;
  double b  = ((y0 - y1) - a * (x0 * x0 - x1 * x1)) / (x0 - x1);
  y         = y0 + a * (x[0] * x[0] - x0 * x0) + b * (x[0] - x0);
  // double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
  // y            = a * ((x[0] - dx) * (x[0] - dx) - x1 * x1) + b_tmp * ((x[0] - dx) - x1) + y1;

  return y;
  // return y * par[0]; // Scale by par[0]
}

template <typename HistType> TGraph *get_peak_graph(HistType *hist_orig, bool is_int = false, int y_cut = 0) {

  double m, b;
  if (is_int) {
    m = (cut_line_pt2[1] - cut_line_pt1[1]) * int_amp_ratio / (cut_line_pt2[0] - cut_line_pt1[0]);
    b = cut_line_pt1[1] * int_amp_ratio - m * cut_line_pt1[0];
  } else {
    m = (cut_line_pt2[1] - cut_line_pt1[1]) / (cut_line_pt2[0] - cut_line_pt1[0]);
    b = cut_line_pt1[1] - m * cut_line_pt1[0];
  }

  HistType *hist = (HistType *)hist_orig->Clone(Form("%s_clone", hist_orig->GetName()));
  // Optionally rebin cloned histogram to reduce noise / speed up peak finding
  if (rebinX > 1 || rebinY > 1) {
    HistType *rebinned = hist;
    // RebinX / RebinY return a new histogram when a new name is supplied; use unique names to avoid collisions
    if (rebinX > 1) {
      rebinned = (HistType *)hist->RebinX(rebinX, Form("%s_rebinX%d", hist->GetName(), rebinX));
    }
    if (rebinY > 1) {
      rebinned = (HistType *)rebinned->RebinY(rebinY, Form("%s_rebinX%d_rebinY%d", hist->GetName(), rebinX, rebinY));
    }
    // If rebin created a new histogram, delete the original clone and use the rebinned one
    if (rebinned != hist) {
      delete hist;
      hist = rebinned;
    }
  }

  // Zero out all bins below the cut line y = m*x + b
  for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
    double x    = hist->GetXaxis()->GetBinCenter(ix);
    double ycut = m * x + b;
    for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
      double y = hist->GetYaxis()->GetBinCenter(iy);
      if (y < ycut) {
        hist->SetBinContent(ix, iy, 0.0);
        hist->SetBinError(ix, iy, 0.0);
      }
      if ((!is_int && x >= front_box_cut[0] && (y <= y_cut || (y >= y_cut * 2 && y_cut > 0))) ||
          (is_int && x >= front_box_cut[0] &&
           (y <= y_cut * int_amp_ratio || (y >= y_cut * 2 * int_amp_ratio && y_cut > 0)))) {
        hist->SetBinContent(ix, iy, 0.0);
        hist->SetBinError(ix, iy, 0.0);
      }
    }
  }

  TGraph *graph = new TGraph();
  int nPoints   = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double x   = hist->GetXaxis()->GetBinCenter(i);
    TH1D *proj = hist->ProjectionY(Form("h_%s_bin_%i", hist->GetName(), i), i, i, "e");
    // proj->Smooth();
    TSpectrum *spec = new TSpectrum();
    spec->Search(proj, 1, "goff", 0.1);
    if (spec->GetNPeaks() >= 1) {
      // Get the highest peak
      double peak        = spec->GetPositionX()[0];
      double peak_height = spec->GetPositionY()[0];
      for (int j = 1; j < spec->GetNPeaks(); ++j) {
        if (spec->GetPositionY()[j] > peak_height) {
          peak        = spec->GetPositionX()[j];
          peak_height = spec->GetPositionY()[j];
        }
      }
      double y_line = m * x + b;
      double diff   = peak - y_line;
      if (diff < 0)
        diff = -diff;
      double tol = line_pt_tol;
      if (is_int)
        tol *= int_amp_ratio;
      if (diff >= tol) {
        graph->SetPoint(nPoints++, x, peak);
      }
      // graph->SetPoint(nPoints++, x, peak);
      // cout << "Peak found at: " << peak << endl;
    } /*else{
             cout << "No peak found for bin " << i << endl;
     }*/
    delete proj;
    delete spec;
  }
  return graph;
}

void fit_lad_edep2() {
  // Set ROOT to batch mode to suppress GUI
  gROOT->SetBatch(kTRUE);
  // TFile *file = TFile::Open("lad_edep_plots_22609_H.root", "READ");
  TFile *file = TFile::Open("lad_edep_plots_FT_post_gain_H.root", "READ");
  // TFile *file = TFile::Open("lad_edep_plots_LD2_setting3_P.root", "READ");
  // TFile *file = TFile::Open("lad_edep_plots_C3_23105_23109_P.root", "READ");
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open ROOT file.\n");
    return;
  }

  TFile *outfile = TFile::Open("fit_lad_edep_FT_post_gain.root", "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    printf("Error: Cannot create output ROOT file.\n");
    file->Close();
    return;
  }
  string hist_names[2] = {"KIN/TDC_DIFF_VS_ADC_AMP/c_TDC_DIFF_VS_ADC_AMP_plane_%s",
                          "KIN/TDC_DIFF_VS_ADC_INT/c_TDC_DIFF_VS_ADC_INT_plane_%s"};

  double raw_adc_max[2][N_PLANES][N_PADDLES] = {0};
  double raw_dx[2][N_PLANES][N_PADDLES]      = {0};
  for (int i_hist = 0; i_hist < 2; ++i_hist) {
    string hist_name = hist_names[i_hist];
    for (int i_plane = 0; i_plane < N_PLANES; ++i_plane) {

      TCanvas *canvas                     = (TCanvas *)file->Get(Form(hist_name.c_str(), plane_names[i_plane].c_str()));
      TList *primitives                   = canvas->GetListOfPrimitives();
      TH2D *hist_arr[N_PLANES][N_PADDLES] = {nullptr};

      TH2D *hist = nullptr;
      int count  = 0;
      for (TObject *obj : *primitives) {
        if (obj->InheritsFrom("TPad")) {
          TPad *pad            = (TPad *)obj;
          TList *padPrimitives = pad->GetListOfPrimitives();
          for (TObject *padObj : *padPrimitives) {
            if (padObj->InheritsFrom("TH2")) {
              hist_arr[0][count++] = (TH2D *)padObj;
              // Now you can use 'hist' as before
            }
          }
        }
      }

      string edep_type = (i_hist == 0) ? "amp" : "int";
      TCanvas *c1      = new TCanvas(Form("edep_%s_plane_%s", edep_type.c_str(), plane_names[i_plane].c_str()),
                                     Form("Edep %s Plane %s", edep_type.c_str(), plane_names[i_plane].c_str()));
      c1->Divide(4, 3);
      c1->Update();

      int pad_idx = 0;
      for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < N_PADDLES; ++j) {
          TH2D *hist = hist_arr[i][j];
          if (!hist)
            continue;
          c1->cd(++pad_idx);
          // TProfile *profile = new TProfile(TString::Format("profile_%d_%d_%d",i_hist,i_plane,j), "Profile above
          // line", hist->GetNbinsX(),
          //                                  hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
          // TProfile *profile = hist->ProfileX(TString::Format("profile_%d_%d_%d", i_hist, i_plane, j), 1, -1, "e");
          if (hist) {
            int maxBinX = 1, maxBinY = 1;
            double maxContent = -1.0;
            for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
              for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
                double c = hist->GetBinContent(ix, iy);
                if (c > maxContent) {
                  maxContent = c;
                  maxBinX    = ix;
                  maxBinY    = iy;
                }
              }
            }
            double x_max                    = hist->GetXaxis()->GetBinCenter(maxBinX);
            double y_max                    = hist->GetYaxis()->GetBinCenter(maxBinY);
            max_edep[i_hist][i_plane][j][0] = x_max; // x position of max bin
            max_edep[i_hist][i_plane][j][1] = y_max; // y position of max bin
          }
          int y_max_box_cut = i_plane % 2 == 0 ? 3 * max_edep[i_hist][i_plane][j][1] : 0;
          TGraph *profile   = get_peak_graph(hist, i_hist, y_max_box_cut);
          if (!profile) {
            cout << "No profile found for plane " << i_plane << ", paddle " << j << endl;
            continue;
          }
          profile->SetLineColor(kRed);
          /*for (int binx = 1; binx <= hist->GetNbinsX(); ++binx) {
            double x = hist->GetXaxis()->GetBinCenter(binx);
          // Calculate y_line as the y value on the line connecting cut_line_pt1 and cut_line_pt2 at position x
          // Calculate the average y value of the histogram for this x bin
          // You can use the TH2::ProjectionY method to get the Y projection for a given X bin,
          // then use GetMean() on the resulting TH1D to get the average Y for that X bin.
          TH1D* projY = hist->ProjectionY("_py", binx, binx);
          double y_mean = (projY->GetEntries() > 0) ? projY->GetMean() : 0;
          delete projY;
          TH1D *projX = hist->ProjectionX("_px", binx, binx);
          double x_mean = (projX->GetEntries() > 0) ? projX->GetMean() : 0;
          delete projX;

          double x1 = cut_line_pt1[0], y1 = cut_line_pt1[1];
          double x2 = cut_line_pt2[0], y2 = cut_line_pt2[1];
          double m      = (y2 - y1) / (x2 - x1);
          double b      = y1 - m * x1;
          double y_line = m * x + b;
          for (int biny = 1; biny <= hist->GetNbinsY(); ++biny) {
          double y = hist->GetYaxis()->GetBinCenter(biny);
          //if (y > y_line) {
          double content = hist->GetBinContent(binx, biny);
          // Fill the profile with y, weighted by the bin content
          for (int k = 0; k < int(content); ++k) {
          profile->Fill(x, y);
          }
          //}
          }
          }*/

          // Optionally rebin before drawing (change factors below as needed)
          TH2D *hist_to_draw  = hist;
          TH2D *hist_rebinned = nullptr;
          if (rebinX > 1 || rebinY > 1) {
            hist_rebinned = (TH2D *)hist->Clone(Form("%s_rebinned", hist->GetName()));
            if (rebinX > 1)
              hist_rebinned->RebinX(rebinX);
            if (rebinY > 1)
              hist_rebinned->RebinY(rebinY);
            hist_to_draw = hist_rebinned;
          }

          hist_to_draw->Draw("colz");

          // Note: keep hist_rebinned alive while the canvas exists to avoid removing the drawn object.
          // If you want to free memory, delete hist_rebinned after writing/saving the canvas.
          if (profile) {
            profile->SetMarkerStyle(2);
            profile->SetMarkerColor(kRed);
            profile->Draw("p");
            TF1 *edep_func;
            if (i_plane % 2 == 0) {
              edep_func = new TF1("edep_fit_func_profile_quad_front", edep_fit_func_profile_quad_front,
                                  fit_func_xmin_front, fit_func_xmax_front, 4);
              // edep_func->SetParameter(2, max_edep[i_hist][i_plane][j][0]); // x0
              // edep_func->SetParameter(1, 6); // y0
              // edep_func->FixParameter(0, -2); // x0

              edep_func->FixParameter(2, max_edep[i_hist][i_plane][j][0]); // x0
              edep_func->FixParameter(3, max_edep[i_hist][i_plane][j][1]); // y0
            } else {
              edep_func = new TF1("edep_fit_func_profile_quad_back", edep_fit_func_profile_quad_back,
                                  fit_func_xmin_back, fit_func_xmax_back, 4);
              edep_func->FixParameter(2, max_edep[i_hist][i_plane][j][0]); // x0
              edep_func->FixParameter(3, max_edep[i_hist][i_plane][j][1]); // y0
            }

            // TF1 *edep_func = new TF1("edep_fit_func_profile", edep_fit_func_profile, fit_func_xmin_back,
            // fit_func_xmax_back, 1); Initial guess for scaling
            profile->Fit(edep_func, "QNR");
            edep_func->SetLineColor(kRed);
            edep_func->SetLineWidth(2);
            edep_func->Draw("same");
            // profile->SetLineColor(kRed);
            // profile->Draw("same");
            // Draw the cut line
            // Calculate slope (m) and intercept (b) for the cut line
            double m, b;
            if (i_hist == 0) {
              m = (cut_line_pt2[1] - cut_line_pt1[1]) / (cut_line_pt2[0] - cut_line_pt1[0]);
              b = cut_line_pt1[1] - m * cut_line_pt1[0];
            } else {
              m = (cut_line_pt2[1] - cut_line_pt1[1]) * int_amp_ratio / (cut_line_pt2[0] - cut_line_pt1[0]);
              b = cut_line_pt1[1] * int_amp_ratio - m * cut_line_pt1[0];
            }
            TLine *yline = new TLine(fit_func_xmin_back, m * fit_func_xmin_back + b, fit_func_xmax_back,
                                     m * fit_func_xmax_back + b);
            yline->SetLineColor(kBlack);
            yline->SetLineStyle(2); // dashed
            yline->SetLineWidth(2);
            yline->Draw("same");
            if (i_plane % 2 == 0) {
              raw_adc_max[i_hist][i_plane][j] = edep_func->GetMaximum(fit_func_xmin_front, fit_func_xmax_front);
            } else {
              raw_adc_max[i_hist][i_plane][j] = edep_func->GetMaximum(fit_func_xmin_back, fit_func_xmax_back);
            }

            raw_dx[i_hist][i_plane][j] = edep_func ? edep_func->GetParameter(1) : 0.0;
            cout << "Plane: " << plane_names[i_plane] << ", Paddle: " << j
                 << ", Raw ADC Max: " << raw_adc_max[i_hist][i_plane][j] << ", dx: " << raw_dx[i_hist][i_plane][j]
                 << endl;
            cout << "Plane: " << plane_names[i_plane] << ", Paddle: " << j << ", Par 0: " << edep_func->GetParameter(0)
                 << ", Par 1: " << edep_func->GetParameter(1) << endl;
          }
          // delete profile; // Clean up the profile after drawing
          if (pad_idx >= 12)
            break;
        }
        if (pad_idx >= 12)
          break;
      }
      outfile->cd();
      // c1->Update();
      c1->Write();
      delete c1; // Clean up the canvas
    }
  }

  // Plot raw_dx vs paddle number for each plane (amp & int on same axes) and save canvases to outfile
  for (int i_plane = 0; i_plane < N_PLANES; ++i_plane) {
    TCanvas *c_dx = new TCanvas(Form("dx_plane_%s", plane_names[i_plane].c_str()),
                                Form("raw_dx Plane %s (amp & int)", plane_names[i_plane].c_str()), 900, 400);
    c_dx->cd();

    TGraph *g_amp = new TGraph();
    TGraph *g_int = new TGraph();
    int nAmp = 0, nInt = 0;
    double ymin = 1e12, ymax = -1e12;

    for (int ip = 0; ip < N_PADDLES; ++ip) {
      double x = ip; // paddle number (0-based)

      // amp (i_hist == 0)
      double y_amp = raw_dx[0][i_plane][ip];
      g_amp->SetPoint(nAmp++, x, y_amp);
      if (y_amp < ymin)
        ymin = y_amp;
      if (y_amp > ymax)
        ymax = y_amp;

      // int (i_hist == 1)
      double y_int = raw_dx[1][i_plane][ip];
      g_int->SetPoint(nInt++, x, y_int);
      if (y_int < ymin)
        ymin = y_int;
      if (y_int > ymax)
        ymax = y_int;
    }

    if (nAmp == 0 && nInt == 0) {
      delete g_amp;
      delete g_int;
      delete c_dx;
      continue;
    }

    // Add small margin
    double ymargin = (ymax - ymin) * 0.12;
    if (ymargin == 0)
      ymargin = 1.0;
    ymin -= ymargin;
    ymax += ymargin;

    // Style graphs
    g_amp->SetMarkerStyle(21);
    g_amp->SetMarkerSize(1.2);
    g_amp->SetLineWidth(2);
    g_amp->SetMarkerColor(kBlue);
    g_amp->SetLineColor(kBlue);

    g_int->SetMarkerStyle(22);
    g_int->SetMarkerSize(1.2);
    g_int->SetLineWidth(2);
    g_int->SetMarkerColor(kRed);
    g_int->SetLineColor(kRed);

    // Draw frame to control axis ranges and labels
    TH1F *frame = new TH1F(Form("frame_dx_%s", plane_names[i_plane].c_str()),
                           Form("Plane %s : amp & int;Paddle (index);dx", plane_names[i_plane].c_str()), N_PADDLES,
                           -0.5, N_PADDLES - 0.5);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->Draw();

    // Draw graphs on same axes
    g_amp->Draw("P SAME");
    g_int->Draw("P SAME");

    // Optional: draw a horizontal zero line
    TLine *zero = new TLine(-0.5, 0.0, N_PADDLES - 0.5, 0.0);
    zero->SetLineStyle(3);
    zero->SetLineColor(kBlack);
    zero->Draw();

    // Legend
    TLegend *leg = new TLegend(0.70, 0.75, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g_amp, "amp", "p");
    leg->AddEntry(g_int, "int", "p");
    leg->Draw();

    // Write and cleanup (keep graphs alive until canvas is written)
    outfile->cd();
    c_dx->Write();

    // delete owned objects (frame, line, legend are not owned by canvas after Write in batch, delete them)
    delete frame;
    delete zero;
    delete leg;
    delete c_dx;
    // g_amp and g_int are drawn on the canvas; they will be cleaned up by ROOT when the file is closed / canvas deleted
  }

  outfile->Close();
  file->Close();

  // Write fit results to param file

  ofstream outPARAM(Form("../../PARAM/LAD/HODO/ladhodo_edep.param"));
  outPARAM << "; LAD Hodoscope Edep Fit Parameters\n";
  outPARAM << endl << endl;

  outPARAM << "lladhodo_adcAmp2MeV = ";
  for (int i_paddle = 0; i_paddle < N_PADDLES; ++i_paddle) {
    for (int i_plane = 0; i_plane < N_PLANES; ++i_plane) {
      if (i_plane == 0) {
        if (i_paddle == 0) {
          outPARAM << true_adc_MeV / raw_adc_max[0][i_plane][i_paddle];
        } else {
          outPARAM << setw(30) << true_adc_MeV / raw_adc_max[0][i_plane][i_paddle];
        }
      } else {
        outPARAM << ", " << setw(15) << true_adc_MeV / raw_adc_max[0][i_plane][i_paddle];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << "lladhodo_adcInt2MeV = ";
  for (int i_paddle = 0; i_paddle < N_PADDLES; ++i_paddle) {
    for (int i_plane = 0; i_plane < N_PLANES; ++i_plane) {
      if (i_plane == 0) {
        if (i_paddle == 0) {
          outPARAM << true_adc_MeV / raw_adc_max[1][i_plane][i_paddle];
        } else {
          outPARAM << setw(30) << true_adc_MeV / raw_adc_max[1][i_plane][i_paddle];
        }
      } else {
        outPARAM << ", " << setw(15) << true_adc_MeV / raw_adc_max[1][i_plane][i_paddle];
      }
    }
    outPARAM << fixed << endl;
  }

  delete file; // Clean up the input file
  return;
}