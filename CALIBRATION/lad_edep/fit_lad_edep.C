#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLine.h>
#include <TProfile.h>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

const char spec_prefix = 'H'; // Spectrometer to replay

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

const double fit_func_xmin = 2.2; // Minimum x value for fit function
const double fit_func_xmax = 6.0; // Maximum x value for fit function

const double quad_fit_param[2] = {5.31637, 12.9101};
double cut_line_pt1[2]   = {2, 120}; // Point 1 for cut line
const double cut_line_pt2[2]   = {6, 0};   // Point 2 for cut line

const double true_adc_MeV = 80; // True ADC value in MeV

double edep_fit_func_profile(double *x, double *par) {
  // x[0]: x, x[1]: y
  // par[0]: amplitude, par[1]: x-mean, par[2]: x-sigma, par[3]: y-mean, par[4]: y-sigma, par[5]: x0, par[6]: constant
  // 2D Gaussian for x < x0, else constant
  // Only one free parameter: par[0] is the y-axis scaling factor
  double y;
  if (x[0] < 3.5) {
    // Line between (1.5,40) and (3.5,180)
    double x1 = 1.5, y1 = 40.0;
    double x2 = 3.5, y2 = 180.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  } else {
    // Line between (3.5,180) and (6,20)
    double x1 = 3.5, y1 = 180.0;
    double x2 = 6.0, y2 = 20.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  }
  return par[0] * y;
}

double edep_fit_func_profile_quad_back(double *x, double *par) {
  double y;
  // Fit two quadratics: one for x < x_mid, one for x >= x_mid
  // par[0]: scaling
  // par[1]: a_left, par[2]: b_left, par[3]: c_left
  // par[4]: a_right, par[5]: b_right, par[6]: c_right
  double x_mid = 3.1;
  double x1, y1, x2, y2, a;
  if (x[0] < x_mid) {
    // Quadratic passing through (1,20) and (3.3,190)
    // One free parameter
    x1 = 1.0, y1 = 20.0;
    x2 = x_mid, y2 = 190.0;
    a = quad_fit_param[0];
  } else {
    // Quadratic passing through (3.3,190) and (7,20)
    x1 = x_mid, y1 = 190.0;
    x2 = 6.0, y2 = 10.0;
    a = quad_fit_param[1];
  }
  // double a = par[0];
  double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
  y            = a * (x[0] * x[0] - x1 * x1) + b_tmp * (x[0] - x1) + y1;

  return y * par[0]; // Scale by par[0]
}

double edep_fit_func_profile_quad_front(double *x, double *par) {
  double y;
  // Fit two quadratics: one for x < x_mid, one for x >= x_mid
  // par[0]: scaling
  // par[1]: a_left, par[2]: b_left, par[3]: c_left
  // par[4]: a_right, par[5]: b_right, par[6]: c_right

  // Quadratic passing through (3.5,180) and (6,20)
  double x1 = 2.6, y1 = 80.0;
  double x2 = 6.5, y2 = 170;
  // double a = quad_fit_param[1];
  double a = -5.82599; // Make a global variable
  // double a = par[0];
  double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
  y            = a * (x[0] * x[0] - x1 * x1) + b_tmp * (x[0] - x1) + y1;

  return y * par[0]; // Scale by par[0]
}

void fit_lad_edep() {
  // Set ROOT to batch mode to suppress GUI
  gROOT->SetBatch(kTRUE);
  // TFile *file = TFile::Open("lad_edep_plots_22609_H.root", "READ");
  TFile *file = TFile::Open("lad_edep_plots_new_timing_H.root", "READ");
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open ROOT file.\n");
    return;
  }

  TFile *outfile = TFile::Open("fit_lad_edep_new_timing.root", "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    printf("Error: Cannot create output ROOT file.\n");
    file->Close();
    return;
  }
  string hist_names[2] = {"KIN/TDC_DIFF_VS_ADC_AMP/c_TDC_DIFF_VS_ADC_AMP_plane_%s",
                          "KIN/TDC_DIFF_VS_ADC_INT/c_TDC_DIFF_VS_ADC_INT_plane_%s"};

  double raw_adc_max[2][N_PLANES][N_PADDLES] = {0};
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
          TProfile *profile = new TProfile("profile", "Profile above line", hist->GetNbinsX(),
                                           hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

          for (int binx = 1; binx <= hist->GetNbinsX(); ++binx) {
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

            double x1 = x_mean, y1 = y_mean*4;
            double x2 = cut_line_pt2[0], y2 = cut_line_pt2[1];
            double m      = (y2 - y1) / (x2 - x1);
            double b      = y1 - m * x1;
            double y_line = m * x + b;
            for (int biny = 1; biny <= hist->GetNbinsY(); ++biny) {
              double y = hist->GetYaxis()->GetBinCenter(biny);
              if (y > y_line) {
                double content = hist->GetBinContent(binx, biny);
                // Fill the profile with y, weighted by the bin content
                for (int k = 0; k < int(content); ++k) {
                  profile->Fill(x, y);
                }
              }
            }
          }

          hist->Draw("colz");
          if (profile) {
            // profile->Draw("same");
            TF1 *edep_func;
            if (i_plane % 2 == 0) {
              edep_func = new TF1("edep_fit_func_profile_quad_front", edep_fit_func_profile_quad_front, fit_func_xmin,
                                  fit_func_xmax, 1);
            } else {
              edep_func = new TF1("edep_fit_func_profile_quad_back", edep_fit_func_profile_quad_back, fit_func_xmin,
                                  fit_func_xmax, 1);
            }

            // TF1 *edep_func = new TF1("edep_fit_func_profile", edep_fit_func_profile, fit_func_xmin, fit_func_xmax,
            // 1); Initial guess for scaling
            profile->Fit(edep_func, "QNR");
            edep_func->SetLineColor(kRed);
            edep_func->SetLineWidth(2);
            edep_func->Draw("same");
            // profile->SetLineColor(kRed);
            // profile->Draw("same");
            // Draw the cut line
            // Calculate slope (m) and intercept (b) for the cut line
            double m     = (cut_line_pt2[1] - cut_line_pt1[1]) / (cut_line_pt2[0] - cut_line_pt1[0]);
            double b     = cut_line_pt1[1] - m * cut_line_pt1[0];
            TLine *yline = new TLine(fit_func_xmin, m * fit_func_xmin + b, fit_func_xmax, m * fit_func_xmax + b);
            yline->SetLineColor(kBlack);
            yline->SetLineStyle(2); // dashed
            yline->SetLineWidth(2);
            yline->Draw("same");

            raw_adc_max[i_hist][i_plane][j] = edep_func->GetMaximum(fit_func_xmin, fit_func_xmax);
            cout << "Plane: " << plane_names[i_plane] << ", Paddle: " << j
                 << ", Raw ADC Max: " << raw_adc_max[i_hist][i_plane][j] << endl;
          }
          delete profile; // Clean up the profile after drawing
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
        }
        else{
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
        }
        else{
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