#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TROOT.h>
#include <string>
const int nPlanes                     = 6;
const std::string planeNames[nPlanes] = {"000", "001", "100", "101", "200", "REFBAR"};
const int nSides                      = 2;
const int maxPaddles[nPlanes]         = {11, 11, 11, 11, 11, 1};
const int run_numbers[]               = {300, 22824, 23249, 23798, 23451};
const std::string run_dates[]         = {"23-Mar-25", "30-May-25", "17-Jun-25", "14-Jul-25", "27-Jun-25"};
const std::string param_names[]       = {"lladhodo_cableFit", "lladhodo_LCoeff"};
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

std::vector<std::vector<double>> get_params_for_run(int run_number, std::string param_name) {
  std::vector<std::vector<double>> params(nPlanes);
  for (int i = 0; i < nPlanes; ++i) {
    params[i].resize(maxPaddles[i], 0.0);
  }

  // Construct filename based on run_number
  std::string filename = "../../PARAM/LAD/HODO/ladhodo_Vpcalib_" + std::to_string(run_number) + ".param";
  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return params;
  }

  std::string line;
  bool found_section = false;
  std::vector<double> values;

  while (std::getline(infile, line)) {
    if (!found_section) {
      if (line.find(param_name) != std::string::npos) {
        found_section = true;
        // Extract values after '=' on this line
        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
          std::string after_eq = line.substr(eq_pos + 1);
          std::istringstream iss(after_eq);
          double val;
          char comma;
          while (iss >> val) {
            values.push_back(val);
            iss >> comma;
          }
        }
      }
    } else {
      // Read values from subsequent lines
      std::istringstream iss(line);
      double val;
      char comma;
      while (iss >> val) {
        values.push_back(val);
        iss >> comma;
      }
      // Stop if we've read enough values
      if (values.size() >= nPlanes * maxPaddles[0])
        break;
    }
  }

  // Fill params vector
  size_t idx        = 0;
  int maxNumPaddles = *std::max_element(maxPaddles, maxPaddles + nPlanes);
  for (int j = 0; j < maxNumPaddles; ++j) {
    for (int i = 0; i < nPlanes; ++i) {
      if (idx < values.size()) {
        if (j < maxPaddles[i]) {
          params[i][j] = values[idx++];
        } else {
          ++idx; // skip unused paddles in the flat list
        }
      } else {
        if (j < maxPaddles[i]) {
          params[i][j] = 0.0;
        }
      }
    }
  }

  return params;
}

void compare_calib() {
  gROOT->SetBatch(kTRUE);
  TFile *outfile           = new TFile("calibration_comparison_all_planes.root", "RECREATE");
  const int nRunsToCompare = sizeof(run_numbers) / sizeof(run_numbers[0]);
  const int nParams        = sizeof(param_names) / sizeof(param_names[0]);
  TCanvas *c[nParams][nPlanes];
  for (int i_param = 0; i_param < nParams; ++i_param) {
    for (int plane = 0; plane < nPlanes; ++plane) {
      c[i_param][plane] = new TCanvas(
          Form("c_%s_%s", param_names[i_param].c_str(), planeNames[plane].c_str()),
          Form("%s vs Run Number, Plane %s", param_names[i_param].c_str(), planeNames[plane].c_str()), 1200, 900);
      c[i_param][plane]->Divide(4, 3);
    }

    const int nRuns         = sizeof(run_numbers) / sizeof(run_numbers[0]);
    const int maxNumPaddles = *std::max_element(maxPaddles, maxPaddles + nPlanes);
    TGraph *graphs[nRuns][nPlanes][maxNumPaddles];
    for (int i = 0; i < nRuns; ++i) {
      for (int plane = 0; plane < nPlanes; ++plane) {
        for (int paddle = 0; paddle < maxPaddles[plane]; ++paddle) {
          graphs[i][plane][paddle] = nullptr;
        }
      }
    }

    for (int i = 0; i < nRuns; ++i) {
      std::vector<std::vector<double>> param = get_params_for_run(run_numbers[i], param_names[i_param]);
      for (int plane = 0; plane < nPlanes; ++plane) {
        for (int paddle = 0; paddle < maxPaddles[plane]; ++paddle) {
          if (param[plane][paddle] != 0.0) {
            if (!graphs[i][plane][paddle]) {
              // std::cout << "Run: " << run_numbers[i] << ", Plane: " << planeNames[plane] << ", Paddle: " << paddle +
              // 1
              //           << ", Param: " << param[plane][paddle] << std::endl;
              graphs[i][plane][paddle] = new TGraph();
              graphs[i][plane][paddle]->SetTitle(Form("Plane %s Paddle %02d", planeNames[plane].c_str(), paddle + 1));
              graphs[i][plane][paddle]->SetMarkerStyle(20 + i);
              graphs[i][plane][paddle]->SetMarkerSize(1);
              graphs[i][plane][paddle]->SetMarkerColor(i + 1);
              graphs[i][plane][paddle]->SetLineColor(i + 1);
            }
            int pointIndex = graphs[i][plane][paddle]->GetN();
            graphs[i][plane][paddle]->SetPoint(pointIndex, run_numbers[i], param[plane][paddle]);
          }
        }
      }
    }
    // Loop over all planes and paddles

    for (int plane = 0; plane < nPlanes; ++plane) {
      int padIdx = 1;
      for (int paddle = 0; paddle < maxPaddles[plane]; ++paddle) {
        c[i_param][plane]->cd(padIdx);

        // Find min/max y for this paddle across all runs
        double minY = 1e9, maxY = -1e9;
        int nPoints = 0;
        for (int i = 0; i < nRuns; ++i) {
          if (graphs[i][plane][paddle] && graphs[i][plane][paddle]->GetN() > 0) {
            for (int j = 0; j < graphs[i][plane][paddle]->GetN(); ++j) {
              double x, y;
              graphs[i][plane][paddle]->GetPoint(j, x, y);
              if (y < minY)
                minY = y;
              if (y > maxY)
                maxY = y;
              ++nPoints;
            }
          }
        }
        if (nPoints == 0) {
          ++padIdx;
          continue;
        }

        // Draw all runs for this paddle, with even x spacing
        double xSpacing = 1.0;
        for (int i = 0; i < nRuns; ++i) {
          if (graphs[i][plane][paddle] && graphs[i][plane][paddle]->GetN() > 0) {
            // Remap x to even spacing
            double y;
            double x_dummy;
            graphs[i][plane][paddle]->GetPoint(0, x_dummy, y);
            graphs[i][plane][paddle]->SetPoint(0, i * xSpacing, y);

            if (i == 0) {
              graphs[i][plane][paddle]->GetYaxis()->SetRangeUser(minY - 0.05 * fabs(minY), maxY + 0.05 * fabs(maxY));
              graphs[i][plane][paddle]->GetXaxis()->SetTitle("Run Index");
              graphs[i][plane][paddle]->Draw("AP");
            } else {
              graphs[i][plane][paddle]->Draw("P SAME");
            }
          }
        }
        // Set consistent x axis for all graphs: ticks only at run indices, labels are run numbers
        c[i_param][plane]->cd(padIdx);

        // Set up a dummy histogram for axis control
        double yMargin = 0.05 * (fabs(maxY) > 1e-6 ? fabs(maxY) : 1.0);
        TH1F *hFrame   = new TH1F(
            Form("hFrame_%s_plane%s_paddle%d", param_names[i_param].c_str(), planeNames[plane].c_str(), paddle + 1), "",
            nRuns, -0.5, nRuns - 0.5);
        hFrame->SetMinimum(minY - yMargin);
        hFrame->SetMaximum(maxY + yMargin);
        hFrame->GetXaxis()->SetTitle("Run Number");
        std::string yTitle = param_names[i_param] + " (ns)";
        hFrame->GetYaxis()->SetTitle(yTitle.c_str());
        // hFrame->GetXaxis()->SetNdivisions(nRuns + 1, kTRUE); // One tick per bin center
        // hFrame->GetXaxis()->SetNdivisions(nRuns, kTRUE); // One tick per run
        hFrame->GetXaxis()->CenterLabels(kTRUE); // Center labels in bins
        hFrame->SetStats(0);
        hFrame->Draw();

        // Set custom x axis labels to run numbers and draw them at 45 degrees
        TAxis *xaxis  = hFrame->GetXaxis();
        double labelY = hFrame->GetMinimum() - 0.01 * (hFrame->GetMaximum() - hFrame->GetMinimum());
        for (int i = 0; i < nRuns; ++i) {
          xaxis->SetBinLabel(i + 1, ""); // Remove default label
          double x = i;
          // TLatex *latex = new TLatex(x, labelY, Form("%d", run_numbers[i]));
          // Split label into two lines: run number and date
          std::string label        = Form("%d", run_numbers[i]);
          std::string date         = run_dates[i];
          std::string twoLineLabel = label + "\n" + date;
          TLatex *latex            = new TLatex(x, labelY, twoLineLabel.c_str());
          latex->SetTextAlign(33); // right, top
          latex->SetTextAngle(45);
          latex->SetTextSize(0.03);
          latex->Draw("SAME");
        }
        xaxis->SetTickLength(0.03); // Make sure ticks are visible at each run

        // Draw all graphs on top, without modifying x axis ticks or labels
        for (int i = 0; i < nRuns; ++i) {
          if (graphs[i][plane][paddle] && graphs[i][plane][paddle]->GetN() > 0) {
            graphs[i][plane][paddle]->Draw("P SAME");
          }
        }
        ++padIdx;
      }
    }
    // Save all canvases to the same output ROOT file
    for (int plane = 0; plane < nPlanes; ++plane) {
      c[i_param][plane]->Write();
    }
  }
  outfile->Close();

  delete outfile;
  return;
}