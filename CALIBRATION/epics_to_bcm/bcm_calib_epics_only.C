#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLine.h>
#include <TROOT.h>
#include <TTree.h>
#include <sstream>

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

double seg_start_end_threshold = 2.0;  // Threshold to start/end a segment
double tol                     = 0.30; // Tolerance for stability check
double min_value               = 0.05; // Minimum value to consider for segments

// Structure to hold segment information
struct segment_info {
  int start_idx;
  int end_idx;
  double avg;
  double rms;
};

int time_index_start = 0;
int time_index_end   = 2500;

double frq2I(double frq, const bcm_calib_param &param) {
  // Convert frequency to current using calibration parameters
  return (frq - param.A) * param.B;
}

double I2frq(double I, const bcm_calib_param &param) {
  // Convert current to frequency using calibration parameters
  return I / param.B + param.A;
}

void bcm_calib_epics_only() {
  // Read the file as a text file
  FILE *fp = fopen("epics_bcm_with_current.dat", "r");
  if (!fp) {
    printf("Error: Cannot open input file epics_bcm_raw.dat\n");
    return;
  }

  char line[1024];
  // Read header line
  if (!fgets(line, sizeof(line), fp)) {
    printf("Error: Cannot read header line\n");
    fclose(fp);
    return;
  }

  std::vector<std::string> headers;
  std::istringstream header_stream(line);
  std::string header;
  while (header_stream >> header) {
    headers.push_back(header);
  }

  // Prepare map from header to vector of strings (for Date) or doubles
  std::map<std::string, std::vector<std::string>> date_map;
  std::map<std::string, std::vector<double>> data_map;
  for (const auto &h : headers) {
    if (h == "Date")
      date_map[h] = std::vector<std::string>();
    else
      data_map[h] = std::vector<double>();
  }

  // Read data lines, only between time_index_start and time_index_end
  int line_index = 0;
  while (fgets(line, sizeof(line), fp)) {
    if (line_index < time_index_start) {
      ++line_index;
      continue;
    }
    if (line_index > time_index_end)
      break;

    std::istringstream data_stream(line);
    for (size_t i = 0; i < headers.size(); ++i) {
      if (headers[i] == "Date") {
        std::string date_part1, date_part2;
        data_stream >> date_part1 >> date_part2;
        std::string full_date = date_part1 + " " + date_part2;
        date_map["Date"].push_back(full_date);
      } else {
        double value;
        data_stream >> value;
        data_map[headers[i]].push_back(value);
      }
    }
    ++line_index;
  }

  fclose(fp);

  // Map from variable name to vector of segment_info
  std::map<std::string, std::vector<segment_info>> segments_map;

  for (const auto &kv : data_map) {
    const std::string &var            = kv.first;
    const std::vector<double> &values = kv.second;
    int n                             = values.size();

    bool in_segment = false;
    bool is_rising  = false;
    int seg_start   = -1;
    std::vector<double> seg_values;

    for (int i = 1; i < n; ++i) {
      double prev = values[i - 1];
      double curr = values[i];

      if (!in_segment && curr > 2.0 * prev && curr > min_value) {
        is_rising = true; // Start of a rising segment
      }

      bool start_data = false;
      if (is_rising) {
        // Check next n points for stability within 5% tolerance
        const int n_check = 3;
        start_data        = true;
        for (int j = 0; j < n_check && (i + j) < n; ++j) {
          double ref = curr;
          double val = values[i + j];
          if (fabs(val - ref) > tol * fabs(ref)) {
            start_data = false;
            break;
          }
        }
      }

      if (start_data) {
        // Start segment
        is_rising  = false; // Reset rising flag
        in_segment = true;
        seg_start  = i;
        seg_values.clear();
        seg_values.push_back(curr);
      } else if (in_segment) {
        seg_values.push_back(curr);

        // End segment if the next n_check points are not within 5% tolerance
        const int n_check = 2;
        bool end_segment  = false;
        for (int j = 0; j < n_check && (i + j) < n; ++j) {
          double ref = curr;
          double val = values[i + j];
          if (fabs(val - ref) > tol * fabs(ref)) {
            end_segment = true;
            break;
          }
        }
        if (end_segment) {
          // End segment
          int seg_end = i;
          // Calculate average and rms
          double sum = 0.0, sum2 = 0.0;
          for (double v : seg_values) {
            sum += v;
            sum2 += v * v;
          }
          double avg = sum / seg_values.size();
          double rms = sqrt(sum2 / seg_values.size() - avg * avg);

          segments_map[var].push_back({seg_start, seg_end, avg, rms});
          in_segment = false;
        }
      }
    }
    // If still in segment at the end, close it
    if (in_segment && !seg_values.empty()) {
      int seg_end = n - 1;
      double sum = 0.0, sum2 = 0.0;
      for (double v : seg_values) {
        sum += v;
        sum2 += v * v;
      }
      double avg = sum / seg_values.size();
      double rms = sqrt(sum2 / seg_values.size() - avg * avg);
      segments_map[var].push_back({seg_start, seg_end, avg, rms});
    }
  }

  // Create a ROOT file to save the plots
  TFile *outfile = new TFile("epics_bcm_plots.root", "RECREATE");

  // For each variable except "Date", make a TGraph vs time index
  for (const auto &kv : data_map) {
    const std::string &var            = kv.first;
    const std::vector<double> &values = kv.second;
    int n                             = values.size();

    // Use the index as a proxy for time (since Date is a string)
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i)
      x[i] = i;

    TGraph *g = new TGraph(n, x.data(), values.data());
    g->SetTitle((var + " vs Time Index").c_str());
    g->GetXaxis()->SetTitle("Time Index");
    g->GetYaxis()->SetTitle(var.c_str());

    // Optionally, draw and save to canvas
    TCanvas *c = new TCanvas();
    g->Draw("ALP");

    if (segments_map.count(var)) {
      for (const auto &seg : segments_map[var]) {
        double x1   = seg.start_idx;
        double x2   = seg.end_idx;
        double y    = seg.avg;
        TLine *line = new TLine(x1, y, x2, y);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw();
      }
    }
    c->Write((var + "_canvas").c_str());
    // g->Write((var + "_graph").c_str());
    delete c;
  }

  outfile->Close();
  delete outfile;

  FILE *seg_fp = fopen("segment_averages.txt", "w");
  if (seg_fp) {
    // Write header
    fprintf(seg_fp, "# Variable    StartIdx    EndIdx    Average        RMS\n");
    fprintf(seg_fp, "# Columns: Variable name, segment start index, segment end index, segment average, segment RMS\n");
    for (const auto &kv : segments_map) {
      const std::string &var                = kv.first;
      const std::vector<segment_info> &segs = kv.second;
      for (const auto &seg : segs) {
        fprintf(seg_fp, "%-12s %8d %8d %14.8f %14.8f\n", var.c_str(), seg.start_idx, seg.end_idx, seg.avg, seg.rms);
      }
    }
    fclose(seg_fp);
  } else {
    printf("Error: Cannot open output file segment_averages.txt\n");
  }
}