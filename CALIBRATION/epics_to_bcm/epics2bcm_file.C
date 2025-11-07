#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLine.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
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

double time_index_start = 400;
double time_index_end   = 2500;

double frq2I(double frq, const bcm_calib_param &param) {
  // Convert frequency to current using calibration parameters
  return (frq - param.A) * param.B;
}

double I2frq(double I, const bcm_calib_param &param) {
  // Convert current to frequency using calibration parameters
  return I / param.B + param.A;
}

void epics2bcm_file() {
  // Read the file as a text file
  FILE *fp = fopen("epics_bcm_raw.dat", "r");
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

  std::vector<double> ibcm1, ibcm2, ibcm4a, ibcm4b, ibcm4c;
  const auto &bcm1  = data_map["ibcm1"];
  const auto &bcm2  = data_map["ibcm2"];
  const auto &bcm4a = data_map["ibcm3H04A"];
  const auto &bcm4b = data_map["ibcm3H04B"];
  const auto &bcm4c = data_map["ibcm3H04C"];
  size_t nrows      = bcm1.size();

  for (size_t i = 0; i < nrows; ++i) {
    ibcm1.push_back(I2frq(bcm1[i], param_bcm1));
    ibcm2.push_back(I2frq(bcm2[i], param_bcm2));
    ibcm4a.push_back(I2frq(bcm4a[i], param_bcm4a));
    ibcm4b.push_back(I2frq(bcm4b[i], param_bcm4b));
    ibcm4c.push_back(I2frq(bcm4c[i], param_bcm4c));
  }

  // Write to new text file
  FILE *outf = fopen("epics_bcm_with_current.dat", "w");
  if (!outf) {
    printf("Error: Cannot open output file epics_bcm_with_current.dat\n");
    return;
  }

  // Write header
  for (const auto &h : headers) {
    fprintf(outf, "%-20s", h.c_str());
  }
  fprintf(outf, "%-20s%-20s%-20s%-20s%-20s\n", "ibcm1(frq)", "ibcm2(frq)", "ibcm3H04A(frq)", "ibcm3H04B(frq)",
          "ibcm3H0C(frq)");

  // Write data
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < headers.size(); ++j) {
      if (headers[j] == "Date") {
        fprintf(outf, "%-20s", date_map["Date"][i].c_str());
      } else {
        fprintf(outf, "%-20.8f", data_map[headers[j]][i]);
      }
    }
    fprintf(outf, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", ibcm1[i], ibcm2[i], ibcm4a[i], ibcm4b[i], ibcm4c[i]);
  }

  fclose(outf);
}