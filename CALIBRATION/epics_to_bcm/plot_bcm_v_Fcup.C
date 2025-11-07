#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

int plot_bcm_v_Fcup() {
  std::ifstream infile("segment_averages.txt");
  if (!infile) {
    std::cerr << "Error opening file segment_averages.txt\n";
    return 1;
  }

  // For each sensor, store a vector of (sum, count) pairs representing different averages
  std::map<std::string, std::vector<std::pair<double, int>>> sensor_groups;
  std::string line;

  while (std::getline(infile, line)) {
    // Skip comments and empty lines
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string sensor;
    int startIdx, endIdx;
    double avg, rms;

    if (!(iss >> sensor >> startIdx >> endIdx >> avg >> rms)) {
      // Try to handle extra spaces in sensor name column
      iss.clear();
      iss.str(line);
      iss >> sensor;
      iss >> startIdx >> endIdx >> avg >> rms;
    }

    auto &groups = sensor_groups[sensor];
    bool added   = false;
    for (auto &group : groups) {
      double group_avg = group.first / group.second;
      // If within 20% of existing group average, add to this group
      if (std::abs(avg - group_avg) <= 0.3 * group_avg) {
        group.first += avg;
        group.second += 1;
        added = true;
        break;
      }
    }
    if (!added) {
      // Start a new group
      groups.emplace_back(avg, 1);
    }
  }

  std::cout << std::fixed << std::setprecision(8);
  for (const auto &kv : sensor_groups) {
    const std::string &sensor = kv.first;
    const auto &groups        = kv.second;
    for (size_t i = 0; i < groups.size(); ++i) {
      double mean = groups[i].first / groups[i].second;
      std::cout << sensor << " average group " << (i + 1) << ": " << mean << " (n=" << groups[i].second << ")"
                << std::endl;
    }
  }

  // Remove groups with 1 or fewer values and sort remaining groups by mean
  for (auto &kv : sensor_groups) {
    auto &groups = kv.second;
    // Remove groups with 1 or fewer values
    groups.erase(std::remove_if(groups.begin(), groups.end(),
                                [](const std::pair<double, int> &group) { return group.second <= 1; }),
                 groups.end());
    // Sort groups by mean (ascending)
    std::sort(groups.begin(), groups.end(), [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
      return (a.first / a.second) < (b.first / b.second);
    });
  }

  std::vector<double> ibcm3H04A_means, ibcm3H04B_means, ibcm3H04C_means;
  std::vector<double> FCupsCORRECTEd_means;

  // Collect means for ibcm3H04B and FCupsCORRECTEd
  auto find_means = [&](const std::string &sensor, std::vector<double> &means) {
    auto it = sensor_groups.find(sensor);
    if (it != sensor_groups.end()) {
      for (const auto &group : it->second) {
        means.push_back(group.first / group.second);
      }
    }
  };

  find_means("ibcm3H04B(frq)", ibcm3H04B_means);
  find_means("FCupsCORRECTED", FCupsCORRECTEd_means);

  // Ensure both vectors have the same size
  size_t n = std::min(ibcm3H04B_means.size(), FCupsCORRECTEd_means.size());
  if (n == 0) {
    std::cerr << "No matching groups found for ibcm3H04B and FCupsCORRECTEd\n";
    return 2;
  }
  // Truncate the larger vector to match the smaller size
  ibcm3H04B_means.resize(n);
  FCupsCORRECTEd_means.resize(n);

  // Create a histogram of ibcm3H04B vs FCupsCORRECTEd
  TFile *fout = new TFile("bcm_vs_fcup.root", "RECREATE");
  TCanvas *c1 = new TCanvas("c1", "ibcm3H04B vs FCupsCORRECTEd", 800, 600);
  // Create a TGraph of ibcm3H04B vs FCupsCORRECTEd
  TGraph *gr = new TGraph(n, FCupsCORRECTEd_means.data(), ibcm3H04B_means.data());
  gr->SetTitle("ibcm3H04B vs FCupsCORRECTED;FCupsCORRECTED;ibcm3H04B");
  gr->SetMarkerStyle(20);
  gr->Draw("AP");
  c1->Write();
  gr->Write("bcm_vs_fcup_graph");
  fout->Close();
  delete fout;

  return 0;
}