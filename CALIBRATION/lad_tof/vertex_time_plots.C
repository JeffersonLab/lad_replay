// vertex_time_plots.C
// Reads a runlist of ROOT files, makes per-run vertex-time histograms,
// fits a Gaussian to each, and produces run-number vs mean and width plots.
// Usage in ROOT: .x vertex_time_plots.C("../files/run-lists/all_C3_runlist.dat");

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void vertex_time_plots(const char *listfile = "../files/run-lists/all_C3_runlist.dat", const char *spec_prefix = "H",
                       double etot_cut = 0.7) {
  gStyle->SetOptFit(1111);
  gROOT->SetBatch(kTRUE);
  std::string sp(spec_prefix);

  std::ifstream infile(listfile);
  if (!infile) {
    std::cerr << "Cannot open runlist: " << listfile << "\n";
    return;
  }

  std::vector<int> runs;
  std::vector<double> means;
  std::vector<double> sigmas;

  TFile *out = TFile::Open(("vertex_time_plots_" + sp + ".root").c_str(), "RECREATE");
  std::string line;
  auto parse_runnum = [](const std::string &fname) -> int {
    // Fast path: look for "LAD_COIN_" followed by 5 digits
    const std::string key = "LAD_COIN_";
    size_t pos            = fname.find(key);
    if (pos != std::string::npos && pos + key.size() + 5 <= fname.size()) {
      std::string s = fname.substr(pos + key.size(), 5);
      bool ok       = true;
      for (char c : s)
        if (!std::isdigit((unsigned char)c)) {
          ok = false;
          break;
        }
      if (ok)
        return std::stoi(s);
    }
    // Fallback: find first run of 5 digits anywhere in the filename
    for (size_t i = 0; i + 5 <= fname.size(); ++i) {
      bool ok = true;
      for (size_t j = 0; j < 5; ++j)
        if (!std::isdigit((unsigned char)fname[i + j])) {
          ok = false;
          break;
        }
      if (ok)
        return std::stoi(fname.substr(i, 5));
    }
    return -1;
  };

  // Read all lines first to compute total and enable progress reporting
  std::vector<std::string> filelist;
  while (std::getline(infile, line)) {
    // trim whitespace
    size_t p = line.find_first_not_of(" \t\r\n");
    if (p == std::string::npos)
      continue;
    std::string fname = line.substr(p);
    if (fname.empty())
      continue;
    filelist.push_back(fname);
  }

  size_t total = filelist.size();
  if (total == 0) {
    std::cerr << "No files in runlist: " << listfile << "\n";
    out->Close();
    return;
  }

  int last_percent = -1;
  for (size_t iline = 0; iline < total; ++iline) {
    std::string fname = filelist[iline];

    int runnum = parse_runnum(fname);
    if (runnum < 0) {
      std::cerr << "Could not parse run number from: " << fname << " -- skipping\n";
      continue;
    }

    // Open the file directly and get the tree to avoid TChain overhead
    TFile *fin = TFile::Open(fname.c_str());
    if (!fin || fin->IsZombie()) {
      std::cerr << "Cannot open file: " << fname << " -- skipping\n";
      if (fin)
        fin->Close();
      continue;
    }
    TTree *tree = (TTree *)fin->Get("T");
    if (!tree) {
      std::cerr << "No tree 'T' in file " << fname << " -- skipping\n";
      fin->Close();
      continue;
    }
    // Speed: disable all branches except the few we actually need
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus((sp + ".ladkin.t_vertex").c_str(), 1);
    tree->SetBranchStatus((sp + ".ladkin.t_vertex_RFcorr").c_str(), 1);
    tree->SetBranchStatus("g.evtyp", 1);
    tree->SetBranchStatus((sp + ".cal.etottracknorm").c_str(), 1);
    Long64_t entries = tree->GetEntries();
    if (entries == 0) {
      std::cerr << "No entries for run " << runnum << " (" << fname << ") -- skipping\n";
      fin->Close();
      continue;
    }

    std::string expr = std::string("(") + sp + ".ladkin.t_vertex-" + sp + ".ladkin.t_vertex_RFcorr)";
    int evtyp_val    = (sp == "P") ? 1 : 2;
    std::string sel_base =
        Form("(g.evtyp==%d)&&(%s.ladkin.t_vertex<100000)&&(%s.ladkin.t_vertex!=%s.ladkin.t_vertex_RFcorr)", evtyp_val,
             sp.c_str(), sp.c_str(), sp.c_str());

    // Names for combined and split (by H.cal.etottracknorm) histograms
    std::string hname_comb = Form("h_%05d", runnum);
    std::string hname_high = Form("h_%05d_high", runnum);
    std::string hname_low  = Form("h_%05d_low", runnum);

    std::string draw_comb = expr + ">>" + hname_comb + "(100,-3,3)";
    std::string draw_high = expr + ">>" + hname_high + "(100,-3,3)";
    std::string draw_low  = expr + ">>" + hname_low + "(100,-3,3)";

    std::string sel_high =
        Form("(g.evtyp==1)&&(%s.ladkin.t_vertex<100000)&&(%s.cal.etottracknorm>%g)", sp.c_str(), sp.c_str(), etot_cut);
    std::string sel_low =
        Form("(g.evtyp==1)&&(%s.ladkin.t_vertex<100000)&&(%s.cal.etottracknorm<%g)", sp.c_str(), sp.c_str(), etot_cut);

    // Draw histograms in the input file directory so ROOT creates them there
    fin->cd();
    tree->Draw(draw_comb.c_str(), sel_base.c_str(), "goff");
    TH1 *h_comb = (TH1 *)fin->Get(hname_comb.c_str());
    if (!h_comb) {
      std::cerr << "Combined histogram not created for run " << runnum << "\n";
      fin->Close();
      continue;
    }

    // Use histogram statistics (mean, RMS) from the combined selection
    double mean  = h_comb->GetMean();
    double sigma = h_comb->GetRMS();

    runs.push_back(runnum);
    means.push_back(mean);
    sigmas.push_back(sigma);

    // Draw split histograms (may be empty for some runs)
    fin->cd();
    tree->Draw(draw_high.c_str(), sel_high.c_str(), "goff");
    TH1 *h_high = (TH1 *)fin->Get(hname_high.c_str());
    fin->cd();
    tree->Draw(draw_low.c_str(), sel_low.c_str(), "goff");
    TH1 *h_low = (TH1 *)fin->Get(hname_low.c_str());
    // progress percentage
    int percent = int(100.0 * (iline + 1) / total + 0.5);
    if (percent != last_percent) {
      std::cout << "Progress: " << percent << "% (" << (iline + 1) << "/" << total << ")\r" << std::flush;
      last_percent = percent;
    }

    // Move/copy histograms into the output file and save a canvas with two pads
    out->cd();

    TH1 *h_high_out = nullptr;
    TH1 *h_low_out  = nullptr;

    if (h_high) {
      h_high_out = (TH1 *)h_high->Clone(hname_high.c_str());
      h_high_out->SetTitle(
          Form("Run %05d (etottracknorm > %.3g);Vertex time diff (H - RFcorr);Counts", runnum, etot_cut));
    } else {
      h_high_out = new TH1D(
          hname_high.c_str(),
          Form("Run %05d (etottracknorm > %.3g);Vertex time diff (H - RFcorr);Counts", runnum, etot_cut), 100, -3, 3);
    }

    if (h_low) {
      h_low_out = (TH1 *)h_low->Clone(hname_low.c_str());
      h_low_out->SetTitle(
          Form("Run %05d (etottracknorm < %.3g);Vertex time diff (H - RFcorr);Counts", runnum, etot_cut));
    } else {
      h_low_out = new TH1D(
          hname_low.c_str(),
          Form("Run %05d (etottracknorm < %.3g);Vertex time diff (H - RFcorr);Counts", runnum, etot_cut), 100, -3, 3);
    }

    // h_high_out->Write();
    // h_low_out->Write();

    // Extract suffix after the run number in the filename (e.g. the "9" in LAD_COIN_23750_9_9_-1.root)
    std::string run_suffix = "0";
    {
      std::string runstr = std::to_string(runnum);
      size_t pos_run     = fname.find(runstr);
      if (pos_run != std::string::npos) {
        size_t pos_after = pos_run + runstr.size();
        if (pos_after < fname.size() && fname[pos_after] == '_') {
          size_t start = pos_after + 1;
          size_t end   = start;
          while (end < fname.size() && fname[end] != '_' && fname[end] != '.')
            ++end;
          if (end > start)
            run_suffix = fname.substr(start, end - start);
        }
      }
    }
    TCanvas *crun = new TCanvas(Form("c_%05d_%s", runnum, run_suffix.c_str()),
                                Form("Run %05d_%s", runnum, run_suffix.c_str()), 900, 600);
    // Draw both histograms on the left pad and etot distribution on the right pad
    h_high_out->SetLineColor(kBlue);
    h_high_out->SetLineWidth(2);
    h_low_out->SetLineColor(kRed);
    h_low_out->SetLineWidth(2);

    crun->Divide(2, 1);
    crun->cd(1);
    h_high_out->Draw();
    h_low_out->Draw("SAME");
    TLegend *leg = new TLegend(0.65, 0.75, 0.9, 0.9);
    leg->AddEntry(h_high_out, Form("%s.cal.etottracknorm > %.3g", sp.c_str(), etot_cut), "l");
    leg->AddEntry(h_low_out, Form("%s.cal.etottracknorm < %.3g", sp.c_str(), etot_cut), "l");
    leg->SetBorderSize(0);
    leg->Draw();

    // Right pad: etottracknorm histogram with cut marker
    crun->cd(2);
    std::string hetot_name = Form("h_etot_%05d", runnum);
    std::string draw_etot  = std::string(sp) + ".cal.etottracknorm>>" + hetot_name + "(100,-0.5,2)";
    fin->cd();
    tree->Draw(draw_etot.c_str(), sel_base.c_str(), "goff");
    TH1 *h_etot     = (TH1 *)fin->Get(hetot_name.c_str());
    TH1 *h_etot_out = nullptr;
    if (h_etot) {
      h_etot_out = (TH1 *)h_etot->Clone(hetot_name.c_str());
      h_etot_out->SetTitle(Form("Run %05d;%s.cal.etottracknorm;Counts", runnum, sp.c_str()));
    } else {
      h_etot_out =
          new TH1D(hetot_name.c_str(), Form("Run %05d;%s.cal.etottracknorm;Counts", runnum, sp.c_str()), 100, -0.5, 2);
    }
    out->cd();
    // h_high_out->Write();
    // h_low_out->Write();
    // h_etot_out->Write();
    crun->cd(2);
    h_etot_out->Draw();
    double ymax    = h_etot_out->GetMaximum();
    TLine *cutline = new TLine(etot_cut, 0, etot_cut, ymax * 1.05);
    cutline->SetLineColor(kGreen + 2);
    cutline->SetLineWidth(2);
    cutline->SetLineStyle(2);
    cutline->Draw();

    crun->Write();

    fin->Close();
  }

  // Sort entries by run number
  size_t n = runs.size();
  if (n == 0) {
    std::cerr << "No runs processed. Exiting.\n";
    out->Close();
    return;
  }

  std::vector<size_t> idx(n);
  for (size_t i = 0; i < n; ++i)
    idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return runs[a] < runs[b]; });

  std::vector<double> xr(n), ym(n), ys(n);
  for (size_t i = 0; i < n; ++i) {
    xr[i] = runs[idx[i]];
    ym[i] = means[idx[i]];
    ys[i] = sigmas[idx[i]];
  }

  TGraph *g_mean = new TGraph(n, &xr[0], &ym[0]);
  g_mean->SetTitle("Run number vs Gaussian mean;Run number;Mean (units of input)");
  g_mean->SetMarkerStyle(20);
  TCanvas *cmean = new TCanvas("run_vs_mean", "Run vs Mean", 900, 600);
  g_mean->Draw("AP");
  cmean->Write();
  g_mean->Write("run_vs_mean");
  cmean->Write();

  TGraph *g_sigma = new TGraph(n, &xr[0], &ys[0]);
  g_sigma->SetTitle("Run number vs Gaussian sigma;Run number;Sigma (units of input)");
  g_sigma->SetMarkerStyle(21);
  TCanvas *csig = new TCanvas("run_vs_sigma", "Run vs Sigma", 900, 600);
  g_sigma->Draw("AP");
  csig->Write();
  g_sigma->Write("run_vs_sigma");
  csig->Write();

  // out->Write();
  out->Close();

  std::cout << "Finished. Created vertex_time_plots.root with canvases, histograms, and graphs." << std::endl;
}
