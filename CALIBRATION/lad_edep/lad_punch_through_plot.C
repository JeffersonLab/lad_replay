// lad_punch_through_plot.C
//
// Reads the output ROOT file from lad_punch_through_cut.C and writes a new
// file with identical canvases, each overlaid with the punch-through cut line:
//   (3, 1000) -> (3, 40) -> (5.5, 0)
//
// Usage:
//   root -l -b -q 'lad_punch_through_plot.C("punch_through_cut.root","punch_through_plot.root")'

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <iostream>

const char *DEFAULT_IN_FILE  = "files/punch_through_cut/punch_through_cut.root";
const char *DEFAULT_OUT_FILE = "files/punch_through_cut/punch_through_plot.root";

// Run from the lad_edep directory:
//   root -l -b -q
//   'lad_punch_through_plot.C("files/punch_through_cut/punch_through_cut.root","files/punch_through_cut/punch_through_plot.root")'

// Draw the two-segment cut line on the current pad using data coordinates.
static void draw_cut_line() {
  TLine *l1 = new TLine(2.8, 1000., 2.8, 35.);
  TLine *l2 = new TLine(2.8, 35., 5., 0.);
  l1->SetLineColor(kRed);
  l1->SetLineWidth(2);
  l2->SetLineColor(kRed);
  l2->SetLineWidth(2);
  l1->Draw();
  l2->Draw();
}

// Recursively copy directories and canvases from src to dst,
// overlaying the cut line on every pad of every TCanvas.
static void process_dir(TDirectory *src, TDirectory *dst) {
  TIter next(src->GetListOfKeys());
  TKey *key = nullptr;
  while ((key = (TKey *)next())) {
    const TString cname = key->GetClassName();

    if (cname == "TDirectoryFile" || cname == "TDirectory") {
      TDirectory *sub_src = (TDirectory *)key->ReadObj();
      TDirectory *sub_dst = dst->mkdir(key->GetName());
      process_dir(sub_src, sub_dst);
      dst->cd();
      continue;
    }

    if (cname == "TCanvas") {
      TCanvas *c = (TCanvas *)key->ReadObj();

      // Iterate over all pads (sub-pads) in the canvas.
      TIter ipad(c->GetListOfPrimitives());
      TObject *pobj = nullptr;
      while ((pobj = ipad())) {
        if (!pobj->InheritsFrom(TPad::Class()))
          continue;
        TPad *pad = (TPad *)pobj;
        pad->cd();
        draw_cut_line();
        pad->Update();
      }

      // If the canvas itself has no sub-pads with histograms but does have
      // primitives (single-pad canvas), add the line directly.
      bool has_subpads = false;
      {
        TIter it(c->GetListOfPrimitives());
        TObject *o = nullptr;
        while ((o = it())) {
          if (o->InheritsFrom(TPad::Class())) {
            has_subpads = true;
            break;
          }
        }
      }
      if (!has_subpads) {
        c->cd();
        draw_cut_line();
        c->Update();
      }

      dst->cd();
      c->Write();
      delete c;
      continue;
    }

    // Copy any other objects (histograms, etc.) as-is.
    TObject *obj = key->ReadObj();
    if (obj) {
      dst->cd();
      obj->Write(key->GetName());
      delete obj;
    }
  }
}

void lad_punch_through_plot(const char *in_file = DEFAULT_IN_FILE, const char *out_file = DEFAULT_OUT_FILE) {
  gROOT->SetBatch(kTRUE);

  TFile fin(in_file, "READ");
  if (fin.IsZombie()) {
    std::cerr << "Cannot open input file: " << in_file << "\n";
    return;
  }

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "Cannot open output file: " << out_file << "\n";
    return;
  }

  process_dir(&fin, &fout);

  fout.Close();
  fin.Close();
  std::cout << "[lad_punch_through_plot] Done. Wrote " << out_file << "\n";
}
