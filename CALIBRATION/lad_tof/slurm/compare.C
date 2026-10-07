// Bin-by-bin comparison of every histogram in two files, including histograms
// drawn inside TCanvas objects (recursing into sub-pads and directories).
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>
#include <TPad.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
static void collect(TObject *o, const std::string &path, std::map<std::string, TH1 *> &m) {
  if (auto *h = dynamic_cast<TH1 *>(o)) { m[path] = h; return; }
  if (auto *p = dynamic_cast<TPad *>(o)) {
    int i = 0;
    for (auto *x : *p->GetListOfPrimitives()) collect(x, path + "/" + x->GetName() + "#" + std::to_string(i++), m);
    return;
  }
  if (auto *d = dynamic_cast<TDirectory *>(o)) {
    for (auto *k : *d->GetListOfKeys()) {
      auto *key = (TKey *)k;
      if (key->GetCycle() != d->GetKey(key->GetName())->GetCycle()) continue;
      collect(key->ReadObj(), path + "/" + key->GetName(), m);
    }
  }
}
void compare(const char *fa, const char *fb) {
  TH1::AddDirectory(false);
  TFile a(fa), b(fb);
  std::map<std::string, TH1 *> ma, mb;
  collect(&a, "", ma);
  collect(&b, "", mb);
  long nh = 0, ndiff = 0, nmiss = 0, nbins = 0;
  double maxrel = 0;
  std::string worst;
  for (auto &[k, ha] : ma) {
    auto it = mb.find(k);
    if (it == mb.end()) { if (nmiss++ < 5) std::cout << "  only in A: " << k << "\n"; continue; }
    TH1 *hb = it->second;
    ++nh;
    if (ha->GetNcells() != hb->GetNcells()) { ++ndiff; std::cout << "  binning differs: " << k << "\n"; continue; }
    bool d = false;
    for (int i = 0; i < ha->GetNcells(); ++i) {
      double x = ha->GetBinContent(i), y = hb->GetBinContent(i), ex = ha->GetBinError(i), ey = hb->GetBinError(i);
      double r = std::max(std::fabs(x - y) / std::max({std::fabs(x), std::fabs(y), 1e-12}),
                          std::fabs(ex - ey) / std::max({std::fabs(ex), std::fabs(ey), 1e-12}));
      ++nbins;
      if (r > maxrel) { maxrel = r; worst = k; }
      if (r > 1e-9) d = true;
    }
    if (d && ++ndiff && getenv("CMP_LIST")) std::cout << "  DIFF " << k << "\n";
  }
  for (auto &[k, hb] : mb) if (!ma.count(k)) { if (nmiss++ < 10) std::cout << "  only in B: " << k << "\n"; }
  std::cout << "[compare] " << fa << " vs " << fb << "\n  histograms compared=" << nh << " cells=" << nbins
            << " differing(>1e-9 rel)=" << ndiff << " missing=" << nmiss << " max rel diff=" << maxrel
            << (worst.empty() ? "" : " (" + worst + ")") << "\n";
}
