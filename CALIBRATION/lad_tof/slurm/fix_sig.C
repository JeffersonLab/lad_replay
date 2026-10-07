// Rewrite the ";runlist=<hash>" part of a merged histogram cache's signature so
// it matches <runlist> (hash computed exactly as the lad_tof macros do).
#include <TFile.h>
#include <TNamed.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
void fix_sig(const char *cache, const char *runlist) {
  std::string datlist, ln;
  std::ifstream fin(runlist);
  while (std::getline(fin, ln)) {
    size_t a = ln.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) continue;
    std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
    if (p.empty() || p[0] == '#') continue;
    datlist += p + "\n";
  }
  TFile f(cache, "UPDATE");
  auto *s = dynamic_cast<TNamed *>(f.Get("signature"));
  if (!s) { std::cerr << "[fix_sig] no signature in " << cache << "\n"; return; }
  std::string sig = s->GetTitle();
  size_t k = sig.rfind(";runlist=");
  if (k == std::string::npos) { std::cerr << "[fix_sig] no runlist field\n"; return; }
  sig = sig.substr(0, k) + ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));
  TNamed n("signature", sig.c_str());
  f.Delete("signature;*");
  n.Write();
  std::cout << "[fix_sig] " << cache << " -> " << sig.substr(k) << "\n";
}
