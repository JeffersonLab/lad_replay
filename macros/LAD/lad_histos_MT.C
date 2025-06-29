// Created by ehingerl on 04/01/25.
// Produces low level histograms for LAD Hodo and GEMs. Had to use this, as DEF-Files failed to make certain histograms.
// 1D & 2D histograms can be defined using the scructs. The script will loop over all leaves, and evaluate cuts on a
// per-leaf basis
// Currently implements multi-threading to optimize performance over lad_histos_MT.C

#include <ROOT/RConfig.hxx> // Include the header for ROOT::EnableThreadSafety
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>
#include <functional> // Ensure std::function is included
#include <iomanip>
#include <iostream>
#include <string>
#include <thread> // Include thread for std::thread::hardware_concurrency
#include <vector>
// #include <sys/resource.h> // Include for struct rlimit and related functions

const int NDDATA_MAX = 20000; // Maximum number of data points to read. GEM's need at least >10k (10k throws errors for
                              // some events). 100k might be overkill, but throws now errors.
const int MAX_THREADS      = 20;               // Maximum number of threads
const int NEVTS_PER_THREAD = 5000;             // Number of events per thread
const int CACHE_SIZE       = 50 * 1024 * 1024; // 50 MB cache size

struct Histo1DCommand {
  TString name;       // Histogram name
  TString title;      // Histogram title
  TString expression; // Expression to fill the histogram
  int bins;           // Number of bins
  double xMin;        // Minimum x-axis value
  double xMax;        // Maximum x-axis value
  TString cut;        // Cut expression
};
struct Histo2DCommand {
  TString name;         // Histogram name
  TString title;        // Histogram title
  TString expression_x; // Expression for x-axis
  TString expression_y; // Expression for y-axis
  int binsX;            // Number of bins in x-axis
  double xMin;          // Minimum x-axis value
  double xMax;          // Maximum x-axis value
  int binsY;            // Number of bins in y-axis
  double yMin;          // Minimum y-axis value
  double yMax;          // Maximum y-axis value
  TString cut;          // Cut expression
};

// Define a structure to store cut information
struct CutInfo {
  std::string var1; // First variable or branch name
  std::string op;   // Operator (e.g., ==, !=, >, <, >=, <=)
  std::string var2; // Second variable or branch name or constant value
};

template <typename T> void add_branch(TTree *tree, const char *branch_name, T *branch_data) {
  // Add a branch to the tree
  tree->SetBranchStatus(branch_name, 1); // Enable the branch
  tree->SetBranchAddress(branch_name, branch_data);
  tree->AddBranchToCache(branch_name, kTRUE);
  return;
}

bool getValues(const TString &var, std::vector<double> &values,
               const std::map<std::string, std::vector<Double_t>> &branchDataDouble,
               const std::map<std::string, std::vector<Int_t>> &branchDataInt, int nEntries_Evt) {
  if (nEntries_Evt < 0)
    return false;
  if (branchDataDouble.find(var.Data()) != branchDataDouble.end()) {
    values = branchDataDouble.at(var.Data());
  } else if (branchDataInt.find(var.Data()) != branchDataInt.end()) {
    for (const auto &val : branchDataInt.at(var.Data())) {
      values.push_back(static_cast<double>(val));
    }
  } else if (var.IsFloat() || var.IsDigit()) {
    values.insert(values.end(), nEntries_Evt, var.Atof());
  } else {
    std::cerr << "Error: Could not find or convert " << var << std::endl;
    return false;
  }
  return true;
}
// Function to parse a cut expression into CutInfo objects
std::vector<CutInfo> parseCutExpression(const TString &cutExpression) {
  std::vector<CutInfo> cuts;
  std::string cutExpr = cutExpression.Data();
  std::vector<std::string> tokens;
  tokens.reserve(10); // Reserve space for tokens to avoid multiple allocations
  size_t start = 0, end = 0;

  // Tokenize the cut expression based on "&&" or "||"
  while ((end = cutExpr.find_first_of("&|", start)) != std::string::npos) {
    if (end > start) {
      tokens.push_back(cutExpr.substr(start, end - start));
    }
    start = cutExpr.find_first_not_of("&|", end);
  }
  if (start < cutExpr.size()) {
    tokens.push_back(cutExpr.substr(start));
  }

  for (const auto &token : tokens) {
    std::string trimmedToken = token;
    trimmedToken.erase(remove(trimmedToken.begin(), trimmedToken.end(), ' '), trimmedToken.end()); // Remove spaces

    // Parse the token into var1, operator, and var2
    std::string var1, op, var2;
    size_t opPos = trimmedToken.find("==");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("!=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find(">=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("<=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find(">");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("<");

    if (opPos != std::string::npos) {
      var1 = trimmedToken.substr(0, opPos);
      if (trimmedToken.substr(opPos, 2) == "==" || trimmedToken.substr(opPos, 2) == "!=" ||
          trimmedToken.substr(opPos, 2) == ">=" || trimmedToken.substr(opPos, 2) == "<=") {
        op = trimmedToken.substr(opPos, 2);
      } else {
        op = trimmedToken.substr(opPos, 1); // Handle single-character operators like > or <
      }
      var2 = trimmedToken.substr(opPos + op.length());
      CutInfo cut;
      cut.var1 = var1;
      cut.op   = op;
      cut.var2 = var2;
      cuts.push_back(cut);
    } else {
      std::cerr << "Error: Invalid cut token format: " << trimmedToken << std::endl;
    }
  }
  return cuts;
}

void processCutExpression(const TString &cutExpression, const std::function<void(const TString &)> &processBranch) {
  if (!cutExpression.IsNull()) {
    auto cutInfos = parseCutExpression(cutExpression);
    for (const auto &cut : cutInfos) {
      processBranch(cut.var1.c_str());
      processBranch(cut.var2.c_str());
    }
  }
}

// Function to process a range of events in a thread
void processEventsRange(TString inputFileName, TString treeName, Long64_t start, Long64_t end,
                        const std::vector<Histo1DCommand> &hist1DCommands,
                        const std::vector<Histo2DCommand> &hist2DCommands, std::map<TString, TH1D *> &localHistograms1D,
                        std::map<TString, TH2D *> &localHistograms2D) {

  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error: Could not open file " << inputFileName << std::endl;
    return;
  }

  TTree *clonedTree = (TTree *)inputFile->Get(treeName);
  if (!clonedTree) {
    std::cerr << "Error: Could not find tree " << treeName << " in file " << inputFileName << std::endl;
    inputFile->Close();
    return;
  }
  clonedTree->SetCacheSize(CACHE_SIZE);
  clonedTree->SetBranchStatus("*", 0); // Disable all branches initially
  // Local branch data for this thread
  std::map<std::string, std::vector<Double_t>> branchDataDouble;
  std::map<std::string, std::vector<Int_t>> branchDataInt;
  std::map<std::string, int> nEntriesMap;

  // Loop through histogram expressions and add each tree name to branchDataDouble or branchDataInt
  auto processBranch = [&](const TString &branchName) {
    if (branchDataDouble.find(branchName.Data()) == branchDataDouble.end() &&
        branchDataInt.find(branchName.Data()) == branchDataInt.end()) {
      TBranch *branch = clonedTree->GetBranch(branchName.Data());
      if (branch) {
        TLeaf *leaf = branch->GetLeaf(branchName.Data());
        if (leaf) {
          TString leafType = leaf->GetTypeName();
          if (leafType == "Double_t" || leafType == "Float_t") {
            branchDataDouble[branchName.Data()] = std::vector<Double_t>(NDDATA_MAX);
          } else if (leafType == "Int_t") {
            branchDataInt[branchName.Data()] = std::vector<Int_t>(NDDATA_MAX);
          } else {
            std::cerr << "Error: Unsupported leaf type " << leafType << " for branch " << branchName << std::endl;
          }
        }
        // delete leaf;
      }
      // delete branch;
    }
  };

  // Process 1D histogram expressions and cuts
  for (const auto &cmd : hist1DCommands) {
    if (!cmd.expression.IsNull()) {
      processBranch(cmd.expression);
    }
    processCutExpression(cmd.cut, processBranch);
  }

  // Process 2D histogram expressions and cuts
  for (const auto &cmd : hist2DCommands) {
    if (!cmd.expression_x.IsNull()) {
      processBranch(cmd.expression_x);
    }
    if (!cmd.expression_y.IsNull()) {
      processBranch(cmd.expression_y);
    }
    processCutExpression(cmd.cut, processBranch);
  }

  for (auto &branchPair : branchDataDouble) {
    TString nEntriesBranchName = "Ndata." + TString(branchPair.first.c_str());
    TBranch *nEntriesBranch    = clonedTree->GetBranch(nEntriesBranchName.Data());
    if (nEntriesBranch) {
      nEntriesMap[branchPair.first] = 0; // Initialize nEntries for this branch
      add_branch(clonedTree, nEntriesBranchName.Data(), &nEntriesMap[branchPair.first]);
    } else {
      nEntriesMap[branchPair.first] = 1;
    }
  }

  for (auto &branchPair : branchDataInt) {
    TString nEntriesBranchName = "Ndata." + TString(branchPair.first.c_str());
    TBranch *nEntriesBranch    = clonedTree->GetBranch(nEntriesBranchName.Data());
    if (nEntriesBranch) {
      nEntriesMap[branchPair.first] = 0; // Initialize nEntries for this branch
      add_branch(clonedTree, nEntriesBranchName.Data(), &nEntriesMap[branchPair.first]);
    }
  }

  // Thread-safe access to the tree
  // auto getEntryThreadSafe = [&](Long64_t entry) {
  //   std::lock_guard<std::mutex> lock(treeMutex);
  //   tree->GetEntry(entry);
  // };
  // Set up branch addresses (similar to the main function)
  for (auto &branchPair : branchDataDouble) {
    add_branch(clonedTree, branchPair.first.c_str(), branchPair.second.data());
  }
  for (auto &branchPair : branchDataInt) {
    add_branch(clonedTree, branchPair.first.c_str(), branchPair.second.data());
  }

  // Process events in the range
  for (Long64_t i = start; i < end; ++i) {

    // Clear branch data for this event
    for (auto &branchPair : branchDataDouble) {
      branchPair.second.clear();
      branchPair.second.resize(NDDATA_MAX, 0);
    }
    for (auto &branchPair : branchDataInt) {
      branchPair.second.clear();
      branchPair.second.resize(NDDATA_MAX, 0);
    }
    // Safely get the entry
    // getEntryThreadSafe(i);

    // treeMutex.lock();
    clonedTree->GetEntry(i);
    // treeMutex.unlock();

    // clonedTree->GetEntry(i);
    // for (const auto &value : branchDataDouble["L.gem.clust.axis"]) {
    //   std::cout << "Value: " << value << std::endl;
    // }

    // Ensure NData values do not exceed their respective limits
    for (auto &branchPair : nEntriesMap) {
      if (branchPair.second > NDDATA_MAX) {
        branchPair.second = NDDATA_MAX;
        std::cout << "Warning: NData for " << branchPair.first << " exceeds maximum limit. Clamping to " << NDDATA_MAX
                  << std::endl;
      }
    }
    // Process 1D histograms
    for (const auto &cmd : hist1DCommands) {
      if (localHistograms1D.find(cmd.name) == localHistograms1D.end()) {
        localHistograms1D[cmd.name] = new TH1D(cmd.name, cmd.title, cmd.bins, cmd.xMin, cmd.xMax);
      }
      TH1D *hist = localHistograms1D[cmd.name];
      // hist->SetDirectory(nullptr); // Prevent ROOT from adding it to the global directory
      // Parse the cut expression into CutInfo objects
      std::vector<CutInfo> cutInfos;
      if (!cmd.cut.IsNull()) {
        cutInfos = parseCutExpression(cmd.cut);
      }
      // Check if the cut expression is valid
      std::vector<bool> passCut(nEntriesMap[cmd.expression.Data()], true);
      for (const auto &cut : cutInfos) {
        std::vector<double> values1, values2;

        // Get the values of var1
        if (!getValues(cut.var1, values1, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
          passCut = {false};
          break;
        }

        // Get the values of var2
        if (!getValues(cut.var2, values2, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
          passCut = {false};
          break;
        }

        // Compare values based on the operator
        for (size_t j = 0; j < nEntriesMap[cmd.expression.Data()]; ++j) {
          if (cut.op == "==") {
            passCut[j] = (values1[j] == values2[j]);
          } else if (cut.op == "!=") {
            passCut[j] = (values1[j] != values2[j]);
          } else if (cut.op == ">") {
            passCut[j] = (values1[j] > values2[j]);
          } else if (cut.op == "<") {
            passCut[j] = (values1[j] < values2[j]);
          } else if (cut.op == ">=") {
            passCut[j] = (values1[j] >= values2[j]);
          } else if (cut.op == "<=") {
            passCut[j] = (values1[j] <= values2[j]);
          } else {
            std::cerr << "Error: Unknown operator " << cut.op << std::endl;
            passCut[j] = false;
            break;
          }
        }
      }
      std::vector<double> valueVec;
      if (!getValues(cmd.expression, valueVec, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression << std::endl;
        continue;
      }
      for (size_t j = 0; j < nEntriesMap[cmd.expression.Data()]; ++j) {
        if (passCut[j]) {
          hist->Fill(valueVec[j]);
        }
      }
    }

    // Process 2D histograms
    for (const auto &cmd : hist2DCommands) {
      if (localHistograms2D.find(cmd.name) == localHistograms2D.end()) {
        localHistograms2D[cmd.name] =
            new TH2D(cmd.name, cmd.title, cmd.binsX, cmd.xMin, cmd.xMax, cmd.binsY, cmd.yMin, cmd.yMax);
      }
      TH2D *hist2D = localHistograms2D[cmd.name];
      hist2D->SetDirectory(nullptr); // Prevent ROOT from adding it to the global directory
      // Parse the cut expression into CutInfo objects
      std::vector<CutInfo> cutInfos;
      if (!cmd.cut.IsNull()) {
        cutInfos = parseCutExpression(cmd.cut);
      }

      // Check if the cut expression is valid
      std::vector<bool> passCut(nEntriesMap[cmd.expression_x.Data()], true);
      for (const auto &cut : cutInfos) {
        std::vector<double> values1, values2;

        // Get the values of var1
        if (!getValues(cut.var1, values1, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression_x.Data()])) {
          passCut = {false};
          break;
        }

        // Get the values of var2
        if (!getValues(cut.var2, values2, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression_x.Data()])) {
          passCut = {false};
          break;
        }

        // Compare values based on the operator
        for (size_t j = 0; j < nEntriesMap[cmd.expression_x.Data()]; ++j) {
          if (cut.op == "==") {
            passCut[j] = (values1[j] == values2[j]);
          } else if (cut.op == "!=") {
            passCut[j] = (values1[j] != values2[j]);
          } else if (cut.op == ">") {
            passCut[j] = (values1[j] > values2[j]);
          } else if (cut.op == "<") {
            passCut[j] = (values1[j] < values2[j]);
          } else if (cut.op == ">=") {
            passCut[j] = (values1[j] >= values2[j]);
          } else if (cut.op == "<=") {
            passCut[j] = (values1[j] <= values2[j]);
          } else {
            std::cerr << "Error: Unknown operator " << cut.op << std::endl;
            passCut[j] = false;
            break;
          }
        }
      }

      std::vector<double> valueVecX, valueVecY;
      if (!getValues(cmd.expression_x, valueVecX, branchDataDouble, branchDataInt,
                     nEntriesMap[cmd.expression_x.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression_x << std::endl;
        continue;
      }
      if (!getValues(cmd.expression_y, valueVecY, branchDataDouble, branchDataInt,
                     nEntriesMap[cmd.expression_y.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression_y << std::endl;
        continue;
      }
      for (size_t j = 0; j < nEntriesMap[cmd.expression_x.Data()]; ++j) {
        if (passCut[j]) {
          hist2D->Fill(valueVecX[j], valueVecY[j]);
        }
      }
    }
    if (start == 0 && i % 10 == 0) {
      double progress = ((double)(i - start) / (end - start)) * 100.0;
      std::cout << "\rThread progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }
  }
  if (start == 0) {
    std::cout << "\rThread progress: 100%" << std::endl;
  }
}
// void gem_histos(const char *inputFileName, const char *treeName) {
void lad_histos_MT(
    const char *inputFileName = "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1.root",
    bool lad_only = true, int num_evts = -1) {

  auto start_time = std::chrono::high_resolution_clock::now();
  // // Set no maximum limit on the stack size
  // struct rlimit rl;
  // if (getrlimit(RLIMIT_STACK, &rl) == 0) {
  //   rl.rlim_cur = RLIM_INFINITY;
  //   setrlimit(RLIMIT_STACK, &rl);
  // }
  // Enable thread safety for ROOT
  ROOT::EnableThreadSafety();
  const char *treeName = "T";
  // Open the ROOT file and get the tree
  TFile *inputFile = TFile::Open(inputFileName, "UPDATE");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error: Could not open file " << inputFileName << std::endl;
    return;
  }

  TTree *tree = (TTree *)inputFile->Get(treeName);
  if (!tree) {
    std::cerr << "Error: Could not find tree " << treeName << " in file " << inputFileName << std::endl;
    inputFile->Close();
    return;
  }

  // Define the list of histogram commands
  std::vector<TString> prefixes;
  if (lad_only) {
    prefixes.push_back("L");
  } else {
    prefixes.push_back("H");
    prefixes.push_back("P");
  }
  std::vector<Histo1DCommand> hist1DCommands;
  std::vector<Histo2DCommand> hist2DCommands;

  for (const auto &prefix : prefixes) {
    // Add 1D histogram commands for each prefix
    hist1DCommands.insert(
        hist1DCommands.end(),
        {{prefix + "_h1_gem_clustWidthU_0", "X cluster size Layer 0; cluster size", prefix + ".gem.clust.nstrip", 10,
          0.5, 10.5, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustWidthV_0", "Y cluster size Layer 0; cluster size", prefix + ".gem.clust.nstrip", 10,
          0.5, 10.5, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustWidthU_1", "X cluster size Layer 1; cluster size", prefix + ".gem.clust.nstrip", 10,
          0.5, 10.5, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustWidthV_1", "Y cluster size Layer 1; cluster size", prefix + ".gem.clust.nstrip", 10,
          0.5, 10.5, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustSampMaxU_0", "Peak time sample Layer 0; Peak time sample (X strip)",
          prefix + ".gem.clust.maxsamp", 6, -0.5, 5.5, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustSampMaxV_0", "Peak time sample Layer 0; Peak time sample (Y strip)",
          prefix + ".gem.clust.maxsamp", 6, -0.5, 5.5, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustSampMaxU_1", "Peak time sample Layer 1; Peak time sample (X strip)",
          prefix + ".gem.clust.maxsamp", 6, -0.5, 5.5, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustSampMaxV_1", "Peak time sample Layer 1; Peak time sample (Y strip)",
          prefix + ".gem.clust.maxsamp", 6, -0.5, 5.5, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustADCMaxU_0", "X cluster Max strip ADC Layer 0; MAX strip ADC",
          prefix + ".gem.clust.maxadc", 1500, 0, 1500, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustADCMaxV_0", "Y cluster Max strip ADC Layer 0; MAX strip ADC",
          prefix + ".gem.clust.maxadc", 1500, 0, 1500, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustADCMaxU_1", "X cluster Max strip ADC Layer 1; MAX strip ADC",
          prefix + ".gem.clust.maxadc", 1500, 0, 1500, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustADCMaxV_1", "Y cluster Max strip ADC Layer 1; MAX strip ADC",
          prefix + ".gem.clust.maxadc", 1500, 0, 1500, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustADCSumU_0", "X cluster ADC sum Layer 0; ADC sum", prefix + ".gem.clust.adc", 1500, 0,
          3000, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustADCSumV_0", "Y cluster ADC sum Layer 0; ADC sum", prefix + ".gem.clust.adc", 1500, 0,
          3000, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer<1"},
         {prefix + "_h1_gem_clustADCSumU_1", "X cluster ADC sum Layer 1; ADC sum", prefix + ".gem.clust.adc", 1500, 0,
          3000, prefix + ".gem.clust.axis<1&&" + prefix + ".gem.clust.layer>0"},
         {prefix + "_h1_gem_clustADCSumV_1", "Y cluster ADC sum Layer 1; ADC sum", prefix + ".gem.clust.adc", 1500, 0,
          3000, prefix + ".gem.clust.axis>0&&" + prefix + ".gem.clust.layer>0"},
         //  {prefix + "_h1_gem_nhits_0", "Number of hits Layer 0; Number of hits", prefix + ".gem.sp.nhits", 50,
         //  -0.5, 49.5,
         //   prefix + ".gem.sp.layer<1"},
         //  {prefix + "_h1_gem_nhits_1", "Number of hits Layer 1; Number of hits", prefix + ".gem.sp.nhits", 50,
         //  -0.5, 49.5,
         //   prefix + ".gem.sp.layer>0"},
         {prefix + "_h1_gem_time_0", "Time Mean; Time Mean", prefix + ".gem.sp.time", 100, 0, 120,
          prefix + ".gem.sp.layer<1"},
         {prefix + "_h1_gem_time_1", "Time Mean; Time Mean", prefix + ".gem.sp.time", 100, 0, 120,
          prefix + ".gem.sp.layer>0"},
         {prefix + "_h1_gem_ADCMean_0", "ADC Mean; ADC Mean", prefix + ".gem.sp.adc", 100, 0, 70000,
          prefix + ".gem.sp.layer<1"},
         {prefix + "_h1_gem_ADCMean_1", "ADC Mean; ADC Mean", prefix + ".gem.sp.adc", 100, 0, 70000,
          prefix + ".gem.sp.layer>0"},
         {prefix + "_h1_gem_ADCAsym_0", "ADC Asymmetry; ADC Asymmetry", prefix + ".gem.sp.asym", 100, -1, 1,
          prefix + ".gem.sp.layer<1"},
         {prefix + "_h1_gem_ADCAsym_1", "ADC Asymmetry; ADC Asymmetry", prefix + ".gem.sp.asym", 100, -1, 1,
          prefix + ".gem.sp.layer>0"},
         {prefix + "_h1_gem_TimeDiff_0", "Time Difference; Time Difference", prefix + ".gem.sp.dt", 100, -40, 40,
          prefix + ".gem.sp.layer<1"},
         {prefix + "_h1_gem_TimeDiff_1", "Time Difference; Time Difference", prefix + ".gem.sp.dt", 100, -40, 40,
          prefix + ".gem.sp.layer>0"},
         {prefix + "_h1_gem_TimeCorr_0", "Corrected Time; Corrected Time", prefix + ".gem.sp.ct", 100, -40, 40,
          prefix + ".gem.sp.layer<1"},
         {prefix + "_h1_gem_TimeCorr_1", "Corrected Time; Corrected Time", prefix + ".gem.sp.ct", 100, -40, 40,
          prefix + ".gem.sp.layer>0"}});

    // Add 2D histogram commands for each prefix
    TString prefix_lower = prefix;
    prefix_lower.ToLower();
    hist2DCommands.insert(
        hist2DCommands.end(),
        {{prefix_lower + "ladhodwTrack_000_good_hit_time_avg", "LAD 000 Good Hit with Track Edep; PMT Number; Hit Edep",
          prefix + ".ladhod.goodhit_paddle_0", prefix + ".ladhod.goodhit_hitedep_0", 11, 0.5, 11.5, 200, 0, 500000,
          prefix + ".ladhod.goodhit_plane_0==0"},
         {prefix_lower + "ladhodwTrack_001_good_hit_time_avg", "LAD 001 Good Hit with Track Edep; PMT Number; Hit Edep",
          prefix + ".ladhod.goodhit_paddle_0", prefix + ".ladhod.goodhit_hitedep_0", 11, 0.5, 11.5, 200, 0, 500000,
          prefix + ".ladhod.goodhit_plane_0==1"},
         {prefix_lower + "ladhodwTrack_100_good_hit_time_avg", "LAD 100 Good Hit with Track Edep; PMT Number; Hit Edep",
          prefix + ".ladhod.goodhit_paddle_0", prefix + ".ladhod.goodhit_hitedep_0", 11, 0.5, 11.5, 200, 0, 500000,
          prefix + ".ladhod.goodhit_plane_0==2"},
         {prefix_lower + "ladhodwTrack_101_good_hit_time_avg", "LAD 101 Good Hit with Track Edep; PMT Number; Hit Edep",
          prefix + ".ladhod.goodhit_paddle_0", prefix + ".ladhod.goodhit_hitedep_0", 11, 0.5, 11.5, 200, 0, 500000,
          prefix + ".ladhod.goodhit_plane_0==3"},
         {prefix_lower + "ladhodwTrack_200_good_hit_time_avg", "LAD 200 Good Hit with Track Edep; PMT Number; Hit Edep",
          prefix + ".ladhod.goodhit_paddle_0", prefix + ".ladhod.goodhit_hitedep_0", 11, 0.5, 11.5, 200, 0, 500000,
          prefix + ".ladhod.goodhit_plane_0==4"},
         {prefix + "_h2_gem_2dhit_0", "Layer 0 XY space points (2Dhits); x; y", prefix + ".gem.sp.posX",
          prefix + ".gem.sp.posY", 100, -65, 65, 100, -35, 35, prefix + ".gem.sp.layer<1"},
         {prefix + "_h2_gem_2dhit_1", "Layer 1 XY space points (2Dhits); x; y", prefix + ".gem.sp.posX",
          prefix + ".gem.sp.posY", 100, -65, 65, 100, -35, 35, prefix + ".gem.sp.layer>0"}});
  }
  int nHistograms = hist1DCommands.size() + hist2DCommands.size();
  int histCount   = 0; // Counter for histograms

  // Loop over the histogram commands and create histograms
  Long64_t nEntries = tree->GetEntries();
  if (num_evts > 0 && num_evts < nEntries) {
    nEntries = num_evts;
  }

  // Declare the map to store 1D histograms
  static std::map<TString, TH1D *> histograms1D;
  static std::map<TString, TH2D *> histograms2D;

  // Divide events among threads
  int numThreads_max = std::thread::hardware_concurrency();
  int numThreads     = nEntries / NEVTS_PER_THREAD;
  numThreads         = std::min(numThreads, numThreads_max);
  numThreads         = std::min(numThreads, MAX_THREADS);
  numThreads         = std::max(numThreads, 1); // Ensure at least one thread

  std::cout << "Running on " << numThreads << " cores. " << numThreads_max << " available." << std::endl;
  Long64_t chunkSize = nEntries / numThreads;
  std::vector<std::thread> threads;
  std::vector<std::map<TString, TH1D *>> threadHistograms1D(numThreads);
  std::vector<std::map<TString, TH2D *>> threadHistograms2D(numThreads);

  for (int t = 0; t < numThreads; ++t) {
    Long64_t start = t * chunkSize;
    Long64_t end   = (t == numThreads - 1) ? nEntries : start + chunkSize;

    // TTree *clonedTree = tree->CloneTree(0);
    threads.emplace_back([&, t, start, end]() {
      processEventsRange(inputFileName, treeName, start, end, hist1DCommands, hist2DCommands, threadHistograms1D[t],
                         threadHistograms2D[t]);
      // Ensure unique histogram names
      for (auto &pair : threadHistograms1D[t]) {
        pair.second->SetDirectory(nullptr); // Prevent ROOT from adding to the global directory
      }
      for (auto &pair : threadHistograms2D[t]) {
        pair.second->SetDirectory(nullptr); // Prevent ROOT from adding to the global directory
      }
    });
  }

  // Wait for all threads to finish
  for (auto &thread : threads) {
    thread.join();
  }

  // Merge histograms from all threads
  std::map<TString, TH1D *> finalHistograms1D;
  std::map<TString, TH2D *> finalHistograms2D;

  for (int t = 0; t < numThreads; ++t) {
    for (auto &pair : threadHistograms1D[t]) {
      if (finalHistograms1D.find(pair.first) == finalHistograms1D.end()) {
        finalHistograms1D[pair.first] = (TH1D *)pair.second->Clone();
      } else {
        finalHistograms1D[pair.first]->Add(pair.second);
      }
      // delete pair.second;
      // std::cout << finalHistograms1D[pair.first]->Integral() << std::endl;
    }
    for (auto &pair : threadHistograms2D[t]) {
      if (finalHistograms2D.find(pair.first) == finalHistograms2D.end()) {
        finalHistograms2D[pair.first] = (TH2D *)pair.second->Clone();
      } else {
        finalHistograms2D[pair.first]->Add(pair.second);
      }
      // delete pair.second;
    }
  }

  inputFile->cd();
  // Write final histograms to file
  for (auto &pair : finalHistograms1D) {
    // pair.second->SetDirectory(inputFile);
    pair.second->Write();
    delete pair.second;
  }
  for (auto &pair : finalHistograms2D) {
    pair.second->Write();
    delete pair.second;
  }

  // inputFile->Write();
  inputFile->Close();

  auto end_time                              = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  std::cout << "Elapsed time: " << elapsed_time.count() << " seconds" << std::endl;
}