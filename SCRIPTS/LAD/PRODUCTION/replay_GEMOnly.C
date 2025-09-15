// Only replay's LAD GEM detector
// Currently set up to replay GEM data from the test lab
// #include "../../LAD/LAD_link_defs.h". Leave this line commented. Used for debugging purposes.

#include "LADFilteredStreamBuf.h"
#include "MultiFileRun.h"
#include "THcLADKine.h" // Include the header for THcLADKine
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "THcRun.h"
#include "THaRun.h"
void load_GEM_CM_PED(int runNumber) {
  std::vector<int> ped_cm_runs;
  int ped_cm_runs_count = 0;

  // Open the file
  std::ifstream infile("PARAM/LAD/GEM/lgem_cm_ped_runs.param");
  if (infile) {
    std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    std::stringstream ss(content);
    std::string token;
    while (std::getline(ss, token, ',')) {
      std::stringstream token_ss(token);
      int value;
      if (token_ss >> value) {
        ped_cm_runs.push_back(value);
      }
    }
    ped_cm_runs_count = ped_cm_runs.size();
  } else {
    std::cerr << "Error: Could not open PARAM/LAD/GEM/lgem_cm_ped_runs.param" << std::endl;
    ped_cm_runs_count = 0;
  }
  std::sort(ped_cm_runs.begin(), ped_cm_runs.end());

  int ped_cm_file_num = ped_cm_runs[ped_cm_runs_count - 1]; // Default to the last run number in the list
  for (int i = ped_cm_runs_count - 1; i > 0; --i) {
    if (ped_cm_runs[i] < runNumber) {
      ped_cm_file_num = ped_cm_runs[i];
      break;
    }
  }

  gHcParms->AddString("lgem_pedfile", Form("PARAM/LAD/GEM/PED/daq_ped_L_gem_run%d.dat", ped_cm_file_num));
  //gHcParms->AddString("lgem_pedfile", Form("PARAM/LAD/GEM/PED/gem_ped_%d.dat", ped_cm_file_num));
  gHcParms->AddString("lgem_cmfile", Form("PARAM/LAD/GEM/CM/db_cmr_L_gem_run%d.dat", ped_cm_file_num));
  //gHcParms->AddString("lgem_cmfile", Form("PARAM/LAD/GEM/CM/CommonModeRange_%d.txt", ped_cm_file_num));
  return;
}

void replay_GEMOnly(int RunNumber = 22820, int MaxEvent = 5000, int run_type = 3, int FirstEvent = 1,
                                int MaxSegment = 0, int FirstSegment = 0) {

  if (RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if (RunNumber <= 0)
      return;
  }
  if (MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if (MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      return;
    }
  }

  // Create file name patterns.
  // fname_prefix, RunNumber, iseg.
  // To replay files with no segment number, use -1 for MaxSegment.
  const char *RunFileNamePattern;
  const char *SummaryFileNamePattern;
  const char *REPORTFileNamePattern;

  if (MaxSegment == -1) {
    RunFileNamePattern = "%s_%d.dat";
  } else {
    RunFileNamePattern = "%s_%d.dat.%u";
  }
  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back("./raw/");
  pathList.push_back("./raw/../raw.copiedtotape");
  pathList.push_back("./cache");
  pathList.push_back("/cache/hallc/c-lad/raw/");

  const char *ROOTFileNamePattern;
  TString ROOTFileName;
  pathList.push_back("/volatile/hallc/c-lad/ehingerl/raw_data/LAD_cosmic");

  const char *fname_prefix;
  switch (run_type) {
  case 0:
    fname_prefix = "lad_Production";
    break;
  case 1:
    fname_prefix = "lad_Production_noGEM";
    break;
  case 2:
    fname_prefix = "lad_LADwGEMwROC2";
    break;
  case 3:
    fname_prefix = "lad_GEMonly";
    break;
  case 4:
    fname_prefix = "lad_LADonly";
    break;
  case 5:
    fname_prefix = "lad_SHMS_HMS";
    break;
  case 6:
    fname_prefix = "lad_SHMS";
    break;
  case 7:
    fname_prefix = "lad_HMS";
    break;
  default:
    cout << "Invalid run type: " << run_type << ". Please enter a valid run type." << endl;
    return;
    break;
  }



  ROOTFileNamePattern = "ROOTfiles/lad_gem_%d_%d.root";
  ROOTFileName = Form(ROOTFileNamePattern, RunNumber, MaxEvent);

  // Params
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // gHcParms->Load("PARAM/LAD/GEM/lgem_geom.param");
  // gHcParms->Load("DB_LAD/lgem_chan.map");
  // gHcParms->Load("DB_LAD/lgem_cuts.param");

  // Load correct GEM common mode and pedestal files
  load_GEM_CM_PED(RunNumber);
  // Detector map
  //  gHcDetectorMap = new THcDetectorMap();
  //  gHcDetectorMap->Load("detector.map");

  THcLADSpectrometer *lad = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(lad);

  THcLADGEM *gem = new THcLADGEM("gem", "gem");
  lad->AddDetector(gem);
  THcAnalyzer *analyzer = new THcAnalyzer;
  THaEvent *event       = new THaEvent;
  
  vector<string> fileNames = {};
  TString codafilename;
  if (MaxSegment == -1) {
    cout << RunFileNamePattern;
    codafilename.Form(RunFileNamePattern, fname_prefix, RunNumber);
    cout << "codafilename = " << codafilename << endl;
    fileNames.emplace_back(codafilename.Data());
  } else {
    for (Int_t iseg = FirstSegment; iseg <= MaxSegment; iseg++) {
      codafilename.Form(RunFileNamePattern, fname_prefix, RunNumber, iseg);
      cout << "codafilename = " << codafilename << endl;
      fileNames.emplace_back(codafilename.Data());
    }
  }
  THcRun* run = new THcRun(pathList, fileNames[0].c_str() );
  //THaRunBase *run = new THaRun(pathList, Form(RunFileNamePattern, RunNumber));
  //auto *run = new Podd::MultiFileRun(pathList, fileNames);

  run->SetEventRange(1, MaxEvent);
  // run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  run->SetDataRequired(0x7);
  run->Print();

  analyzer->SetEvent(event);
  analyzer->SetCountMode(2); // 2 = counter is event number


  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");
  analyzer->SetOutFile(ROOTFileName.Data());
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/lstackana_production_gem.def");
  analyzer->SetCutFile( "DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts_gem.def" );
  //  analyzer->SetCrateMapFileName("cratemap.dat");

  analyzer->Process(run);
}
