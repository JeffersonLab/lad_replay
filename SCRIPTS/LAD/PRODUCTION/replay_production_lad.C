// replays both the LAD hodoscope and GEM detectors (not standard spectrometers)
#include "../../LAD/LAD_link_defs.h" // Leave this line commented. Used for debugging purposes only.
#include "LADFilteredStreamBuf.h"

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

  gHcParms->AddString("lgem_pedfile", Form("PARAM/LAD/GEM/PED/gem_ped_%d.dat", ped_cm_file_num));
  gHcParms->AddString("lgem_cmfile", Form("PARAM/LAD/GEM/CM/CommonModeRange_%d.txt", ped_cm_file_num));
  return;
}

void replay_production_lad(Int_t RunNumber = 0, Int_t MaxEvent = 0, int run_type = 0) {

  // Get RunNumber and MaxEvent if not provided.
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
      exit;
    }
  }

  // Set the run type based on the run number
  const char *RunFileNamePattern;
  const char *ROOTFileNamePattern;
  TString ROOTFileName;
  ROOTFileNamePattern = "ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_%d_%d.root";

  switch (run_type) {
  case 0:
    RunFileNamePattern = "lad_Production_%02d.dat.0";
    break;
  case 1:
    RunFileNamePattern = "lad_Production_noGEM_%02d.dat.0";
    break;
  case 2:
    RunFileNamePattern = "lad_LADwGEMwROC2_%02d.dat.0";
    break;
  case 3:
    RunFileNamePattern = "lad_GEMonly_%02d.dat.0";
    break;
  case 4:
    RunFileNamePattern = "lad_LADonly_%02d.dat.0";
    break;
  case 5:
    RunFileNamePattern = "lad_SHMS_HMS_%02d.dat.0";
    break;
  case 6:
    RunFileNamePattern = "lad_SHMS_%02d.dat.0";
    break;
  case 7:
    RunFileNamePattern = "lad_HMS_%02d.dat.0";
    break;
  case 8:
    RunFileNamePattern = "lad_LADwROC2_%02d.dat.0";
    break;
  case 9:
    RunFileNamePattern = "shms_cosmic_%02d.dat.0";
    break;
  default:
    cout << "Invalid run type: " << run_type << ". Please enter a valid run type." << endl;
    return;
    break;
  }
  ROOTFileName = Form(ROOTFileNamePattern, RunNumber, MaxEvent);

  vector<TString> pathList;
  pathList.push_back(".");
  // pathList.push_back("./ROOTfiles/COSMICS/raw/");
  pathList.push_back("/cache/hallc/c-lad/raw/");
  pathList.push_back("/volatile/hallc/c-lad/ehingerl/raw_data/LAD_cosmic");

  // const char *ROOTFileNamePattern = "ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_%d_%d.root";
  // const char *ROOTFileNamePattern = "ROOTfiles/COSMICS/LAD_wREF_cosmic_hall_%d_%d.root";

  // Load Global parameters
  // Add variables to global list.
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  gHcParms->Load("PARAM/TRIG/tshms.param");

  // Load correct GEM common mode and pedestal files
  load_GEM_CM_PED(RunNumber);
  // Load fadc debug parameters
  // gHcParms->Load("PARAM/HMS/GEN/h_fadc_debug.param");

  // const char* CurrentFileNamePattern = "low_curr_bcm/bcmcurrent_%d.param";
  // gHcParms->Load(Form(CurrentFileNamePattern, RunNumber));

  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();
  if (RunNumber < 22590) {
    gHcDetectorMap->Load("MAPS/LAD/DETEC/HODO/lhodo_laser.map");
  }
  else{
    gHcDetectorMap->Load("MAPS/LAD/DETEC/HODO/lhodo_laser_May14.map");
  }

  THaApparatus *TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);
  THcTrigDet *shms = new THcTrigDet("shms", "SHMS Trigger Information");
  TRG->AddDetector(shms);

  // Add LAD detector
  THcLADSpectrometer *LAD = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(LAD);

  THcLADHodoscope *lhod = new THcLADHodoscope("ladhod", "LAD Hodoscope");
  LAD->AddDetector(lhod);

  THcLADGEM *gem = new THcLADGEM("gem", "gem");
  LAD->AddDetector(gem);

  // // Add handler for prestart event 125.
  // THcConfigEvtHandler *ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  // gHaEvtHandlers->Add(ev125);
  // // Add handler for EPICS events
  // THaEpicsEvtHandler *hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 181");
  // gHaEvtHandlers->Add(hcepics);
  // // Add handler for scaler events
  // THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 2");
  // hscaler->AddEvtType(2);
  // hscaler->AddEvtType(129);
  // hscaler->SetDelayedType(129);
  // hscaler->SetUseFirstEvent(kTRUE);
  // gHaEvtHandlers->Add(hscaler);

  // Add event handler for DAQ configuration event
  // THcConfigEvtHandler *hconfig = new THcConfigEvtHandler("hconfig", "Hall C configuration event handler");
  // gHaEvtHandlers->Add(hconfig);

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Acpparatus's and PhysicsModules,
  // and executes the output routines.
  THcAnalyzer *analyzer = new THcAnalyzer;

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent *event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THcRun *run = new THcRun(pathList, Form(RunFileNamePattern, RunNumber));

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");

  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();

  // Define the analysis parameters
  analyzer->SetCountMode(2); // 0 = counter is # of physics triggers
                             // 1 = counter is # of all decode reads
                             // 2 = counter is event number

  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(181);
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap_lad.dat");
  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define output DEF-file
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/lstackana_production_all.def");
  // Define cuts file
  // analyzer->SetCutFile("DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts_cosmic.def"); // optional
  analyzer->SetCutFile("DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts.def");
  // File to record cuts accounting information for cuts
  // analyzer->SetSummaryFile(
  //     Form("REPORT_OUTPUT/HMS/PRODUCTION/summary_production_%d_%d.report", RunNumber, MaxEvent)); // optional
  // Start the actual analysis.

  // Add a buffer filter to handle specific errors
  LADFilteredStreamBuf ladbuf(std::clog);
  ladbuf.includeRootErrorMessages(true); // Include ROOT error messages
  ladbuf.addFilterString("Module L.react does not exist");
  // ladbuf.addFilterString("crate 6");

  // Start the actual analysis
  analyzer->Process(run);
  // Create report file from template.
  // analyzer->PrintReport("TEMPLATES/LAD/PRODUCTION/lstackana_production.template",
  //                       Form("REPORT_OUTPUT/LAD/PRODUCTION/replay_lad_production_%d_%d.report", RunNumber,
  //                       MaxEvent));
}
