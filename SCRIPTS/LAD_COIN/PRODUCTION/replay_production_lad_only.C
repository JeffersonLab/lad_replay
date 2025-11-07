#include "LADFilteredStreamBuf.h"
#include "MultiFileRun.h"
#include "THcLADKine.h" // Include the header for THcLADKine
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../LAD/LAD_link_defs.h" //Leave this line commented. Used for debugging purposes only.

// Returns 1 if runNumber has a GEM CM PED file, 0 otherwise
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

void replay_production_lad_only(int RunNumber = 0, int MaxEvent = 0, int run_type = 1, int FirstEvent = 1,
                                int MaxSegment = 0, int FirstSegment = 0) {

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
    RunFileNamePattern = "%s_%05d.dat";
  } else {
    RunFileNamePattern = "%s_%05d.dat.%u";
  }
  vector<string> pathList;
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

  //(CA)
  TString REPORTFileName;
  REPORTFileNamePattern =
      "REPORT_OUTPUT/LAD_COIN/PRODUCTION/replayReport_LAD_coin_production_lad_only_%d_%d_%d_%d.report";
  REPORTFileName = Form(REPORTFileNamePattern, RunNumber, FirstSegment, MaxSegment, MaxEvent);

  ROOTFileNamePattern = "ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_lad_only_%d_%d_%d_%d.root";
  ROOTFileName        = Form(ROOTFileNamePattern, RunNumber, FirstSegment, MaxSegment, MaxEvent);
  // ROOTFileNamePattern = "ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_production_hall_%d_%d.root";
  // ROOTFileNamePattern = "ROOTfiles/LAD_COIN/CALIBRATION/LAD_COIN_calibration_hall_%d_%d.root";

  // LHE.
  // Consider changing naming scheme for root files. Also, consider adding summary file (in addition to report file).

  // Load global parameters
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD_COIN/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // Load parameters for HMS and SHMS trigger configuration
  gHcParms->Load("PARAM/TRIG/tshms.param");
  gHcParms->Load("PARAM/TRIG/thms.param");
  // Load fadc debug parameters
  // gHcParms->Load("PARAM/SHMS/GEN/p_fadc_debug.param");

  if (RunNumber > 22589)
    gHcParms->Load("PARAM/LAD/HODO/lhodo_cuts_May16.param");
  // Load the Hall C detector map
  // Load map depending on whether run is before or after SHMS DC swap
  gHcDetectorMap = new THcDetectorMap();
  if (RunNumber < 22157)
    gHcDetectorMap->Load("MAPS/LAD_COIN/DETEC/coin_lad_3pass_swapped_DC.map");
  else if (RunNumber < 22590)
    gHcDetectorMap->Load("MAPS/LAD_COIN/DETEC/coin_lad_5pass.map");
  else
    gHcDetectorMap->Load("MAPS/LAD_COIN/DETEC/coin_lad_5pass_May14.map");

  // load_GEM_CM_PED(RunNumber);
  // Add the dec data class for debugging
  // Podd::DecData *decData = new Podd::DecData("D", "Decoder Raw Data");
  // gHaApps->Add(decData);

  // Add trigger apparatus
  THaApparatus *TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);

  //////////////////////////////////////////////////////////////////////////
  //      SHMS
  //////////////////////////////////////////////////////////////////////////

  // Add trigger detector to trigger apparatus
  THcTrigDet *shms = new THcTrigDet("shms", "SHMS Trigger Information");
  // shms->SetSpectName("P");
  TRG->AddDetector(shms);

  // Set up the equipment to be analyzed
  THcHallCSpectrometer *SHMS = new THcHallCSpectrometer("P", "SHMS");
  SHMS->SetEvtType(1);
  SHMS->AddEvtType(4);
  SHMS->AddEvtType(5);
  SHMS->AddEvtType(6);
  SHMS->AddEvtType(7);
  gHaApps->Add(SHMS);
  // // Add Noble Gas Cherenkov to SHMS apparatus
  // THcCherenkov *ngcer = new THcCherenkov("ngcer", "Noble Gas Cherenkov");
  // SHMS->AddDetector(ngcer);
  // // Add drift chambers to SHMS apparatus
  // THcDC *pdc = new THcDC("dc", "Drift Chambers");
  // SHMS->AddDetector(pdc);
  // // Add hodoscope to SHMS apparatus
  // THcHodoscope *phod = new THcHodoscope("hod", "Hodoscope");
  // SHMS->AddDetector(phod);
  // // Add Heavy Gas Cherenkov to SHMS apparatus
  // THcCherenkov *phgcer = new THcCherenkov("hgcer", "Heavy Gas Cherenkov");
  // SHMS->AddDetector(phgcer);
  // // Add Aerogel Cherenkov to SHMS apparatus
  // THcAerogel *aero = new THcAerogel("aero", "Aerogel");
  // SHMS->AddDetector(aero);
  // // Add calorimeter to SHMS apparatus
  // THcShower *pcal = new THcShower("cal", "Calorimeter");
  // SHMS->AddDetector(pcal);

  THcLADHodoscope *shms_lhod = new THcLADHodoscope("ladhod", "LAD Hodoscope");
  SHMS->AddDetector(shms_lhod);

  // THcLADGEM *shms_gem = new THcLADGEM("gem", "gem");
  // shms_gem->setIgnoreVertex(true);
  // SHMS->AddDetector(shms_gem);

  // Add rastered beam apparatus
  THaApparatus *pbeam = new THcRasteredBeam("P.rb", "Rastered Beamline");
  gHaApps->Add(pbeam);
  // Add physics modules
  // Calculate reaction point
  // THcReactionPoint *prp = new THcReactionPoint("P.react", "SHMS reaction point", "P", "P.rb");
  // gHaPhysics->Add(prp);
  // // Calculate extended target corrections
  // THcExtTarCor *pext = new THcExtTarCor("P.extcor", "HMS extended target corrections", "P", "P.react");
  // gHaPhysics->Add(pext);
  // // Calculate golden track quantites
  // THaGoldenTrack *gtr = new THaGoldenTrack("P.gtr", "SHMS Golden Track", "P");
  // gHaPhysics->Add(gtr);
  // // Calculate primary (scattered beam - usually electrons) kinematics
  // THcPrimaryKine *pkin = new THcPrimaryKine("P.kin", "SHMS Single Arm Kinematics", "P", "P.rb");
  // gHaPhysics->Add(pkin);
  // // Calculate the hodoscope efficiencies
  // THcHodoEff *peff = new THcHodoEff("phodeff", "SHMS hodo efficiency", "P.hod");
  // gHaPhysics->Add(peff);

  // Can I leave this??
  // THcLADKine *ladkin_p = new THcLADKine("P.ladkin", "LAD Kinematics", "P", "P.kin", "P.react");
  // gHaPhysics->Add(ladkin_p);

  // Add event handler for scaler events
  THcScalerEvtHandler *pscaler = new THcScalerEvtHandler("P", "Hall C scaler event type 1");
  pscaler->AddEvtType(1);
  pscaler->AddEvtType(2);
  pscaler->AddEvtType(3);
  pscaler->AddEvtType(4);
  pscaler->AddEvtType(129);
  pscaler->SetDelayedType(129);
  pscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(pscaler);

  //////////////////////////////////////////////////////////////////////////
  //      HMS
  //////////////////////////////////////////////////////////////////////////
  THcTrigDet *hms = new THcTrigDet("hms", "HMS Trigger Information");
  // hms->SetSpectName("H");
  TRG->AddDetector(hms);

  THcHallCSpectrometer *HMS = new THcHallCSpectrometer("H", "HMS");
  HMS->SetEvtType(2);
  HMS->AddEvtType(4);
  HMS->AddEvtType(5);
  HMS->AddEvtType(6);
  HMS->AddEvtType(7);
  gHaApps->Add(HMS);
  // Add drift chambers to HMS apparatus
  // THcDC *hdc = new THcDC("dc", "Drift Chambers");
  // HMS->AddDetector(hdc);
  // // Add hodoscope to HMS apparatus
  // THcHodoscope *hhod = new THcHodoscope("hod", "Hodoscope");
  // HMS->AddDetector(hhod);
  // // Add Cherenkov to HMS apparatus
  // THcCherenkov *hcer = new THcCherenkov("cer", "Heavy Gas Cherenkov");
  // HMS->AddDetector(hcer);
  // Add Aerogel Cherenkov to HMS apparatus
  // THcAerogel* aero = new THcAerogel("aero", "Aerogel");
  // HMS->AddDetector(aero);
  // Add calorimeter to HMS apparatus
  // THcShower *hcal = new THcShower("cal", "Calorimeter");
  // HMS->AddDetector(hcal);

  THcLADHodoscope *hms_lhod = new THcLADHodoscope("ladhod", "LAD Hodoscope");
  HMS->AddDetector(hms_lhod);

  // THcLADGEM *hms_gem = new THcLADGEM("gem", "gem");
  // hms_gem->setIgnoreVertex(true);
  // HMS->AddDetector(hms_gem);

  // Add rastered beam apparatus
  // THaApparatus *hbeam = new THcRasteredBeam("H.rb", "Rastered Beamline");
  // gHaApps->Add(hbeam);
  // // Add physics modules
  // // Calculate reaction point
  // THcReactionPoint *hrp = new THcReactionPoint("H.react", "HMS reaction point", "H", "H.rb");
  // gHaPhysics->Add(hrp);
  // // Calculate extended target corrections
  // THcExtTarCor *hext = new THcExtTarCor("H.extcor", "HMS extended target corrections", "H", "H.react");
  // gHaPhysics->Add(hext);
  // // Calculate golden track quantities
  // THaGoldenTrack *hgtr = new THaGoldenTrack("H.gtr", "HMS Golden Track", "H");
  // gHaPhysics->Add(hgtr);

  // Calculate primary (scattered beam - usually electrons) kinematics
  // THcPrimaryKine *hkin = new THcPrimaryKine("H.kin", "HMS Single Arm Kinematics", "H", "H.rb");
  // gHaPhysics->Add(hkin);

  // // Calculate the hodoscope efficiencies
  // THcHodoEff *heff = new THcHodoEff("hhodeff", "HMS hodo efficiency", "H.hod");
  // gHaPhysics->Add(heff);

  // Can I leave this?
  // THcLADKine *ladkin_h = new THcLADKine("H.ladkin", "LAD Kinematics", "H", "H.kin", "H.react");
  // gHaPhysics->Add(ladkin_h);

  // Add event handler for prestart event 125.
  THcConfigEvtHandler *ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);
  // Add event handler for EPICS events
  THaEpicsEvtHandler *hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 181");
  gHaEvtHandlers->Add(hcepics);
  // Add event handler for scaler events
  THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 1");
  hscaler->AddEvtType(1);
  hscaler->AddEvtType(2);
  hscaler->AddEvtType(3);
  hscaler->AddEvtType(4);
  hscaler->AddEvtType(131);
  hscaler->SetDelayedType(131);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  // Add event handler for DAQ configuration event
  THcConfigEvtHandler *config = new THcConfigEvtHandler("config", "Hall C configuration event handler");
  gHaEvtHandlers->Add(config);

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
  // THcRun* run = new THcRun( pathList, Form(RunFileNamePattern, RunNumber) );
  // Could lead to an infinite loop, all segments in range analyzed.

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

  auto *run = new Podd::MultiFileRun(pathList, fileNames);
  // THcRun *run =
  //     new THcRun(pathList, Form(RunFileNamePattern, RunNumber)); // FIXME: Ultimately will want to use MiltiFileRun

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");

  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();

  // Moved file naming from here to top of script for transparency.

  analyzer->SetCountMode(2); // 0 = counter is # of physics triggers
                             // 1 = counter is # of all decode reads
                             // 2 = counter is event number
  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(181);
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");
  // analyzer->SetCrateMapFileName("MAPS/db_cratemap_lad.dat");// Temp set it to only LAD to avoid error

  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define DEF-file
  analyzer->SetOdefFile("DEF-files/LAD_COIN/PRODUCTION/coin_production_lad_only.def");
  // Define cuts file
  analyzer->SetCutFile("DEF-files/LAD_COIN/PRODUCTION/CUTS/coin_production_cuts_lad_only.def"); // optional
  // File to record accounting information for cuts
  // analyzer->SetSummaryFile(SummaryFileName.Data()); // optional

  LADFilteredStreamBuf ladbuf(std::cout);
  ladbuf.addFilterString("HitList(event");
  // Start the actual analysis.
  analyzer->Process(run);
  // Create report file from template

  // analyzer->PrintReport("TEMPLATES/SHMS/PRODUCTION/pstackana_production.template", REPORTFileName.Data()); //

  analyzer->PrintReport("TEMPLATES/LAD/PRODUCTION/lstackana_production.template", REPORTFileName.Data()); //

  // optional
}
