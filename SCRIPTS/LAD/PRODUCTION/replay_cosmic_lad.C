// #include "../LAD_link_defs.h" used for debugging. Not necessary.
void replay_cosmic_lad(Int_t RunNumber = 0, Int_t MaxEvent = 0) {

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

  // Create file name patterns.
  // const char *RunFileNamePattern = "lad_esb_%01d.evio.0";
  const char *RunFileNamePattern = "ladvme1_%03d.dat.0";
  vector<TString> pathList;
  pathList.push_back(".");
  // pathList.push_back("./ROOTfiles/COSMICS/raw/");
  pathList.push_back("/cache/hallc/c-lad/raw/");


  const char *ROOTFileNamePattern = "ROOTfiles/COSMICS/LAD_cosmic_hall_%d_%d.root";

  // Load Global parameters
  // Add variables to global list.
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/COSMIC/standard_cosmic.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);

  // Load fadc debug parameters
  // gHcParms->Load("PARAM/HMS/GEN/h_fadc_debug.param");

  // const char* CurrentFileNamePattern = "low_curr_bcm/bcmcurrent_%d.param";
  // gHcParms->Load(Form(CurrentFileNamePattern, RunNumber));

  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();
  gHcDetectorMap->Load("MAPS/LAD/DETEC/HODO/lhodo.map");
  // gHcDetectorMap->Load("MAPS/LAD/DETEC/HODO/lhodo_cosmic.map");

  // Add LAD detector
  THcLADSpectrometer *LAD = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(LAD);

  THcLADHodoscope *lhod = new THcLADHodoscope("hod", "LAD Hodoscope");
  LAD->AddDetector(lhod);

  // THcBCMCurrent* hbc = new THcBCMCurrent("H.bcm", "BCM current check");
  // gHaPhysics->Add(hbc);

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
  TString ROOTFileName = Form(ROOTFileNamePattern, RunNumber, MaxEvent);
  analyzer->SetCountMode(2); // 0 = counter is # of physics triggers
                             // 1 = counter is # of all decode reads
                             // 2 = counter is event number

  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(181);
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap_cosmic.dat");
  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define output DEF-file
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/lstackana_production_all.def");
  // Define cuts file
  // analyzer->SetCutFile("DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts_cosmic.def"); // optional
  // File to record cuts accounting information for cuts
  // analyzer->SetSummaryFile(
  //     Form("REPORT_OUTPUT/HMS/PRODUCTION/summary_production_%d_%d.report", RunNumber, MaxEvent)); // optional
  // Start the actual analysis.
  analyzer->Process(run);
  // Create report file from template.
  analyzer->PrintReport("TEMPLATES/LAD/PRODUCTION/lstackana_production.template",
                        Form("REPORT_OUTPUT/LAD/PRODUCTION/replay_lad_production_%d_%d.report", RunNumber, MaxEvent));
}
