#include "MultiFileRun.h"

void replay_production_LAD_HMS(int RunNumber=0, int MaxEvent=0, int FirstEvent = 1, int MaxSegment = 1, int FirstSegment = 0, const char* fname_prefix = "nps_coin")
{

  // Get RunNumber and MaxEvent if not provided.
  if(RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if( RunNumber<=0 ) return;
  }
  if(MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if(MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      exit;
    }
  }

  // Create file name patterns.
  //  const char* RunFileNamePattern="NPS_3crate_%d.evio.0"; Need to change to LAD
  const char* RunFileNamePattern="%s_%d.dat.%u"; 
  vector<string> pathList;
  pathList.push_back(".");
  pathList.push_back("./raw");
  pathList.push_back("./raw/../raw.copiedtotape");
  pathList.push_back("./cache");
  pathList.push_back("/net/cdaq/cdaql1data/coda/data/raw");

  const char* ROOTFileNamePattern;
  if (MaxEvent == 50000 && FirstEvent == 1){
    ROOTFileNamePattern = "ROOTfiles/COIN/50k/lad_hms_%d_%d_%d.root"; // Current run_gen will have the first event = 0 so it doesnt interfere here.
  }
  else if (MaxEvent == -1 && (FirstSegment - MaxSegment) == 0) {
    ROOTFileNamePattern = "ROOTfiles/COIN/PRODUCTION/lad_hms_coin_%d_%d_%d_%d.root";
  }  else{
    ROOTFileNamePattern = "ROOTfiles/COIN/PRODUCTION/lad_hms_coin_%d_%d_%d.root";
  }
  
  
  // Add variables to global list.
  gHcParms->Define("gen_run_number", "Run Number", RunNumber); 
  gHcParms->AddString("g_ctp_database_filename", "DBASE/NPS/standard_coin.database");  // FIXME: DBASE FOR HMS+NPS;
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_det_calib_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_bcm_calib_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_bpm_calib_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_optics_filename"));
  // Load parameters for SHMS trigger configuration
  gHcParms->Load(gHcParms->GetString("g_ctp_trig_config_filename"));
  // Load hpcentral momentum offset 
  //gHcParms->Load("PARAM/HMS/GEN/hpcentral_function_sp18.param");
  // Load fadc debug parameters
  gHcParms->Load("PARAM/HMS/GEN/h_fadc_debug_sp18.param");


   
  // Load the Hall C style detector map 
  gHcDetectorMap = new THcDetectorMap();
  gHcDetectorMap->Load("MAPS/NPS/DETEC/pcal_nps_coin.map"); //FIXME: CHANGE TO COIN MAP
  //gHcDetectorMap->Load("MAPS/NPS/DETEC/pcal_nps_standard.map");

  /*// Dec data
  gHaApps->Add(new Podd::DecData("D","Decoder raw data")); //FIXME: NEED THIS ONE?*/

    // Load BCM values


  //=:=:=
  // HMS 
  //=:=:=

  THcHallCSpectrometer* HMS = new THcHallCSpectrometer("H", "HMS");
  /*HMS->SetEvtType(2);
  HMS->AddEvtType(4);
  HMS->AddEvtType(5);
  HMS->AddEvtType(6);
  HMS->AddEvtType(7);*/
  gHaApps->Add(HMS);
  // Add drift chambers to HMS apparatus
  THcDC* hdc = new THcDC("dc", "Drift Chambers");
  HMS->AddDetector(hdc);
  // Add hodoscope to HMS apparatus
  THcHodoscope* hhod = new THcHodoscope("hod", "Hodoscope");
  HMS->AddDetector(hhod);
  // Add Cherenkov to HMS apparatus
  THcCherenkov* hcer = new THcCherenkov("cer", "Heavy Gas Cherenkov");
  HMS->AddDetector(hcer);
  // Add Aerogel Cherenkov to HMS apparatus
  //THcAerogel* haero = new THcAerogel("aero", "Aerogel");
  //HMS->AddDetector(haero);
  // Add calorimeter to HMS apparatus
  THcShower* hcal = new THcShower("cal", "Calorimeter");
  HMS->AddDetector(hcal);

  // Add rastered beam apparatus
  THaApparatus* hbeam = new THcRasteredBeam("H.rb", "Rastered Beamline");
  gHaApps->Add(hbeam);


  // Add physics modules
  // Calculate reaction point
  THcReactionPoint* hrp = new THcReactionPoint("H.react", "HMS reaction point", "H", "H.rb");
  gHaPhysics->Add(hrp);
  // Calculate extended target corrections
  THcExtTarCor* hext = new THcExtTarCor("H.extcor", "HMS extended target corrections", "H", "H.react");
  gHaPhysics->Add(hext);
  // Calculate golden track quantities
  THaGoldenTrack* gtr = new THaGoldenTrack("H.gtr", "HMS Golden Track", "H");
  gHaPhysics->Add(gtr);
  // Calculate primary (scattered beam - usually electrons) kinematics
  THcPrimaryKine* hkin = new THcPrimaryKine("H.kin", "HMS Single Arm Kinematics", "H", "H.rb");
  gHaPhysics->Add(hkin);
  // Calculate the hodoscope efficiencies
  THcHodoEff* heff = new THcHodoEff("hhodeff", "HMS hodo efficiency", "H.hod");
  gHaPhysics->Add(heff);

  // Add handler for scaler events
  THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 2");  
  hscaler->AddEvtType(1);
  //  hscaler->AddEvtType(140);
  hscaler->AddEvtType(129);
  hscaler->AddEvtType(130);
  hscaler->SetDelayedType(129);
  hscaler->SetDelayedType(130);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  // Add event handler for helicity scalers
  THcHelicityScaler *hhelscaler = new THcHelicityScaler("H", "Hall C helicity scaler");
  //hhelscaler->SetDebugFile("HHelScaler.txt");
  hhelscaler->SetROC(5);
  hhelscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hhelscaler);

  // Add event handler for DAQ configuration event
  THcConfigEvtHandler *hconfig = new THcConfigEvtHandler("hconfig", "Hall C configuration event handler");
  gHaEvtHandlers->Add(hconfig);
  
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Kinematics Modules
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // Add Physics Module to calculate primary (scattered electrons) beam kinematics
  //THcPrimaryKine* hkin_primary = new THcPrimaryKine("H.kin.primary", "HMS Single Arm Kinematics", "H", "H.rb");
  //gHaPhysics->Add(hkin_primary);
  // Add Physics Module to calculate secondary (scattered hadrons) beam kinematics
  THcNPSSecondaryKine* pkin_secondary = new THcNPSSecondaryKine("NPS.kin.secondary", "NPS Secondary Kinematics", "P", "H.kin",0.0,"H.react");
  gHaPhysics->Add(pkin_secondary); 


  //=:=:=
  // NPS
  //=:=:=

  //Add NPS spectrometer apparatus
  THaApparatus* NPS = new THcNPSApparatus("NPS","NPS");
  gHaApps->Add(NPS);

  //Add NPS Calorimeter to NPS apparatus
  THcNPSCalorimeter* cal = new THcNPSCalorimeter("cal", "Calorimeter");
  NPS->AddDetector(cal);

  THcNPSTrackInfo* nps_trk = new THcNPSTrackInfo("NPS.cal.trk", "NPS Track Info", "NPS", "H.react");
  gHaPhysics->Add(nps_trk);

  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Global Objects & Event Handlers
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // Add trigger apparatus
  THaApparatus* TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);

  // Add trigger detector to trigger apparatus
  /*THcTrigDet* coin = new THcTrigDet("coin", "Coincidence Trigger Information");
  // Suppress missing reference time warnings for these event types
  coin->SetEvtType(1);
  coin->AddEvtType(2);
  TRG->AddDetector(coin);*/

  THcTrigDet* hms = new THcTrigDet("hms", "HMS Trigger Information");
  TRG->AddDetector(hms); 

  THcHelicity* helicity = new THcHelicity("helicity", "Helicity Detector");
  TRG->AddDetector(helicity); // check later

  //Add coin physics module THcNPSCoinTime::THcNPSCoinTime (const char *name, const char* description, const char* elecArmName, 
  // const char* hadArmName, const char* coinname) :
  THcNPSCoinTime* coinTime = new THcNPSCoinTime("CTime", "Coincidence Time Determination", "NPS", "H", "T.hms");
  gHaPhysics->Add(coinTime);

  // Add event handler for EPICS events
  THaEpicsEvtHandler* hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 180");
  gHaEvtHandlers->Add(hcepics);
  
  // Add handler for prestart event 125.
  THcConfigEvtHandler* ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);

  THcNPSConfigEvtHandler* npsconfig = new THcNPSConfigEvtHandler("npsconfig", "NPS configuration event handler");
  // Four parameters are added by default (see THcNPSConfigEvtHandler::Init) 
  // To add more parameters, one can use: npsconfig->AddParameter("parameter name", "key name")
  // e.g: npsconfig->AddParameter("fadc250_nsa", "FADC250_NSA");
  gHaEvtHandlers->Add(npsconfig);

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Acpparatus's and PhysicsModules,
  // and executes the output routines.
  
  THcAnalyzer* analyzer = new THcAnalyzer;

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  //Define the run(s) that we want to analyze.
  //THcRun* run = new THcRun( pathList, Form(RunFileNamePattern, RunNumber) );
  vector<string> fileNames = {};
  for(Int_t iseg = FirstSegment; iseg <= MaxSegment; iseg++) {
    TString codafilename;
    codafilename.Form(RunFileNamePattern, fname_prefix, RunNumber, iseg);
    cout << "codafilename = " << codafilename << endl;
    fileNames.emplace_back(codafilename.Data());
  }
  auto* run = new Podd::MultiFileRun( pathList, fileNames);

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");
  run->SetEventRange(FirstEvent, MaxEvent);    
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();
  
  // Define the analysis parameters
  TString ROOTFileName;
  if(MaxEvent == -1 && (FirstSegment - MaxSegment) == 0) {
    ROOTFileName = Form(ROOTFileNamePattern, RunNumber, FirstSegment, FirstEvent, MaxEvent);
  } else {
     ROOTFileName = Form(ROOTFileNamePattern, RunNumber, FirstEvent, MaxEvent);
  }

  // Define the analysis parameters
  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(180);
  analyzer->AddEpicsEvtType(181);
  analyzer->AddEpicsEvtType(182);
  analyzer->AddEpicsEvtType(183);
  analyzer->AddEpicsEvtType(184);
  analyzer->AddEpicsEvtType(185);
  analyzer->AddEpicsEvtType(186);
  analyzer->AddEpicsEvtType(187);
  analyzer->AddEpicsEvtType(188);
  analyzer->AddEpicsEvtType(189);
  analyzer->AddEpicsEvtType(190);
  analyzer->AddEpicsEvtType(191);
  
  analyzer->SetCountMode(2);  // 0 = counter is # of physics triggers
                              // 1 = counter is # of all decode reads
                              // 2 = counter is event number

  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/NPS/CRATE/db_cratemap_coin.dat") ; //FIXME: CHANGE
  // Define DEF-file+
  analyzer->SetOdefFile("DEF-files/NPS/NPS_coin.def"); //FIXME: CHANGE
  // Define cuts file
  analyzer->SetCutFile("DEF-files/NPS/NPS_cuts_coin.def"); //FIXME: CHANGE
  // File to record accounting information for cuts
  //analyzer->SetSummaryFile(Form("REPORT_OUTPUT/COIN/PRODUCTION/summary_production_%d_%d.report", RunNumber, MaxEvent));  // optional //FIXME: CHANGE
  // start the actual analysis
  analyzer->Process(run);     
  // Create report file from template.
  if(MaxEvent == -1 && (FirstSegment - MaxSegment) == 0) {
    analyzer->PrintReport("TEMPLATES/NPS/NPS_coin.template",
			  Form("REPORT_OUTPUT/COIN/coin_NPS_HMS_report_%d_%d_%d_%d.report",
			  RunNumber, FirstSegment, FirstEvent, MaxEvent)); //FIXME:CHANGE
  } else {
    analyzer->PrintReport("TEMPLATES/NPS/NPS_coin.template",
			  Form("REPORT_OUTPUT/COIN/coin_NPS_HMS_report_%d_%d.report",
			  RunNumber, MaxEvent)); //FIXME:CHANGE
  }
  
}
