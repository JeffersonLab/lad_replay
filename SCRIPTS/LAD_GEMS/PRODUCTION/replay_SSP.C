void replay_SSP(int RunNumber=111, int MaxEvent=100)
{

  // datafile name pattern
  const char* RunFileNamePattern = "ssp_gem_%d.evio.0";
  vector<TString> pathList;
  pathList.push_back("./");
  pathList.push_back("./ROOTfiles/GEM_EVIO_COSMIC/");
  
  // Params
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  gHcParms->Load("PARAM/LAD/GEM/lgem_geom.param");
  // gHcParms->Load("DB_LAD/lgem_chan.map");
  // gHcParms->Load("DB_LAD/lgem_cuts.param");

  // Detector map
  //  gHcDetectorMap = new THcDetectorMap();
  //  gHcDetectorMap->Load("detector.map");
  
  THcLADSpectrometer* lad = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(lad);

  THcLADGEM* gem = new THcLADGEM("gem", "gem");
  lad->AddDetector(gem);  

  THcAnalyzer* analyzer = new THcAnalyzer;
  THaEvent* event = new THaEvent;  

  //  THcRun* run = new THcRun(pathList, Form(RunFileNamePattern, RunNumber) );
  THaRun* run = new THaRun(pathList, Form(RunFileNamePattern, RunNumber) );

  run->SetEventRange(1, MaxEvent);
  run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  run->Print();

  analyzer->SetEvent(event);
  analyzer->SetCountMode(2); // 2 = counter is event number

  analyzer->SetCrateMapFileName("MAPS/db_cratemap_gem_test.dat");
  analyzer->SetOutFile("ROOTfiles/test_ssp.root");
  analyzer->SetOdefFile( "DEF-files/LAD/PRODUCTION/GEM/test.def" );
  //  analyzer->SetCutFile( "test_cuts.def" );
  //  analyzer->SetCrateMapFileName("cratemap.dat");

  analyzer->Process(run);


}
