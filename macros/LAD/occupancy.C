#include <algorithm>

void occupancy( int module, double Lx, double Ly, int axis=0,  double dt_ns = 150.0, const char* spectrometer = "bb", const char* db_append = "gem" ){
  //TTree *T = (TTree*) gFile->Get("T");

  //gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(53);

  gStyle->SetOptStat(0);
  
  gPad->cd();

  //TH1D *hnstripskeep = new TH1D("hnstripskeep",

  double area_cm2 = Lx*Ly; 

  TString histname;
  if( axis == 0 ){
    histname.Form("hClusterMultiplicity_%s_%s_m%d_U", spectrometer, db_append, module);
  } else {
    histname.Form("hClusterMultiplicity_%s_%s_m%d_V", spectrometer, db_append, module);
  }

  TH1D *htemp;
  gFile->GetObject( histname.Data(), htemp );

  double mean = htemp->GetMean();
  double rms = htemp->GetRMS();


  double xmin = (mean - 10.0*rms > -0.5 ) ? mean - 10.0*rms : -0.5;
  double xmax = mean + 10.0*rms;
  
  htemp->GetXaxis()->SetRangeUser( xmin, xmax );
  
  htemp->Draw();

  double ymax = htemp->GetBinContent(htemp->GetMaximumBin());

  htemp->GetYaxis()->SetRangeUser(0, 1.5*ymax );

  //  double xmin = std::max( -0.5,mean-10.0*rms);
  //double xmax = mean+10.0*rms;

  TPaveText *p = new TPaveText( xmin, 1.01*ymax, xmax, 1.5*ymax, "br" );

  p->SetFillStyle(0);
  p->SetBorderSize(0);

  TString text;

  text.Form("Module: %i", module);

  p->AddText( text.Data() );

  double rate = mean / area_cm2 / dt_ns * 1.e6;
  
  text.Form( "Average clusters/event = %6.2f", mean );

  p->AddText( text.Data() );
  
  text.Form( "Hit rate = %6.2f (kHz/cm^{2})", rate );

  p->AddText( text.Data() );

  p->Draw();

}
