#include <algorithm>

//Plot raw strip multiplicity and occupancy:
void stripmult( int layer=0, int axis=0, int nstrips=3840) {
  //gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(53);

  gStyle->SetOptStat(0);
  gPad->cd();

  TH2D *htemp;
  if( axis == 0 ){
    gFile->GetObject( "h2_gem_NstripsU_layer", htemp );
  } else {
    gFile->GetObject( "h2_gem_NstripsV_layer", htemp );
  }

  TString hnametemp;
  TH1D *htemp_yproj = htemp->ProjectionY( hnametemp.Format( "htemp_yproj_layer%d_axis%d", layer, axis ), layer+1, layer+1 );

  TString htitle;
  if( axis == 0 ){
    htemp_yproj->SetTitle( htitle.Format( "Layer %d U strips; Number of strips fired;", layer ) );
  } else {
    htemp_yproj->SetTitle( htitle.Format( "Layer %d V strips; Number of strips fired;", layer ) );
  }
  
  double mean = htemp_yproj->GetMean();
  double rms = htemp_yproj->GetRMS();

  double xmin = (mean - 10.0*rms > -0.5 ) ? mean - 10.0*rms : -0.5;
  double xmax = mean + 10.0*rms;
  
  htemp_yproj->GetXaxis()->SetRangeUser( xmin, xmax );
  
  htemp_yproj->Draw();

  double ymax = htemp_yproj->GetBinContent(htemp_yproj->GetMaximumBin());

  htemp_yproj->GetYaxis()->SetRangeUser(0, 1.5*ymax );

  //double xmin = std::max( -0.5,mean-10.0*rms);
  //double xmax = mean+10.0*rms;

  TPaveText *p = new TPaveText( xmin, 1.01*ymax, xmax, 1.5*ymax, "br" );

  p->SetFillStyle(0);
  p->SetBorderSize(0);

  TString text;

  //double rate = mean / area_cm2 / dt_ns * 1.e6;
  //double nstripmean = mean;
  double occupancy = mean/double(nstrips);
  
  text.Form( "Raw strip multiplicity = %6.1f", mean );

  p->AddText( text.Data() );
  
  text.Form( "Raw occupancy = %6.2f%%", occupancy * 100.0);

  p->AddText( text.Data() );

  p->Draw();

}
