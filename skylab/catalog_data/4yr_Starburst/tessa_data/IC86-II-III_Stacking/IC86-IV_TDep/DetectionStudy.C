
TCanvas* DetectionStudy(Ark& ark, AnalysisFn& llh, DiscoveryPotential& disco)
{
  double meanSrcEv_ForFlux = ark.psData->GetMeanSrcNev();
  cout << " Mean Number of Source events expected for source model: ";
  cout << meanSrcEv_ForFlux << endl;

  llh.SetAnalysisSet(ark.psData);
  llh.SetSearchCoord(ark.mySrcLocation);

  disco.SetAnalysisSet(ark.psData);
  disco.SetAnalysisFn(&llh);

  disco.AnalyzeDiscoveryPotential();
  
  //
  // PLOT OUTPUT
  //

  // only make new canvas if it didn't already exist:
  canDetection = dynamic_cast<TCanvas*> 
    ( gROOT->GetListOfCanvases()->FindObject("canDetection") );
  if ( !canDetection ) {
    canDetection = new TCanvas("canDetection","canDetection",20,20,1000,500);
    canDetection->Divide(2,1);
  }
  canDetection->cd(1);
  disco.h_nSrc_logP->Draw("box");
  int min_nSrcBins = disco.nsMaxSampled_;

  if (min_nSrcBins <= 20) {
    disco.h_nSrc_logP->GetYaxis()->SetRange(0,min_nSrcBins+1);
    disco.h_nSrc_logP->GetYaxis()->CenterLabels(1);
    disco.h_nSrc_logP->GetYaxis()->SetNdivisions(min_nSrcBins+1);
  }
  
  canDetection->Update();
  vline(-log10(disco.GetDetectionSignificance() ),kRed)->Draw();

  canDetection->cd(2);
  disco.hDiscoveryFraction->GetYaxis()->SetRangeUser(0,1);
  disco.hDiscoveryFraction->Draw();
  if (min_nSrcBins <= 20) {
    disco.hDiscoveryFraction->GetXaxis()->SetRange(0,min_nSrcBins+2);
    //    hDiscoveryFraction->GetXaxis()->CenterLabels(1);
    disco.hDiscoveryFraction->GetXaxis()->SetNdivisions(min_nSrcBins+2);
  }
  canDetection->Update();
  

  TGraph *g = disco.detectionRate_vs_Mean;
  g->SetMarkerStyle(6);
  g->Draw("P");
  cout << "Poisson Mean nSrcEv for " << disco.GetDetectionPower()*100;
  cout << "+/-" << disco.powerUncertainty_*100;
  cout << " % detection probability: ";
  cout  << disco.MeanSrcEv_ForDetection_ << endl;

  canDetection->Update();
  vline(disco.MeanSrcEv_ForDetection_,4,2)->Draw();
  hline(disco.GetDetectionPower(),4,2)->Draw();

  cout << "Threshold significance for detection: ";
  cout  << disco.GetDetectionSignificance() << endl;

  double meanSrcEv_ForDetection = disco.MeanSrcEv_ForDetection_;
  double fluxScale = 
    ark.psData->GetFluxScaleForNev(meanSrcEv_ForDetection);

  cout << "  Flux scale relative to model is: " << fluxScale << endl;


  // can use this if you know you are dealing with a power-law, and
  // know the fluxConstant
  /*
  double meanFlux_ForDetection = fluxScale*fluxConstant;
  cout << "Mean Flux for " << disco.detectionPower_*100;
  cout << "% detection probability: ";
  cout << "  Flux: " << meanFlux_ForDetection << " GeV^-1 cm^-2 s^-1\n";
  */

  return canDetection;
}

double DetectionStudy_d(I3Ark& ark, AnalysisFn& llh, DiscoveryPotential& disco, double &fluxScale = 0)
{
  double meanSrcEv_ForFlux = ark.psData->GetMeanSrcNev();
  cout << " Mean Number of Source events expected for source model: ";
  cout << meanSrcEv_ForFlux << endl;

  llh.SetAnalysisSet(ark.psData);
  llh.SetSearchCoord(ark.mySrcLocation);

  disco.SetAnalysisSet(ark.psData);
  disco.SetAnalysisFn(&llh);

  disco.AnalyzeDiscoveryPotential();
  
  int min_nSrcBins = disco.nsMaxSampled_;

  cout << "Threshold significance for detection: ";
  cout  << disco.GetDetectionSignificance() << endl;
  double meanSrcEv_ForDetection = disco.MeanSrcEv_ForDetection_;
  fluxScale = ark.psData->GetFluxScaleForNev(meanSrcEv_ForDetection);

  cout << "  Flux scale relative to model is: " << fluxScale << endl;


  // can use this if you know you are dealing with a power-law, and
  // know the fluxConstant
  /*
  double meanFlux_ForDetection = fluxScale*fluxConstant;
  cout << "Mean Flux for " << disco.detectionPower_*100;
  cout << "% detection probability: ";
  cout << "  Flux: " << meanFlux_ForDetection << " GeV^-1 cm^-2 s^-1\n";
  */

  //return meanSrcEv_ForDetection;

  return meanSrcEv_ForDetection;
}

double DetectionStudyGauss_d(MultiArk& ark, MultiGaussAnalysisFn& llh, DiscoveryPotential& disco, double &fluxScale = 0){
     double meanSrcEv_ForFlux = ark.psData->GetMeanSrcNev();
  cout << " Mean Number of Source events expected for source model: ";
  cout << meanSrcEv_ForFlux << endl;

  llh.SetAnalysisSet(ark.psData);
  llh.SetSearchCoord(ark.mySrcLocation);

  disco.SetAnalysisSet(ark.psData);
  disco.SetAnalysisFn(&llh);

  disco.AnalyzeDiscoveryPotential();
  
  int min_nSrcBins = disco.nsMaxSampled_;

  cout << "Threshold significance for detection: ";
  cout  << disco.GetDetectionSignificance() << endl;
  double meanSrcEv_ForDetection = disco.MeanSrcEv_ForDetection_;
  fluxScale = ark.psData->GetFluxScaleForNev(meanSrcEv_ForDetection);

  cout << "  Flux scale relative to model is: " << fluxScale << endl;


  // can use this if you know you are dealing with a power-law, and
  // know the fluxConstant
  /*
  double meanFlux_ForDetection = fluxScale*fluxConstant;
  cout << "Mean Flux for " << disco.detectionPower_*100;
  cout << "% detection probability: ";
  cout << "  Flux: " << meanFlux_ForDetection << " GeV^-1 cm^-2 s^-1\n";
  */

  //return meanSrcEv_ForDetection;

  return meanSrcEv_ForDetection;
 
}
