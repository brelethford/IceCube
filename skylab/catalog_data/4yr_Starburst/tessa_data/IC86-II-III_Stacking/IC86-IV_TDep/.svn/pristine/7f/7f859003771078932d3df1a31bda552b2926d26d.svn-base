
void load_ark_ic79_Final_tdep_at(I3Ark& ark, bool opt_UseRealData) {
  cout << "\n-----  LOADING IC-79 BDT  -----\n";

  bool opt_UpScale = false;

  //
  // CONFIGURE EVENTLOADER WITH BACKGROUND / REAL DATA SAMPLE
  //

  ark.evLoader.SetBkgTree( LoadTree_IC79_GoodRuns_BDT_316d_small() );
  
  cout << "Loaded tree" << endl;
    

  if(opt_UpScale) {
    double datalivetime  = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");

    ark.livetime = 318.059 * 86400.; // This is what we want to simulate

    double upScaleFactor = ark.livetime / datalivetime;
    
    cout << "Upscaling factor: " << upScaleFactor << endl;
    
    ark.evLoader.SetBkgLoadMethod_PoissonSamplePlus(TStringify(upScaleFactor));
    ark.evLoader.SetZenithRangeDeg(0., 180., 2.0);
  }
  else
    {
      ark.livetime = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");
      ark.tmin = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmin");
      ark.tmax = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmax");
      ark.evLoader.SetBkgLoadMethod_Exact();

    }
  

  if (opt_UseRealData) {
    ark.evLoader.SetTimeMethod_Actual("timeMJD"); 
    cout << "Will load  * * * * U N B L I N D E D * * * *  Base Data!\n";
  }
  else {
    
    ark.evLoader.SetTimeMethod_Scramble();  // for blindness!
    cout << "Will load scrambled times for Base Data.\n";
  }
  
  cout << "Livetime (Days): " << ark.livetime/86400. << "\n";
  
  
  //
  // CONFIGURE EVENTLOADER WITH SIGNAL SIMULATION DATA SAMPLE
  //
  // evLoader will load signal events for any specified declination

  cout << "Configuring Source Event Sample...\n";

  ark.sourceZenWidthDeg = 0.5;

  ark.evLoader.SetSourceTree( LoadTree_IC79_nugen_numu_6308_9981f_small() );
  ark.evLoader.SetSourceZenWidthDeg(ark.sourceZenWidthDeg);


  //
  // SPECIFY CUTS
  //

  
  // Cuts can be added or reset any time, and events re-loaded with new cuts
  // (cuts are also evaluated in order, for possible efficiency gain)

  
  
  //TCut ic79_Cut = "BDT1.bdt_1 > 0.05 && BDT2.bdt_2 > 0.05 && BDT3.bdt_1 > 0.05 && BDT4.bdt_2 > 0.05";
    
  TCut level3 = "((L3CutBool == 1 && SplitDecision == 0)||(SplitDecision == 1 && (L3CutBool_TT0 == 1 || L3CutBool_TT1 == 1 || L3CutBool_TT2 == 1)))";



  TCut TrainUp = "(MPE_FINAL_MUON_Zd >= 125 &&  MPE_FINAL_MUON_rlogl < 12 && MPE_Pb_FINAL_MUON_SigmaDeg < 10 && MPE_FINAL_MUON_LDirC > 100 && MPE_Pb_PbfStatus == 0) ";

  TCut TrainHo = "(MPE_FINAL_MUON_Zd >= 85 && MPE_FINAL_MUON_Zd < 125 &&  MPE_FINAL_MUON_rlogl < 12 && MPE_Pb_FINAL_MUON_SigmaDeg < 10 && MPE_FINAL_MUON_LDirC > 100 && MPE_Pb_PbfStatus == 0) ";


  TCut BDTPS1 = (level3&&TrainUp&&"(BDT1.bdt_1 > 0.25 || BDT3.bdt_1 > 0.25)");
  TCut BDTPS2 = (level3&&TrainHo&&"(BDT2.bdt_2 > 0.25 || BDT4.bdt_2 > 0.25)");

  TCut BDTPS = (BDTPS1||BDTPS2)&&"(BDT1.bdt_1 > 0.05 && BDT2.bdt_2 > 0.05 && BDT3.bdt_1 > 0.05 && BDT4.bdt_2 > 0.05)";
  
  TCut ic79_Cut = "((((((L3CutBool == 1 && SplitDecision == 0)||(SplitDecision == 1 && (L3CutBool_TT0 == 1 || L3CutBool_TT1 == 1 || L3CutBool_TT2 == 1))))&&((MPE_FINAL_MUON_Zd >= 125 &&  MPE_FINAL_MUON_rlogl < 12 && MPE_Pb_FINAL_MUON_SigmaDeg < 10 && MPE_FINAL_MUON_LDirC > 100 && MPE_Pb_PbfStatus == 0) ))&&((BDT1.bdt_1 > 0.27 || BDT3.bdt_1 > 0.25)))||(((((L3CutBool == 1 && SplitDecision == 0)||(SplitDecision == 1 && (L3CutBool_TT0 == 1 || L3CutBool_TT1 == 1 || L3CutBool_TT2 == 1))))&&((MPE_FINAL_MUON_Zd >= 85 && MPE_FINAL_MUON_Zd < 125 &&  MPE_FINAL_MUON_rlogl < 12 && MPE_Pb_FINAL_MUON_SigmaDeg < 10 && MPE_FINAL_MUON_LDirC > 100 && MPE_Pb_PbfStatus == 0) ))&&((BDT2.bdt_2 > 0.23 || BDT4.bdt_2 > 0.23))))&&((BDT1.bdt_1 > 0.07 && BDT2.bdt_2 > 0.07 && BDT3.bdt_1 > 0.05 && BDT4.bdt_2 > 0.05 && (BDT2.bdt_2 > 0.14 || (mZd < 120 || mZd > 140))))";


  TCut ic79_Cut = "1";
  
  ark.evLoader.AddCut(ic79_Cut);
  
  
  
  //
  // SET NAMES OF VARIABLES TO BE USED FROM ROOT FILES
  //

  TString recoZenRadName = "mZr"; // needed for eProb below
  ark.evLoader.SetName_recoZenith_rad(recoZenRadName);
  ark.evLoader.SetName_recoAzimuth_rad("mAr");
  ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC79(mpbSigmaDeg,mueEn)");
  //ark.evLoader.SetName_sigmaDeg("mpfSigmaDegRescaledIC79V2(mpbSigmaDeg,mueEn)");
  //ark.evLoader.SetName_sigmaDeg("kkk(mpbSigmaDeg,mueEn)");
  ark.evLoader.SetName_runID("Run");
  ark.evLoader.SetName_eventID("EventID");
  // EXAMPLE: How to load different vars for source and bkg, if necessary
  //  dsBkg->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg")
  //  dsSrc->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg*1.4");
  //  evLoader.SetName_sigmaDeg("myRecoSigmaDeg");
  

  // The choice here has consequences for how to fill eProb, below
  TString energyVar = "log10(mueEn)";
  //TString energyVar = "log10(TruncatedMean_FINAL_MUON_Bins_En)";

  ark.evLoader.SetName_energyValue(energyVar);


  //
  // LOAD EVENTS, SET UP ENERGY PDFs
  // 

  cout << "Loading Background Events: " << ark.evLoader.GetBkgTree()->GetTitle() << endl;
  cout << "Using cut: " << ark.evLoader.GetCuts().GetTitle() << endl;
  ark.evLoader.SetMonitor(true);
  ark.evLoader.LoadBkgEvents(ark.baseEvents);

  

  ZenithEnergyProb* zen_eProb = new ZenithEnergyProb();
  
 
  cout << "Filling Energy PDFs:\n";
  zen_eProb->SetSourceZenWidthDeg( ark.sourceZenWidthDeg);
  zen_eProb->SetName_recoZenith_rad(recoZenRadName);

  vector<double> zenMinDegVect;
      
  // Fill vector with bin edges that match "CutDMS" bin edges
  // 0 - 90  zen added
  double tempBot=1;

  //Old binning up to 90 deg                                                                                                                                                                                                                                                                                                                                                                                                            
  /*                                                                                                                                                                                                                                                                                                                                                                                                                                    
    for(int i=0; i<21;i++){                                                                                                                                                                                                                                                                                                                                                                                                             
    //cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;                                                                                                                                                                                                                                                                                                                                                                
    zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());                                                                                                                                                                                                                                                                                                                                                                           
    tempBot-=0.05;                                                                                                                                                                                                                                                                                                                                                                                                                      
  }                                                                                                                                                                                                                                                                                                                                                                                                                                     
  */


  //New binning up 84 deg and then a much finner bin                                                                                                                                                                                                                                                                                                                                                                                    

  for(int i=0; i < 19;i++){
    //cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;                                                                                                                                                                                                                                                                                                                                                                
    zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
    tempBot-=0.05;
  }

 
  //zenMinDegVect.push_back(85.00); 
  zenMinDegVect.push_back(86.00); 
  //zenMinDegVect.push_back(87.00); 
  zenMinDegVect.push_back(88.00); 
  //zenMinDegVect.push_back(89.00); 
  zenMinDegVect.push_back(90.00);
 

  // A reasonable split for the upgoing range
  //zenMinDegVect.push_back(110.00); // ~acos(-0.33)
  //zenMinDegVect.push_back(132.00); // ~acos(-0.67)

  zenMinDegVect.push_back(97.); // Follow BDT regions
  zenMinDegVect.push_back(125.);
  zenMinDegVect.push_back(180.);

  zen_eProb->SetZenithBandsDeg(zenMinDegVect);

  zen_eProb->SetLoadModeNew(true); // true is now faster

  if (1) { } // new default is to constrain signal inside histogram
  if (0) { zen_eProb->SetConstrainSignal(false); }  // the old way

  if (energyVar == "log10(mueEn)") {
  //  if (energyVar == "log10(TruncatedMean_FINAL_MUON_Bins_En)") {
    int nBackFill = 35; // don't backfill previous bins
    zen_eProb->SetEnergyGammaRangeAndBackFill(40,2.,9., 30,1.,4., nBackFill);
  }

  zen_eProb->SetTableBkg(ark.baseEvents);
  TStopwatch ts;
  zen_eProb->SetTableGamma(ark.evLoader.GetSourceTree(), ark.evLoader.GetCuts(), energyVar);
  ts.Print();

  ark.eProb = zen_eProb;

  
  // This seems to be okay for ~ 5000 events or more:
  // 1. deg smoothing and 180 declination maps (e.g. 1 deg binning)
  
  ark.decBkgProb.Initialize(180, 1);
  ark.decBkgProb.SetBaseDecMap(ark.baseEvents);
  
  ark.lcBkgProb.Initialize(12.,90.,true); //so why doesn't this work together?
  ark.lcBkgProb.FillLCBkgHisto(ark.baseEvents);
   
  psData = new I3Analysis();

  psData->SetBkgSpaceProb(ark.decBkgProb);
  psData->SetBaseEvents(ark.baseEvents);
  psData->SetEnergyProb(*(ark.eProb));

  

  EventTimeModule * ic79times = new EventTimeModuleDiscrete();
  

  if (opt_UseRealData) {   
    psData->UseRealData();  // the event set is now exactly  to the data set (i.e. no scrambling, no fake signal added.)
    ic79times->SetTimesWithinRange(ark.tmin, ark.tmax);
    psData->SetEventTimeModulePtr(ic79times);
    cout << "UNBLINDING: Using real times" << endl;
  }
  
  else { 
    psData->SetRandomizeBase(true);
    ic79times->SetTimesFromMJDFile("/atlas/users/christov/psLab/macro_llh/IC86-IV_TDep/HugeListOfIC79TimesNewMethod.txt");
    psData->SetEventTimeModulePtr(ic79times);
    psData->GenerateDataSet_with_nSrcEvents(0);
 }
 
 ark.psData = psData;

 cout << "----- IC-79 BDT Loaded -----\n";
 
}
