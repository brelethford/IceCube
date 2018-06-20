
void load_ark_ic86_II_III_IV(I3Ark& ark, bool opt_UseRealData, TString RecoName) {

  cout << "\n-----  LOADING IC-86 2, 3 and 4 -----\n";

  bool opt_UpScale = false;

  //
  // CONFIGURE EVENTLOADER WITH BACKGROUND / REAL DATA SAMPLE
  //

  ark.evLoader.SetBkgTree( LoadTree_IC86_II_III_IV_GoodRuns_Full() );
  
  cout << "Loaded tree" << endl;
    

  if(opt_UpScale) {
    double datalivetime  = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");
    
    double ic86_II_days = 330.38, ic86_III_days = 359.95, ic86_IV_days = 367.21;	
    ark.livetime = (ic86_II_days + ic86_III_days + ic86_IV_days)*86400.0; //340  * 86400.; // This is what we want to simulate

    double upScaleFactor = ark.livetime / datalivetime;
    
    cout << "Upscaling factor: " << upScaleFactor << endl;
    
    ark.evLoader.SetBkgLoadMethod_PoissonSamplePlus(TStringify(upScaleFactor));
    // jake - i changed the min zenith from 0 to 85...correct?
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

  ark.evLoader.SetSourceTree( LoadTree_IC86_II_III_IV_nugen_numu_E_1() );
  ark.evLoader.SetSourceZenWidthDeg(ark.sourceZenWidthDeg);

//ark.evLoader.SetSourceZenWidthDegRange(0.5, 80., 90.);

  //
  // SPECIFY CUTS
  //

  TCut IC86_Cut = "1";

  // Cuts can be added or reset any time, and events re-loaded with new cuts
  // (cuts are also evaluated in order, for possible efficiency gain)

  ark.evLoader.AddCut(IC86_Cut);


  //
  // SET NAMES OF VARIABLES TO BE USED FROM ROOT FILES
  //

  TString recoZenRadName = RecoName+"Zr"; // needed for eProb below
  ark.evLoader.SetName_recoZenith_rad(recoZenRadName);
  ark.evLoader.SetName_recoAzimuth_rad(RecoName+"Ar");
  if (RecoName.CompareTo("MPEFit_TT")==0) {
    cout << "Rescaling using MPE...\n";
    ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC86_MPE(mpbSigmaDeg,muexEn)");
  }
  else if (RecoName.CompareTo("MuEXAngular4")==0) {
    cout << "Rescaling using MuEXAngular4...\n";
    ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC86_MuEX(mMuEXSigmaDeg,muexEn)");
  }
  else if (RecoName.CompareTo("SplineMPE") == 0) {
    cout << "Rescaling using SplineMPE...\n";
    ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC86_II_III_IV_SplineMPE(SplinePbSigmaDeg,muexEn, SplineMPEZd)");
  }
  else cout <<"AHHHHHHHHHHHHHHHHHHHH";

  ark.evLoader.SetName_runID("RunID");
  ark.evLoader.SetName_eventID("EventID");
  // EXAMPLE: How to load different vars for source and bkg, if necessary
  //  dsBkg->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg")
  //  dsSrc->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg*1.4");
  //  evLoader.SetName_sigmaDeg("myRecoSigmaDeg");


  // The choice here has consequences for how to fill eProb, below
  TString energyVar = "log10(muexEn)";

  ark.evLoader.SetName_energyValue(energyVar);


  //
  // LOAD EVENTS, SET UP ENERGY PDFs
  // 

  cout << "Loading Background Events: " << ark.evLoader.GetBkgTree()->GetTitle() << endl;

  ark.evLoader.SetMonitor(true);
  ark.evLoader.LoadBkgEvents(ark.baseEvents);

  ZenithEnergyProb* zen_eProb = new ZenithEnergyProb();

 
  //Adding mc cut
//   ark.evLoader.AddCut(RecoName+"DelAng < 7. && log10(muexEn) >= 0");
//   cout << "Using cut: " << ark.evLoader.GetCuts().GetTitle() << endl;
  
  cout << "Filling Energy PDFs:\n";
  zen_eProb->SetSourceZenWidthDeg( ark.sourceZenWidthDeg);
  zen_eProb->SetName_recoZenith_rad(recoZenRadName);

  vector<double> zenMinDegVect;
      
  // Fill vector with bin edges that match "CutDMS" bin edges
  // 0 - 90  zen added
  zenMinDegVect.push_back(acos(1.)*TMath::RadToDeg());
  zenMinDegVect.push_back(acos(0.97)*TMath::RadToDeg());
  zenMinDegVect.push_back(acos(0.94)*TMath::RadToDeg());
  double tempBot=0.9;
  //for(int i=0; i<21;i++){
    for(int i=0; i<17;i++){
    cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
    zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
    tempBot-=0.05;
  }
  
  //zenMinDegVect.push_back(85.00); 
  zenMinDegVect.push_back(86.00); 
  //zenMinDegVect.push_back(87.00); 
  zenMinDegVect.push_back(88.00); 
  //zenMinDegVect.push_back(89.00); 
  //zenMinDegVect.push_back(90.00);

  // A reasonable split for the upgoing range
  //zenMinDegVect.push_back(105.00); //? Follow Chad's example
  //jake - try more zen bins, spacing of 0.2 in cos(zen)
  //zenMinDegVect.push_back(110.00); // ~acos(-0.33)
  //zenMinDegVect.push_back(140.00); // ~acos(-0.67)
  //zenMinDegVect.push_back(160.00); // ~acos(-0.67)
  tempBot-=0.05;
  //for(int i=0; i<9;i++){
  for(int i=0; i<5;i++){
    cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
    zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
    //tempBot-=0.1;
    tempBot-=0.2;
  }
  zenMinDegVect.push_back(180.00);

  zen_eProb->SetZenithBandsDeg(zenMinDegVect);

  zen_eProb->SetLoadModeNew(true); // true is now faster

  if (1) { } // new default is to constrain signal inside histogram
  if (0) { zen_eProb->SetConstrainSignal(false); }  // the old way

  if (energyVar == "log10(muexEn)") {
    int nBackFill = 35; // don't backfill previous bins
    cout << "backfill" << endl;
    zen_eProb->SetEnergyGammaRangeAndBackFill(40,2.,9., 30,1.,4., nBackFill);
  }

  zen_eProb->SetTableBkg(ark.baseEvents);
  TStopwatch ts;
  cout << "set table gamma" << endl;
  zen_eProb->SetTableGamma(ark.evLoader.GetSourceTree(), ark.evLoader.GetCuts(), energyVar);
  ts.Print();

  ark.eProb = zen_eProb;
  cout<<" TEST 1 "<<endl;

  // This seems to be okay for ~ 5000 events or more:
  // 1. deg smoothing and 180 declination maps (e.g. 1 deg binning)
    
  ark.decBkgProb.Initialize(180, 1);
  ark.decBkgProb.SetBaseDecMap(ark.baseEvents);
  
  //ark.lcBkgProb.Initialize(12.,90.,true); //so why doesn't this work together?
  //ark.lcBkgProb.FillLCBkgHisto(ark.baseEvents);
  cout<<" TEST 2 "<<endl; 
  psData = new I3Analysis();

  psData->SetBkgSpaceProb(ark.decBkgProb);
  psData->SetBaseEvents(ark.baseEvents);
  psData->SetEnergyProb(*(ark.eProb));
  cout<<" TEST 3 "<<endl;
  /*
  EventTimeModule * ic79times = new EventTimeModuleDiscrete();

  if (opt_UseRealData) {
    psData->UseRealData();  // the event set is now exactly equal
    // to the data set (i.e. no scrambling, no fake signal added.)
  } else {
    psData->SetRandomizeBase(true);
    ic79times->SetTimesFromMJDFile("/net/user/mfbaker/proto_multi/macro_llh/ic79/HugeListOfIC79TimesNewMethod.txt");
    psData->SetEventTimeModulePtr(ic79times);
    //psData->GenerateDataSet_with_nSrcEvents(0); // needs an I3SignalGenerator first??
    } 
  */



    
  
   
if (opt_UseRealData) {   
  psData->UseRealData();  // the event set is now exactly  to the data set (i.e. no scrambling, no fake signal added.)
  cout<<" TEST 4 "<<endl;
 }
 
 else { 
   //psData->SetRandomizeBase(true);
   //ic79times->SetTimesFromMJDFile("/net/user/mfbaker/proto_multi/macro_llh/ic79/HugeListOfIC79TimesNewMethod.txt");
   // psData->SetEventTimeModulePtr(ic79times);
   psData->GenerateDataSet_with_nSrcEvents(0);
   cout<<" TEST 5 "<<endl;
 }
 
 ark.psData = psData;
 cout<<" TEST FINAL LOADARK"<<endl;
 cout << "----- IC-86 2 and 3 Loaded -----\n";
 
}
