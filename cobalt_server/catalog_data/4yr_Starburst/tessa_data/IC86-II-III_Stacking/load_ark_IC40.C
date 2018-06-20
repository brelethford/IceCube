
void load_ark_ic40(I3Ark& ark, bool opt_UseRealData) {
  // For now, this is where the meta-TTree loading scripts are:

  cout << "\n-----  LOADING IC-40  -----\n";


  //
  // CONFIGURE EVENTLOADER WITH BACKGROUND / REAL DATA SAMPLE
  //

  ark.evLoader.SetBkgTree( LoadTree_IC40_CutA6_data_Fix_GoodRuns_small() );
  ark.livetime = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");
  //ark.tmin = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmin");
  //ark.tmax = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmax");

  ark.tmin = 54562.372858;
  ark.tmax = 54971.150254;
  

  ark.evLoader.SetBkgLoadMethod_Exact();

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

  ark.evLoader.SetSourceTree( LoadTree_IC40_CutA6_nugen_numu_8k_small() );
  ark.evLoader.SetSourceZenWidthDeg(ark.sourceZenWidthDeg);


  //
  // SPECIFY CUTS
  //

  TCut ic40_Cut = "1";

  // Cuts can be added or reset any time, and events re-loaded with new cuts
  // (cuts are also evaluated in order, for possible efficiency gain)

  ark.evLoader.AddCut(ic40_Cut);


  //
  // SET NAMES OF VARIABLES TO BE USED FROM ROOT FILES
  //

  TString recoZenRadName = "mZr"; // needed for eProb below
  ark.evLoader.SetName_recoZenith_rad(recoZenRadName);
  ark.evLoader.SetName_recoAzimuth_rad("mAr");
//  ark.evLoader.SetName_sigmaDeg("mpfSigmaDegRescaled2(mpfSigmaDeg,mmueEn)");
  ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC40(mpfSigmaDeg,mmueEn)");
  ark.evLoader.SetName_runID("RunID");
  ark.evLoader.SetName_eventID("EventID");
  // EXAMPLE: How to load different vars for source and bkg, if necessary
  //  dsBkg->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg")
  //  dsSrc->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg*1.4");
  //  evLoader.SetName_sigmaDeg("myRecoSigmaDeg");


  // The choice here has consequences for how to fill eProb, below
  TString energyVar = "log10(mueEn)";

  ark.evLoader.SetName_energyValue(energyVar);


  //
  // LOAD EVENTS, SET UP ENERGY PDFs
  // 

  cout << "Loading Background Events: " << ark.evLoader.GetBkgTree()->GetTitle() << endl;
  cout << "Using cut: " << ark.evLoader.GetCuts().GetTitle() << endl;
  ark.evLoader.LoadBkgEvents(ark.baseEvents);


  ZenithEnergyProb* zen_eProb = new ZenithEnergyProb();

 
  cout << "Filling Energy PDFs:\n";
  zen_eProb->SetSourceZenWidthDeg( ark.sourceZenWidthDeg);
  zen_eProb->SetName_recoZenith_rad(recoZenRadName);

  vector<double> zenMinDegVect;
      
  // Fill vector with bin edges that match "CutDMS" bin edges
  // 0 - 90  zen added
  double tempBot=1;
  for(int i=0; i<21;i++){
    //cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
    zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
    tempBot-=0.05;
  }
  // A reasonable split for the upgoing range
  //zenMinDegVect.push_back(105.00); //? Follow Chad's example
  zenMinDegVect.push_back(110.00); // ~acos(-0.33)
  zenMinDegVect.push_back(132.00); // ~acos(-0.67)
  zenMinDegVect.push_back(180.);

  zen_eProb->SetZenithBandsDeg(zenMinDegVect);

  zen_eProb->SetLoadModeNew(true); // true is now faster

  if (1) { } // new default is to constrain signal inside histogram
  if (0) { zen_eProb->SetConstrainSignal(false); }  // the old way

  if (energyVar == "log10(mueEn)") {
    int nBackFill = 35; // don't backfill previous bins
    zen_eProb->SetEnergyGammaRangeAndBackFill(40,2.,9., 30,1.,4., nBackFill);
  }

  zen_eProb->SetTableBkg(ark.baseEvents);
  zen_eProb->SetTableGamma(ark.evLoader.GetSourceTree(), ark.evLoader.GetCuts(), energyVar);

  ark.eProb = zen_eProb;


  // This seems to be okay for ~ 5000 events or more:
  // 1. deg smoothing and 180 declination maps (e.g. 1 deg binning)
  ark.decBkgProb.Initialize(180, 1);

  ark.decBkgProb.SetBaseDecMap(ark.baseEvents);


  psData = new I3Analysis();

  psData->SetBkgSpaceProb(ark.decBkgProb);
  psData->SetBaseEvents(ark.baseEvents);
  psData->SetEnergyProb(*(ark.eProb));
  
  EventTimeModule * ic40times = new EventTimeModuleDiscrete();

  if (opt_UseRealData) {
    psData->UseRealData();  // the event set is now exactly equal
    // to the data set (i.e. no scrambling, no fake signal added.)
  } else {
    //psData->SetRandomizeBase(true);
    //ic40times->SetTimesFromMJDFile("../ic59/HugeListOfIC40TimesNewMethod1.txt");
    //psData->SetEventTimeModulePtr(ic40times);
    psData->GenerateDataSet_with_nSrcEvents(0); // needs an I3SignalGenerator first??
  } //*/

  ark.psData = psData;

  cout << "----- IC-40 Loaded -----\n";

}
