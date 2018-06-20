void load_ark_IC79_to_86IV_MESE_tdep_lc_folded(I3Ark& ark, bool opt_UseRealData, TString RecoName, double NbinsZen=12., double NbinsAz=90.) {
    cout << "\n-----  LOADING MESE IC-79, 86 1 to 4 with folded lcBkgProb-----\n";
    ark.evLoader.SetBkgTree( LoadTree_IC79_86_I_to_IV_MESEGoodRuns_Full() );
  
    cout << "Loaded tree MESEGoodRuns" << endl;
    
    ark.livetime = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");
    ark.tmin = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmin");
    ark.tmax = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmax");
    
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
  
    cout << "Configuring Source Event Sample...\n";

    ark.sourceZenWidthDeg = 0.5;

    ark.evLoader.SetSourceTree( LoadTree_IC86_nugen_numu() );
    ark.evLoader.SetSourceZenWidthDeg(ark.sourceZenWidthDeg);
    
    TCut ICMESE_Cuta = "TMath::Min(Millipede_FirstLoss_XYDist.value, Millipede_FirstLoss_ZDist.value) > -81.*TMath::Log10(SplineMPEMuEXDifferential.energy) + 426.";
    TCut ICMESE_Cutb = "(acos(sin(LineFit.zenith)*sin(SplineMPE.zenith)*cos(LineFit.azimuth - SplineMPE.azimuth)+cos(LineFit.zenith)*cos(SplineMPE.zenith))*TMath::RadToDeg()<10.**1.62 )";
    TCut ICMESE_Cutc  = "(SplineMPE.zenith*TMath::RadToDeg()<85.)";
    TCut ICMESE_86_3and4 = "timeMJD>56414.416";
    TCut ICMESE_nugen = "timeMJD==55744.986";
 
    TCut MESE_Cut =  ICMESE_Cuta&&ICMESE_Cutb&&ICMESE_Cutc;
    //TCut MESE_Cut =  ICMESE_Cutb&&ICMESE_Cutc;
    ark.evLoader.AddCut(MESE_Cut);

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
        ark.evLoader.SetName_sigmaDeg("RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(SplinePbSigmaDeg,muexEn)");
    }
    else cout <<"AHHHHHHHHHHHHHHHHHHHH";

    ark.evLoader.SetName_runID("RunID");
    ark.evLoader.SetName_eventID("EventID");

    TString energyVar = "log10(muexEn)";

    ark.evLoader.SetName_energyValue(energyVar);

    // LOAD EVENTS, SET UP ENERGY PDFs
    cout << "Loading Background Events: " << ark.evLoader.GetBkgTree()->GetTitle() << endl;

    ark.evLoader.SetMonitor(true);
    ark.evLoader.LoadBkgEvents(ark.baseEvents);

    ZenithEnergyProb* zen_eProb = new ZenithEnergyProb();
 
    cout << "Filling Energy PDFs:\n";
    zen_eProb->SetSourceZenWidthDeg( ark.sourceZenWidthDeg);
    zen_eProb->SetName_recoZenith_rad(recoZenRadName);

    vector<double> zenMinDegVect;
    //reduced energy binning according to Jon Dumm
    zenMinDegVect.push_back(0.0);
    zenMinDegVect.push_back(86.0);
    zenMinDegVect.push_back(130.0);
    zenMinDegVect.push_back(180.00);
      
    zen_eProb->SetZenithBandsDeg(zenMinDegVect);

    zen_eProb->SetLoadModeNew(true); // true is now faster

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

    // This seems to be okay for ~ 5000 events or more:
    // 1. deg smoothing and 180 declination maps (e.g. 1 deg binning)
    
    ark.decBkgProb.Initialize(180, 1);
    ark.decBkgProb.SetBaseDecMap(ark.baseEvents);
  
    ark.lcBkgProb.Initialize(NbinsZen, NbinsAz,true); //so why doesn't this work together?
    ark.lcBkgProb.FillLCBkgHisto_folded(ark.baseEvents,0.,90.);
          
    psData = new I3Analysis();
    psData->SetBkgSpaceProb(ark.decBkgProb);
    psData->SetBaseEvents(ark.baseEvents);
    psData->SetEnergyProb(*(ark.eProb));

  
    EventTimeModule * mese_times = new EventTimeModuleDiscrete();

    if (opt_UseRealData) {
        psData->UseRealData();  // the event set is now exactly equal
        // to the data set (i.e. no scrambling, no fake signal added.)
    } else {
        psData->SetRandomizeBase(true);
        mese_times->SetTimesFromMJDFile("/home/christov/psLab/macro_llh/IC86-IV_TDep/");
        psData->SetEventTimeModulePtr(mese_times);
        psData->GenerateDataSet_with_nSrcEvents(0); // needs an I3SignalGenerator first??
    } 
  
 
    ark.psData = psData;
    
    cout << "----- IC86 2 to 4  Loaded -----\n";
}
