TString LOADTREE  = "/net/user/tcarver/Data/PS/MESE/Data/";
  
TTree* LoadTree_IC86_nugen_numu() {
    /*
    gSystem->Load("libneutrinoflux");
    gROOT->SetMacroPath(TString(gROOT->GetMacroPath())+"/opt/icetray/meta-projects/icerec/V04-11-02/src/neutrinoflux/public/:");
    gInterpreter->AddIncludePath("/opt/icetray/meta-projects/icerec/V04-11-02/src/neutrinoflux/public/");
    gROOT->ProcessLine(".L neutrinoflux/NeutrinoFlux.h+");
    gROOT->ProcessLine(".L neutrinoflux/NuFluxFn.h+");
    NeutrinoFlux *myFlux0 = new AtmosphericNeutrinoFlux("honda2006_numu","sarcevic_std_numu");
    NuFluxFnSet(myFlux0, 0);
    */
    
    TChain *tree  =new TChain("MasterTree");
    TChain *trF =new TChain("MESECut");
    tree->Add("/net/user/tcarver/Data/PS/MESE/Nugen/MergedNugenMillipede_3999.root");
    tree->GetEntries();
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
    
    //tree->SetAlias("atmoWeight","NuFluxFn(MCPrimary1.pdg_encoding,MCPrimary1.energy,cos(MCPrimary1.zenith),0)*I3MCWeightDict.OneWeight/(I3MCWeightDict.NEvents*3999/2.)");

    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
    tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
    tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
    tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
    tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",               "I3EventHeader.Event");

    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
    }
  
    tree->SetAlias("mcTotalGeneratedEvents", "I3MCWeightDict.NEvents * 3999");
    return tree;
}


TTree* LoadTree_IC79_86_I_to_IV_MESEGoodRuns_Full() {
  
    TString filesIC79 = LOADTREE+"IC79_MESE.root";
    TString files2011 = LOADTREE+"IC86_2011_MESE_j.root";
    TString files2012 = LOADTREE+"IC86_2012_MESE_j.root";
    TString files2013 = LOADTREE+"MESE_noCut_IC86.2013_merged.root";
    TString files2014 = LOADTREE+"MESE_noCut_IC86.2014_merged.root";

    TChain *tree = new TChain("MasterTree");
    tree->Add(filesIC79);
    tree->Add(files2011);
    tree->Add(files2012);
    tree->Add(files2013);
    tree->Add(files2014);
    tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  
    // SET META-INFORMATION:
    // These aliases store meta-information for the tree
    
    tree->SetAlias("livetimeTotal","1.0*148203648.0"); // livetime in sec (II+III+IV=  91371456.0)
    tree->SetAlias("tmin", "1.0* 55347.3");//in MJD for ic79
    tree->SetAlias("tmax", "1.0* 57160.0410"); // max for ic86 IV (for III = 56783.5781)";
    
    //--------------------//
    tree->SetAlias("timeMJD","I3EventHeader.time_start_mjd_day + (I3EventHeader.time_start_mjd_sec + I3EventHeader.time_start_mjd_ns*1.e-9)/86400.");
    // OTHER, EXTRA ALIASES:
        
    tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon.energy");
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
    tree->SetAlias("NChan",                 "HitMultiplicityValuesIC.n_hit_doms");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
    tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",                 "I3EventHeader.Event");

    return tree;
}

TTree* LoadTree_86_II_to_IV_MESEGoodRuns_Full() {
  

    TString files2012 = LOADTREE+"IC86_2012_MESE_j.root";
    TString files2013 = LOADTREE+"MESE_noCut_IC86.2013_merged.root";
    TString files2014 = LOADTREE+"MESE_noCut_IC86.2014_merged.root";

    TChain *tree = new TChain("MasterTree");

    tree->Add(files2012);
    tree->Add(files2013);
    tree->Add(files2014);
    tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  
    // SET META-INFORMATION:
    // These aliases store meta-information for the tree
    
    tree->SetAlias("livetimeTotal","1.0*91371456.0"); // livetime in sec (II+III+IV=  91371456.0)
    tree->SetAlias("tmin", "1.0* 56062.4207");//in MJD for ic79
    tree->SetAlias("tmax", "1.0* 57160.0410"); // max for ic86 IV (for III = 56783.5781)";
    
    //--------------------//
    tree->SetAlias("timeMJD","I3EventHeader.time_start_mjd_day + (I3EventHeader.time_start_mjd_sec + I3EventHeader.time_start_mjd_ns*1.e-9)/86400.");
    // OTHER, EXTRA ALIASES:
        
    tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon.energy");
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
    tree->SetAlias("NChan",                 "HitMultiplicityValuesIC.n_hit_doms");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
    tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",                 "I3EventHeader.Event");

    return tree;
}

