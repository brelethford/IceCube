TString LOADTREE_86_II  = "/data/IceCube/IC86-II/IC86_2012/";
TString LOADTREE_86_III = "/data/IceCube/IC86-III/IC86_2013/";
TString LOADTREE_86_IV = "/data/IceCube/IC86-IV/IC86_2014/";
TString LOADTREE_nugen = "/net/user/tcarver/Data/PS/IC86_2012/";

TTree* LoadTree_IC86_II_III_IV_nugen_numu_E_1() {
    vector <TString> nufiles;
    nufiles.push_back("Merged_11069_upgoing_6089filesIC86.2012.root");
    nufiles.push_back("Merged_11069_downgoing_6089filesIC86.2012.root");
    nufiles.push_back("Merged_11070_upgoing_4829filesIC86.2012.root");
    nufiles.push_back("Merged_11070_downgoing_4829filesIC86.2012.root");
    return LoadTree_IC86_II_III_IV_nugen_numu(nufiles);
}
  
TTree* LoadTree_IC86_II_III_IV_nugen_numu(vector <TString> nufiles) {
    TChain *tree  =new TChain("MasterTree");
    /*
    gSystem->Load("libneutrinoflux");
    gROOT->SetMacroPath(TString(gROOT->GetMacroPath())+"/opt/icetray/meta-projects/icerec/V04-11-02/src/neutrinoflux/public/:");
    gInterpreter->AddIncludePath("/opt/icetray/meta-projects/icerec/V04-11-02/src/neutrinoflux/public/");
    gROOT->ProcessLine(".L neutrinoflux/NeutrinoFlux.h+");
    gROOT->ProcessLine(".L neutrinoflux/NuFluxFn.h+");
    NeutrinoFlux *myFlux0 = new AtmosphericNeutrinoFlux("honda2006_numu","sarcevic_std_numu");
    NuFluxFnSet(myFlux0, 0);
*/
    for (unsigned int j=0;j<nufiles.size();j++){
        tree->Add(LOADTREE_nugen+nufiles[j]);
    }
    
    tree->GetEntries();
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
    //tree->SetAlias("atmoWeight","NuFluxFn(MCPrimary1.pdg_encoding,MCPrimary1.energy,cos(MCPrimary1.zenith),0)*I3MCWeightDict.OneWeight/(I3MCWeightDict.NEvents*3999/2.)");

    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
    tree->SetAlias("SplinePbSigmaDeg", "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status==0) + SplineMPEBootstrapVectStats.median*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status!=0)");
    tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
    tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
    tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",               "I3EventHeader.Event");
    tree->SetAlias("mcOneWeight", "I3MCWeightDict.OneWeight / (I3MCWeightDict.NEvents*((TMath::FloorNint(RunID/100000)==11069)*6089 + (TMath::FloorNint(RunID/100000)==11070)*4829))");

    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
    }
  
    tree->SetAlias("mcTotalGeneratedEvents", "1.*1.");
    return tree;
}


TTree* LoadTree_IC86_II_III_IV_nugen_numu_11029() {
  TString upfilesnugen = LOADTREE_nugen + "Merged_11029_upgoing_5170filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_nugen + "Merged_11029_downgoing_5170filesIC86.2012.root";
  TChain *tree = new TChain("MasterTree");
  tree->Add(upfilesnugen);
  tree->Add(downfilesnugen);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles = 5170;
  
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
  tree->SetAlias("SplinePbSigmaDeg", "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status==0) + SplineMPEBootstrapVectStats.median*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status!=0)");
  tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
  tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
  // changed to JointWeight
  //tree->SetAlias("mcOneWeight",           "JointWeight.value");
  tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",               "I3EventHeader.Event");



  for (int i=0; i<2; i++) {
    tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
  }
  
  TString mcTotalGeneratedEvents = Form("I3MCWeightDict.NEvents * %d",nFiles);


  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  return tree;
}

TTree* LoadTree_IC86_II_III_IV_nugen_numu_11069() {
  TString upfilesnugen = LOADTREE_nugen + "Merged_11069_upgoing_6089filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_nugen + "Merged_11069_downgoing_6089filesIC86.2012.root";
  TChain *tree = new TChain("MasterTree");
  tree->Add(upfilesnugen);
  tree->Add(downfilesnugen);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles = 6089;
  
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
  tree->SetAlias("SplinePbSigmaDeg", "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status==0) + SplineMPEBootstrapVectStats.median*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status!=0)");
  tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
  tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
  // changed to JointWeight
  //tree->SetAlias("mcOneWeight",           "JointWeight.value");
  tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",               "I3EventHeader.Event");



  for (int i=0; i<2; i++) {
    tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
  }
  TString mcTotalGeneratedEvents = Form("I3MCWeightDict.NEvents * %d",nFiles);

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  return tree;
}

TTree* LoadTree_IC86_II_III_IV_nugen_numu_11070() {
  TString upfilesnugen = LOADTREE_nugen + "Merged_11070_upgoing_4829filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_nugen + "Merged_11070_downgoing_4829filesIC86.2012.root";
  TChain *tree = new TChain("MasterTree");
  tree->Add(upfilesnugen);
  tree->Add(downfilesnugen);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles = 4829;
  
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
  tree->SetAlias("SplinePbSigmaDeg", "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status==0) + SplineMPEBootstrapVectStats.median*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status!=0)");
  tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
  tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
  tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",               "I3EventHeader.Event");



  for (int i=0; i<2; i++) {
    tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
  }
  
  TString mcTotalGeneratedEvents = Form("I3MCWeightDict.NEvents * %d",nFiles);

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);

  return tree;
}


TTree* LoadTree_IC86_II_III_IV_GoodRuns_Full() {
  
  TString upfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_upgoing.root";
  TString downfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_downgoing.root";
  
  TString upfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_upgoing.root";
  TString downfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_downgoing.root";  
  
  TString upfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_upgoing.root";
  TString downfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_downgoing.root";  
  
  TChain *tree = new TChain("MasterTree");
  tree->Add(downfiles2012);
  tree->Add(upfiles2012);
  tree->Add(downfiles2013);
  tree->Add(upfiles2013);
  tree->Add(upfiles2014);
  tree->Add(downfiles2014);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  TString livetimeTotalStr = " 1. * 91371456.0"; // livetime in sec (II+III = 59644512.0)
  TString tmin = " 1. * 56062.4207"; //in MJD
  TString tmax = " 1. * 57160.0410"; // for III = 56783.5781";

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", livetimeTotalStr);
  tree->SetAlias("tmin", tmin);
  tree->SetAlias("tmax", tmax);
  
  //--------------------//
  tree->SetAlias("timeMJD","timeMJD.value");
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
  tree->SetAlias("err1","SplineMPEParaboloidFitParams.err1");
  tree->SetAlias("err2","SplineMPEParaboloidFitParams.err2");
  tree->SetAlias("boot_median","SplineMPEBootstrapVectStats.median");
  tree->SetAlias("SplineMPEParaboloidstatus","SplineMPEParaboloidFitParams.status");
  tree->SetAlias("SplinePbSigmaDeg", "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status==0) + SplineMPEBootstrapVectStats.median*TMath::RadToDeg()*(SplineMPEParaboloidFitParams.status!=0)");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",                 "I3EventHeader.Event");


  return tree;
}
