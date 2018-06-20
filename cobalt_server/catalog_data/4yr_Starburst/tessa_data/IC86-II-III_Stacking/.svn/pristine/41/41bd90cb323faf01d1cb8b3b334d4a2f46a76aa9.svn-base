TString LOADTREE_86_II  = "/data/IceCube/IC86-II/";
TString LOADTREE_86_III = "/data/IceCube/IC86-III/";

TTree* LoadTree_IC862plus3_nugen_numu_E_1() {
    vector <TString> nufiles;
    nufiles.push_back("IC86_2012_NuGen_11069_Upgoing_6089files.root");
    nufiles.push_back("IC86_2012_NuGen_11069_Downgoing_6089files.root");
    nufiles.push_back("IC86_2012_NuGen_11070_Upgoing_4829files.root");
    nufiles.push_back("IC86_2012_NuGen_11070_Downgoing_4829files.root");
    return LoadTree_IC862plus3_nugen_numu(nufiles);
}
  
TTree* LoadTree_IC862plus3_nugen_numu(vector <TString> nufiles) {
    TChain *tree  =new TChain("MasterTree");
    TChain *trF =new TChain("NFiles");

    for (unsigned int j=0;j<nufiles.size();j++){
        tree->Add(LOADTREE_86_II+nufiles[j]);
        trF->Add(LOADTREE_86_II+nufiles[j].ReplaceAll(".root","_NfilesTree.root"));
    }
    tree->AddFriend(trF);
    
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
    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
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
  
    tree->SetAlias("mcTotalGeneratedEvents", "I3MCWeightDict.NEvents * NFiles.NFiles");
    return tree;
}


TTree* LoadTree_IC862plus3_nugen_numu_11029() {
  TString upfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11029_Upgoing_5170files.root";
  TString downfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11029_Downgoing_5170files.root";
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
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
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
  //let's see if this method of combining datasets works...
  //TString mcTotalGeneratedEvents = "NEvents.value/NEvents.value * 1.0";

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));
  
  //tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  return tree;
}

TTree* LoadTree_IC862plus3_nugen_numu_11069() {
  TString upfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11069_Upgoing_6089files.root";
  TString downfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11069_Downgoing_6089files.root";
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
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
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
  //let's see if this method of combining datasets works...
  //TString mcTotalGeneratedEvents = "NEvents.value/NEvents.value * 1.0";

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));
  
  //tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  return tree;
}

TTree* LoadTree_IC862plus3_nugen_numu_11070() {
  TString upfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11070_Upgoing_4829files.root";
  TString downfilesnugen = LOADTREE_86_II + "IC86_2012_NuGen_11070_Downgoing_4829files.root";
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
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
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
  //let's see if this method of combining datasets works...
  //TString mcTotalGeneratedEvents = "NEvents.value/NEvents.value * 1.0";

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));
  
  //tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  return tree;
}


TTree* LoadTree_IC862plus3_GoodRuns_Full() {
  
  TString upfiles2012 = LOADTREE_86_II+"IC86_2012_Data_Upgoing.root";
  TString downfiles2012 = LOADTREE_86_II+"IC86_2012_Data_Downgoing.root";
  
  TString upfiles2013 = LOADTREE_86_III+"IC86_2013_Data_Upgoing.root";
  TString downfiles2013 = LOADTREE_86_III+"IC86_2013_Data_Downgoing.root";  
  
  TChain *tree = new TChain("MasterTree");
  tree->Add(downfiles2012);
  tree->Add(upfiles2012);
  tree->Add(downfiles2013);
  tree->Add(upfiles2013);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  TString livetimeTotalStr = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 59644512.0"; // livetime in sec
  TString tmin = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 56062.4207"; //in MJD
  TString tmax = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 56783.5781";

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
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",                 "I3EventHeader.Event");

  //tree->SetAlias("z","cos(SplineMPE.zenith)");

  return tree;
}
