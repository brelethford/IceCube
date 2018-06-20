TString LOADTREE_86_II  = "/net/user/tcarver/Data/PS/IC86_2012/";
TString LOADTREE_86_III = "/net/user/tcarver/Data/PS/IC86_2013/";
TString LOADTREE_86_IV = "/net/user/tcarver/Data/PS/IC86_2014/";
TString LOADTREE_86_IV_tmp = "/net/user/tcarver/Data/PS/IC86-IV_root/";


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
    //TChain *trF =new TChain("NFiles");

    for (unsigned int j=0;j<nufiles.size();j++){
        tree->Add(LOADTREE_86_II+nufiles[j]);
    //    trF->Add(LOADTREE_86_II+nufiles[j].ReplaceAll(".root","_NfilesTree.root"));
    }
  //  tree->AddFriend(trF);
    
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
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    //tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
    tree->SetAlias("EventID",               "I3EventHeader.Event");
    tree->SetAlias("mcOneWeight", "I3MCWeightDict.OneWeight / (I3MCWeightDict.NEvents*((TMath::FloorNint(RunID/100000)==11069)*6089 + (TMath::FloorNint(RunID/100000)==11070)*4829))");


    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
    }
  
    //tree->SetAlias("mcTotalGeneratedEvents", "I3MCWeightDict.NEvents * NFiles.NFiles");
    tree->SetAlias("mcTotalGeneratedEvents", "1.*1.");
    return tree;
}


TTree* LoadTree_IC86_II_III_IV_nugen_numu_11029() {
  TString upfilesnugen = LOADTREE_86_II + "Merged_11029_upgoing_5170filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_86_II + "Merged_11029_downgoing_5170filesIC86.2012.root";
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

TTree* LoadTree_IC86_II_III_IV_nugen_numu_11069() {
  TString upfilesnugen = LOADTREE_86_II + "Merged_11069_upgoing_6089filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_86_II + "Merged_11069_downgoing_6089filesIC86.2012.root";
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

TTree* LoadTree_IC86_II_III_IV_nugen_numu_11070() {
  TString upfilesnugen = LOADTREE_86_II + "Merged_11070_upgoing_4829filesIC86.2012.root";
  TString downfilesnugen = LOADTREE_86_II + "Merged_11070_downgoing_4829filesIC86.2012.root";
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


TTree* LoadTree_IC86_II_III_IV_GoodRuns_Full() {
  
  TString upfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_upgoing.root";
  TString downfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_downgoing.root";
  
  TString upfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_upgoing.root";
  TString downfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_downgoing.root";  
  
  //TString totalfiles2014 = LOADTREE_86_IV+"IC86_2014_total.root";  
  TString upfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_upgoing.root";
  TString downfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_downgoing.root";  
  
  TChain *tree = new TChain("MasterTree");
  tree->Add(downfiles2012);
  tree->Add(upfiles2012);
  tree->Add(downfiles2013);
  tree->Add(upfiles2013);
  //tree->Add(LOADTREE_86_IV_tmp+"*.root");
  tree->Add(downfiles2014);
  tree->Add(upfiles2014);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  TString livetimeTotalStr = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 91371456.0"; // livetime in sec (II+III = 59644512.0)
  TString tmin = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 56062.4207"; //in MJD
  TString tmax = " HitMultiplicityValuesIC.n_hit_doms/HitMultiplicityValuesIC.n_hit_doms * 57160.0410"; // for III = 56783.5781";

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
