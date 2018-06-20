TString LOADTREE_86DIR = "/atlas/users/christov/data/IC86-I/";


TTree* LoadTree_IC86_nugen_numu_joint() {
  TString upfilesnugen = LOADTREE_86DIR + "PrunedNugenUpgoing4866.root";
  TString downfilesnugen = LOADTREE_86DIR + "PrunedNugenDowngoing4866.root";
  TChain *tree = new TChain("MasterTree");
  tree->Add(upfilesnugen);
  tree->Add(downfilesnugen);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles = 4866;
  
  tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon.energy");
  tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
  tree->SetAlias("NChan",                 "HitMultiplicityValues.n_hit_doms");
  string recos[3]={"MPEFit_TT","SplineMPE","MuEXAngular4"};
  for (int i=0; i<3; i++) {
    tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
    tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
  tree->SetAlias("mpbSigmaDeg",           "PbSigmaDeg.value");
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



  for (int i=0; i<3; i++) {
    tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
  }
  
    TString mcTotalGeneratedEvents = Form("NEvents.value * %d",nFiles);
  //let's see if this method of combining datasets works...
  //TString mcTotalGeneratedEvents = "NEvents.value/NEvents.value * 1.0";

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));
  
  //tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  return tree;
}

TTree* LoadTree_IC86_nugen_numu_9366() {
  TString files = LOADTREE_86DIR + "MonteCarlo/9366/BookedV2/Trained*.root";
  TChain *tree = new TChain("MasterTree");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles = 506;
  
  tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon.energy");
  tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
  tree->SetAlias("NChan",                 "HitMultiplicityValues.n_hit_doms");
  string recos[3]={"MPEFit_TT","SplineMPE","MuEXAngular4"};
  for (int i=0; i<3; i++) {
    tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
    tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
  tree->SetAlias("mpbSigmaDeg",           "PbSigmaDeg.value");
  tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
  tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
  tree->SetAlias("SplinePbSigmaDeg",      "SplinePbSigmaDeg.value");
  tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
  tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
  tree->SetAlias("mcOneWeight",           "I3MCWeightDict.OneWeight");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",               "I3EventHeader.Event");
 
  TString mcTotalGeneratedEvents = Form("NEvents.value * %d",nFiles);
  
  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));
  
  //tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  return tree;
}

TTree* LoadTree_IC86_GoodRuns_Burn() {
  
  //TString files = LOADTREE_86DIR+"data/BookedV2/Burn/Trained*.root";
  TString files = LOADTREE_86DIR+"data/Downgoing/Score/NoBundle/Trained*.root";
  
  TChain *tree = new TChain("MasterTree");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  TString livetimeTotalStr = " 1. * 2728913.2"; // livetime in sec
  TString tmin = " 1. * 55695.661898"; //in MJD
  TString tmax = " 1. * 56061.487604";

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", livetimeTotalStr);
  tree->SetAlias("tmin", tmin);
  tree->SetAlias("tmax", tmax);
  
  //--------------------//

  // OTHER, EXTRA ALIASES:
  
  tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon.energy");
  tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
  tree->SetAlias("NChan",                 "HitMultiplicityValues.n_hit_doms");
  string recos[3]={"MPEFit_TT","SplineMPE","MuEXAngular4"};
  for (int i=0; i<3; i++) {
    tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
    tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
  tree->SetAlias("mpbSigmaDeg",           "PbSigmaDeg.value");
  tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
  tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",                 "I3EventHeader.Event");

  //tree->SetAlias("z","cos(SplineMPE.zenith)");

  return tree;
}

TTree* LoadTree_IC86_GoodRuns_Full() {
  
  TString upfiles = LOADTREE_86DIR+"PrunedDataUpgoing_fixed.root";
  TString downfiles = LOADTREE_86DIR+"PrunedDataDowngoing.root";
  
  TChain *tree = new TChain("MasterTree");
  tree->Add(downfiles);
  tree->Add(upfiles);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  TString livetimeTotalStr = " 1. * 28737504.0"; // livetime in sec
  TString tmin = " 1. * 55694.99085730"; //in MJD
  TString tmax = " 1. * 56062.41792913";

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", livetimeTotalStr);
  tree->SetAlias("tmin", tmin);
  tree->SetAlias("tmax", tmax);
  
  //--------------------//
  tree->SetAlias("timeMJD","timeMJD.value");
  // OTHER, EXTRA ALIASES:
  
  tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon.energy");
  tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
  tree->SetAlias("NChan",                 "HitMultiplicityValues.n_hit_doms");
  string recos[3]={"MPEFit_TT","SplineMPE","MuEXAngular4"};
  for (int i=0; i<3; i++) {
    tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
    tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
    tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
  tree->SetAlias("mpbSigmaDeg",           "PbSigmaDeg.value");
  tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
  tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
  tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
  tree->SetAlias("RunID",                 "I3EventHeader.Run");
  tree->SetAlias("EventID",                 "I3EventHeader.Event");

  //tree->SetAlias("z","cos(SplineMPE.zenith)");

  return tree;
}
