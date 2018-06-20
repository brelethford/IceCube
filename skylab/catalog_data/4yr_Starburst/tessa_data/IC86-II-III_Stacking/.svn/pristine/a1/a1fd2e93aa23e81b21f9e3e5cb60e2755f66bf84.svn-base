TString LOADTREE  = "/net/user/tcarver/Data/PS/MESE/Data/";//"/home/tcarver/psLab_MESE/files/data/";
//TString LOADTREE_MESE_III  = "/data/IceCube/IC86-III/MESE_no_cut/";

TTree* LoadTree_IC86_nugen_numu_E_1() {
    vector <TString> nufiles;
    nufiles.push_back("/net/user/tcarver/Data/PS/MESE/Nugen/MergedNugenMillipede_3999.root");
    return LoadTree_IC86_nugen_numu_E_1(nufiles);
}
  
TTree* LoadTree_IC86_nugen_numu() {
    TChain *tree  =new TChain("MasterTree");
    TChain *trF =new TChain("MESECut");
    tree->Add("/net/user/tcarver/Data/PS/MESE/Nugen/MergedNugenMillipede_3999.root");
    trF->Add("/net/user/tcarver/Data/PS/MESE/Nugen/MergedNugenMillipede_3999_friend.root");
    tree->AddFriend(trF);
    tree->GetEntries();
//    tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon.energy");
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
//    tree->SetAlias("NChan",                 "HitMultiplicityValuesIC.n_hit_doms");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
//    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
 //   tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
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


TTree* LoadTree_IC79_86_I_to_IV_GoodRuns_Full() {
  cout<<"TEST 1"<<endl;
  
  TString filesIC79 = LOADTREE+"IC79_MESE.root";
  TString files2011 = LOADTREE+"IC86_2011_MESE_j.root";
  TString files2012 = LOADTREE+"IC86_2012_MESE_j.root";
  TString files2013 = LOADTREE+"MESE_noCut_IC86.2013_merged.root";
  TString files2014 = LOADTREE+"MESE_noCut_IC86.2014_merged.root";
 // TString files2013 = LOADTREE+"IC86_2013_MESE.root";  
 // TString files2014 = LOADTREE+"IC86_2014_MESE.root";  
  cout<<"TEST 2"<<endl;  

  TChain *tree = new TChain("MasterTree");
  TChain *trF = new TChain("MESECut");
  tree->Add(filesIC79);
  //trF->Add(LOADTREE+"IC79_MESE_friend.root");
  tree->Add(files2011);
  //trF->Add(LOADTREE+"IC86_2011_MESE_j_friend.root");
  tree->Add(files2012);
  //trF->Add(LOADTREE+"IC86_2012_MESE_j_friend.root");
  tree->Add(files2013);
  //trF->Add(LOADTREE+"IC86_2013_MESE_friend.root");
  tree->Add(files2014);/*
  trF->Add(LOADTREE+"IC86_2014_MESE_friend.root");
  tree->AddFriend(trF);*/
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  cout<<"TEST 3"<<endl;  

  TString livetimeTotalStr = " SplineMPEMuEXDifferential.energy/SplineMPEMuEXDifferential.energy * 148203648.0"; // livetime in sec (726.78+988.54)*86400.0 (II+III+IV=  988.54 days)
  TString tmin = " SplineMPEMuEXDifferential.energy/SplineMPEMuEXDifferential.energy * 55347.3"; //in MJD for ic79
  TString tmax = " SplineMPEMuEXDifferential.energy/SplineMPEMuEXDifferential.energy * 57160.0410"; // max for ic86 IV (for III = 56783.5781)";
  
  cout<<livetimeTotalStr<<endl;
  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal","1.0*148203648.0"); // livetime in sec (II+III+IV=  91371456.0)
  tree->SetAlias("tmin", "1.0*55347.3");//in MJD for ic79
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

  //tree->SetAlias("z","cos(SplineMPE.zenith)");
  cout<<"TEST TREELOADER in"<<endl;
  return tree;
}
