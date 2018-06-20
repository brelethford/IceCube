TString LOADTREE_IC59_FINAL_DIR = "/atlas/users/christov/data/IC59";
TString LOADTREE_IC59_FINAL_SC_DIR = "/atlas/users/christov/data/IC59";

TTree* LoadTree_IC59_nugen_numu_6471_small() {
  
  TString file1 = LOADTREE_IC59_FINAL_DIR+"/final.nugen_numu_9.6kf_6471_E2_E27.root";
  
  TChain *tree = new TChain("Level3");

  tree->Add(file1);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles =  96 * 100;
  //int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  // int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;

  
  // SET HANDLES:
  // These aliases provide consistent handles across different 
  // types root files (analysis-tree, flat-ntuple, simple-extractor, small-tree),
  // so set accordingly:

  tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
  tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "mmueEn");
  tree->SetAlias("NChan",                 "NCh");

  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  TString mcTotalGeneratedEvents = "NEvents * " + TStringify(nFiles);
  


  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  tree->SetAlias("NFiles",TStringify(nFiles));

  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}


// Recommended signal dataset for analysis (L2b):
// [10k E^-1 files, 5000 events/file]
// 10x the statistics of original unblinding sample
TTree* LoadTree_IC59_nugen_numu_4175_small() {
  
  TString file1 = LOADTREE_IC59_FINAL_DIR+"/final.nugen_numu_93f_4175_E2_E27_v4.root";
  
  TChain *tree = new TChain("Level3");

  tree->Add(file1);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles =  9300;
  //  int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  //int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;

  
  // SET HANDLES:
  // These aliases provide consistent handles across different 
  // types root files (analysis-tree, flat-ntuple, simple-extractor, small-tree),
  // so set accordingly:

  tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
  tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "mmueEn");
  tree->SetAlias("NChan",                 "NCh");

  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree
  
  TString mcTotalGeneratedEvents = "NEvents * " + TStringify(nFiles);

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  tree->SetAlias("NFiles",TStringify(nFiles));

  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}


TTree* LoadTree_IC59_nugen_numu_4175_SC_small() {
  
  TString file1 = LOADTREE_IC59_FINAL_SC_DIR+"/prunedNugen4175_cutB4_IT2B_new_small.root";
  
  TChain *tree = new TChain("Level3");

  tree->Add(file1);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles =  100 * 100;
  int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;

  
  // SET HANDLES:
  // These aliases provide consistent handles across different 
  // types root files (analysis-tree, flat-ntuple, simple-extractor, small-tree),
  // so set accordingly:

  tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
  tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "mmueEn");
  tree->SetAlias("NChan",                 "NCh");

  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("mcTotalGeneratedEvents", TStringify(mcTotalGeneratedEvents));
  tree->SetAlias("NFiles",TStringify(nFiles));

  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}

TTree* LoadTree_IC59_GoodRuns_final_small() {
  
  TString filespec10 = LOADTREE_IC59_FINAL_DIR+"/final.GoodRunList_E2_E27_v4.root";
  
  TChain *tree = new TChain("Level3");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  
  //  double livetimeTotal = 2930922.36; // seconds

  double livetimeTotal = 30079165.2; //seconds

  double tmin = 54971.154549; //in MJD
  double tmax = 55347.284392;

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", "30079165.2 * NCh/NCh");
  tree->SetAlias("tmin", "54971.154549  * NCh/NCh");
  tree->SetAlias("tmax", "55347.284392 * NCh/NCh");
  //--------------------//

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");

  tree->SetAlias("NChan","NCh");
  tree->SetAlias("mueEn","mmueEn");
  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;
}


TTree* LoadTree_IC59_GoodRuns_final_SC_small() {
  
  TString filespec10 = LOADTREE_IC59_FINAL_SC_DIR+"/IC59_grl_small_CutB4_IT2B_new.root";
  TChain *tree = new TChain("Level3");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  
  //  double livetimeTotal = 2930922.36; // seconds

  double livetimeTotal = 30079165.2; //seconds

  double tmin = 54971.154549; //in MJD
  double tmax = 55347.284392;

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", TStringify(livetimeTotal));

  tree->SetAlias("tmin", "54971.154549  * NCh/NCh");                                                                                                                                                                                        
  tree->SetAlias("tmax", "55347.284392 * NCh/NCh");      
  //--------------------//

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");

  tree->SetAlias("NChan","NCh");
  tree->SetAlias("mueEn","mmueEn");
  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;
}

TTree* LoadTree_IC59_nugen_numu_4175_onetenth_small() {

  cout << "\n Only using 1/10th the full simulation statistics! \n\n";

  //TString file1 = LOADTREE_IC59_FINAL_DIR+"/final.nugen_numu_93f_4175_E2_E27_v3.root";
  //TString file1 = "/net/user/aguilar/work/IceCube/scripts/IC59/final.nugen_numu_93f_4175_E2_E27_v4.root";
  TString file1 = "/net/user/aguilar/work/IceCube/scripts/IC59/final.nugen_numu_10f_4175_E2_E27.root";
  TChain *tree = new TChain("Level3");

  tree->Add(file1);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles =  10 * 100;
  int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;


  // SET HANDLES:
  // These aliases provide consistent handles across different 
  // types root files (analysis-tree, flat-ntuple, simple-extractor, small-tree),
  // so set accordingly:

  tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
  tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "mmueEn");
  tree->SetAlias("NChan",                 "NCh");

  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("mcTotalGeneratedEvents", TStringify(mcTotalGeneratedEvents));


  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}

