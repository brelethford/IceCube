TString LOADTREE_IC40_FINAL_DIR = "/net/user/mrameez/trunk/rootfiles/IC40/";//"/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/";

// Original signal dataset for analysis (L2b) [1000 E^-1 files, 5000 events/file]
// Now 10x more statistics available!!!
TTree* LoadTree_IC40_CutA6_nugen_numu_2075_small() {

  char *files = LOADTREE_IC40_FINAL_DIR + "/CutA6_v1_nugen_numu_2075_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; 
  //int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  //  int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;


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

  TString mcTotalGeneratedEvents = " NEvents * " + TStringify(nFiles);

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);


  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}

// Recommended signal dataset for analysis (L2b) [1000 E^-1 files, 5000 events/file]
// 10x the statistics of original unblinding sample
TTree* LoadTree_IC40_CutA6_nugen_numu_8k_small() {

  char *file1 = LOADTREE_IC40_FINAL_DIR + "/CutA6_v1_nugen_numu_2075_01000_small.root";
  char *file2 = LOADTREE_IC40_FINAL_DIR + "/CutA6_v1_nugen_numu_2082_01000_small.root";
  char *file3 = LOADTREE_IC40_FINAL_DIR + "/CutA6_v1_nugen_numu_3311_08000_small.root";

  TChain *tree = new TChain("Level2");
  //tree->Add(file1);
  //tree->Add(file2);
  tree->Add(file3);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 8000; 
  //int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
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


  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}

// For Testing Only: 90% DOM Eff (L2b) [1000 E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_nugen_numu_3264_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_3264_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; 
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

// For Testing Only: 110% DOM Eff (L2b) [1000 E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_nugen_numu_3265_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_3265_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; 
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

// 
TTree* LoadTree_IC40_CutA6_nugen_numu__small() {

  char *files = "";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9775; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFF1
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFF1_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFF1_09993_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9993; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP70
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP70_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP70_09950_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9950; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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
// SPICE2-PLUS-EFFP85
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP85_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP85_08394_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 8394; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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
// SPICE2-PLUS-EFFP93
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP93_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP93_08394_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 8394; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP95
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP95_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP95_09781_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9781; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP95-SAM10
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP95_SAM10_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP95-SAM10_09517_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9517; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP95-SAMP10
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP95_SAMP10_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP95-SAMP10_09837_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9837; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP95-SAP10
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP95_SAP10_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP95-SAP10_09839_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9839; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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

// SPICE2-PLUS-EFFP95-SAPM10
TTree* LoadTree_IC40_CutA6_nugen_numu_SPICE2_PLUS_EFFP95_SAPM10_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2972_SPICE2-PLUS-EFFP95-SAPM10_09835_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 9835; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

  tree->SetAlias("mpfSigmaDeg",           "mpfSigmaDeg/sqrt(sqrt(2.))");

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


// For testing (L2b) [5000 E^-2 files, 200000 events/file]
TTree* LoadTree_IC40_CutA6_nugen_numu_2326_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2326_05000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 5000; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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
// For testing ONLY (L2b) 90% DOM efficiency [2000 E^-2 files, 200000 events/file]
TTree* LoadTree_IC40_CutA6_nugen_numu_3247_small() {

  cout << "\n\n WARNING: You are loading a test sample!!!\n\n";

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_3247_02000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 2000; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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
// For testing ONLY: 110% DOM efficiency (L2b) [2000 E^-2 files, 200000 events/file]
TTree* LoadTree_IC40_CutA6_nugen_numu_3248_small() {

  cout << "\n\n WARNING: You are loading a test sample!!!\n\n";

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_3248_02000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 2000; 
  int generatedEventsPerFile = 200000; // NuMu + NuMuBar per file
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

// For testing ONLY (L2a) [1000 E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_L2a_nugen_numu_1882_small() {

  cout << "\n\n WARNING: You are loading a test sample!!!\n\n";

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v2/CutA6_v2_Nugen1882_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; 
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

// For testing ONLY (L2b) [1000 E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_L2b_nugen_numu_1882_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_1882_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; 
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


// Same time period as IC40 6-month dataset, but L2 reprocessed
TTree* LoadTree_IC40_CutA6_data_6months_GoodRuns_small() {

  char *filespec10 = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_6months_GoodRuns_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  double livetimeTotal = 175.55 * 86400.; // seconds

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", TStringify(livetimeTotal));

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

// Full IC40 dataset, L2_fix reprocessed
// CutA6 introduced the "level2_fix" reprocessing
TTree* LoadTree_IC40_CutA6_data_Fix_GoodRuns_small() {

  char *filespec10 = LOADTREE_IC40_FINAL_DIR+"/CutA6_v1_Full_GoodRuns_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  double livetimeTotal = 32446569.599999998; //375.539 * 86400.; // seconds

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", "32446569.599999998 * NCh / NCh");

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

// Original 1-month burn sample for IC40 dataset, L2 reprocessed
TTree* LoadTree_IC40_CutA6_data_Burn_GoodRuns_small() {

   char *filespec10 = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_Burn_GoodRuns_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  double livetimeTotal = 29.79 * 86400.; // seconds

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", TStringify(livetimeTotal));

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

// Original signal dataset for analysis with NuMu:NuTau = 1:1 
// Now 10x more statistics available
//  (L2b) [NuMu:  1000 E^-1 files, 5000 events/file]
//  (L2b) [NuTau: 5000 E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_numu_2075_plus_nutau_2183_small() {

 char *file1 = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2075_01000_small.root";
 char *file2 = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_nutau_2183_05000_small.root";

 TChain *tree = new TChain("Level2");
 tree->Add(file1);
 tree->Add(file2);
 tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


 // SET DATASET SELECTORS:

 TString DataSetSelector_NuMu  = "(Type==68 || Type==69)";
 TString DataSetSelector_NuTau = "(Type==133 || Type==134)";


 // These default scale factors can be changed at a later time by the user.
 tree->SetAlias("NuMuFluxScaleFactor", "1.");
 tree->SetAlias("NuTauFluxScaleFactor","1.");


 // SET META-INFORMATION ALIAS(ES):

 int NuMu_nFiles = 1000;
 int NuMu_generatedEventsPerFile = 5000; // Nu + NuBar per file
 int NuMu_mcTotalGeneratedEvents = NuMu_nFiles * NuMu_generatedEventsPerFile;

 int NuTau_nFiles = 5000; 
 int NuTau_generatedEventsPerFile = 5000; // Nu + NuBar per file
 int NuTau_mcTotalGeneratedEvents=NuTau_nFiles * NuTau_generatedEventsPerFile;

 tree->SetAlias("mcOneWeight",
		 "("+
		 DataSetSelector_NuMu+"*(NuMuFluxScaleFactor*OW/"+
		 TStringify(NuMu_mcTotalGeneratedEvents)+
		 ")+"+
		 DataSetSelector_NuTau+"*(NuTauFluxScaleFactor*OW/"+
		 TStringify(NuTau_mcTotalGeneratedEvents)+
		 "))"
		 );

 tree->SetAlias("mcTotalGeneratedEvents","1."); // Hidden in mcOneWeight now

 // SET HANDLES:
 // These aliases provide consistent handles across different 
 // types root files (analysis-tree, flat-ntuple, simple-extractor),
 // so set accordingly:

 tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
 tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
 tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");


 //--------------------//

 // OTHER, EXTRA ALIASES:

 tree->SetAlias("mueEn",                 "mmueEn");
 tree->SetAlias("NChan",                 "NCh");
 tree->SetAlias("s32Zr","s32Zd*TMath::DegToRad()");
 tree->SetAlias("s32Ar","s32Ad*TMath::DegToRad()");
 tree->SetAlias("mZr","mZd*TMath::DegToRad()");
 tree->SetAlias("mAr","mAd*TMath::DegToRad()");
 //  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");
 //  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

 return tree;
}

// Recommended signal dataset for analysis with NuMu:NuTau = 1:1 
//  (L2b) [NuMu:  10000 E^-1 files, 5000 events/file]
//  (L2b) [NuTau: 5000  E^-1 files, 5000 events/file]
TTree* LoadTree_IC40_CutA6_numu_10k_plus_nutau_2183_small() {

 char *numuFile1 = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2075_01000_small.root";
 char *numuFile2 = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_2082_01000_small.root";
 char *numuFile3 = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_numu_3311_08000_small.root";
 char *nutauFile = 
   "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA6_v1/CutA6_v1_nugen_nutau_2183_05000_small.root";

 TChain *tree = new TChain("Level2");
 tree->Add(numuFile1);
 tree->Add(numuFile2);
 tree->Add(numuFile3);
 tree->Add(nutauFile);
 tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


 // SET DATASET SELECTORS:

 TString DataSetSelector_NuMu  = "(Type==68 || Type==69)";
 TString DataSetSelector_NuTau = "(Type==133 || Type==134)";


 // These default scale factors can be changed at a later time by the user.
 tree->SetAlias("NuMuFluxScaleFactor", "1.");
 tree->SetAlias("NuTauFluxScaleFactor","1.");


 // SET META-INFORMATION ALIAS(ES):

 int NuMu_nFiles = 10000;
 int NuMu_generatedEventsPerFile = 5000; // Nu + NuBar per file
 int NuMu_mcTotalGeneratedEvents = NuMu_nFiles * NuMu_generatedEventsPerFile;

 int NuTau_nFiles = 5000; 
 int NuTau_generatedEventsPerFile = 5000; // Nu + NuBar per file
 int NuTau_mcTotalGeneratedEvents=NuTau_nFiles * NuTau_generatedEventsPerFile;

 tree->SetAlias("mcOneWeight",
		 "("+
		 DataSetSelector_NuMu+"*(NuMuFluxScaleFactor*OW/"+
		 TStringify(NuMu_mcTotalGeneratedEvents)+
		 ")+"+
		 DataSetSelector_NuTau+"*(NuTauFluxScaleFactor*OW/"+
		 TStringify(NuTau_mcTotalGeneratedEvents)+
		 "))"
		 );

 tree->SetAlias("mcTotalGeneratedEvents","1."); // Hidden in mcOneWeight now

 // SET HANDLES:
 // These aliases provide consistent handles across different 
 // types root files (analysis-tree, flat-ntuple, simple-extractor),
 // so set accordingly:

 tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
 tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
 tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");


 //--------------------//

 // OTHER, EXTRA ALIASES:

 tree->SetAlias("mueEn",                 "mmueEn");
 tree->SetAlias("NChan",                 "NCh");
 tree->SetAlias("s32Zr","s32Zd*TMath::DegToRad()");
 tree->SetAlias("s32Ar","s32Ad*TMath::DegToRad()");
 tree->SetAlias("mZr","mZd*TMath::DegToRad()");
 tree->SetAlias("mAr","mAd*TMath::DegToRad()");
 //  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");
 //  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

 return tree;
}



