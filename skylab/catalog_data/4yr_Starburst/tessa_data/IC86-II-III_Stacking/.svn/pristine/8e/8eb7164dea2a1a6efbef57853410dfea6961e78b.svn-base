TString LOADTREE_DIR = "/net/user/aguilar/work/IceCube/scripts/IC79/";



TTree* LoadTree_IC79_nugen_numu_6308_Sirin_small() {
  TString file1 = "/net/user/sodrowski/IC79NeutrinoLevelData/sim/L6_up1.1_downv4.2_DataSet_006308_9981Files_v13.root";
  TChain *tree = new TChain("CutTree");
  tree->Add(file1);

  int nFiles = 9982;

  tree->SetAlias("mcPrimary_Zenith_rad",  "MC_Zd * 1.74532925199432955e-02");

  tree->SetAlias("mcPrimary_Azimuth_rad", "MC_Ad * 1.74532925199432955e-02");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MC_En");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "MuE_FINAL_MUON_En");
  tree->SetAlias("muexEn",                "MuEX_FINAL_MUON_En");
  tree->SetAlias("NChan",                 "NCh_FINAL_MUON");
  tree->SetAlias("mZd",                   "MPE_FINAL_MUON_Zd");
  tree->SetAlias("mAd",                   "MPE_FINAL_MUON_Ad");
  tree->SetAlias("mpbSigmaDeg",           "MPE_Pb_FINAL_MUON_SigmaDeg");

  TString mcTotalGeneratedEvents = "NEvents * " + TStringify(nFiles);

  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;



}



TTree* LoadTree_IC79_nugen_numu_6308_9981f_small() {

  TString file1 = LOADTREE_DIR+"/final.v19.nugen_numu_6308_9981f_RunID_new.root";

  TChain *tree = new TChain("L3CutsTree");

  tree->Add(file1);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  int nFiles =  9981;

  /*
  TString nufluxfile = LOADTREE_DIR+"/NuFlux_final_up.v4.nugen_numu_6308_9982f_v2.root";
  TChain *nufluxtree = new TChain("NuFlux");
  nufluxtree->Add(nufluxfile);

  assert(tree->GetNtrees() == nufluxtree->GetNtrees());
  tree->AddFriend(nufluxtree);
  */


  tree->SetAlias("mcPrimary_Zenith_rad",  "MC_Zr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "MC_Ar");
  tree->SetAlias("mcPrimary_Energy_GeV",  "MC_En");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "MuE_FINAL_MUON_En");
  tree->SetAlias("muexEn",                "MuEX_FINAL_MUON_En");
  tree->SetAlias("NChan",                 "NCh_FINAL_MUON");
  tree->SetAlias("mZd",                   "MPE_FINAL_MUON_Zd");
  tree->SetAlias("mAd",                   "MPE_FINAL_MUON_Ad");
  tree->SetAlias("mpbSigmaDeg",           "MPE_Pb_FINAL_MUON_SigmaDeg");


  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  //tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  //tree->SetAlias("BartolFluxWeight",      "NuFlux.bartol_naumov");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  TString mcTotalGeneratedEvents = "NEvents * " + TStringify(nFiles);



  tree->SetAlias("mcTotalGeneratedEvents", mcTotalGeneratedEvents);
  //tree->SetAlias("NFiles",TStringify(nFiles));

  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}




TTree *LoadTree_IC79_GoodRuns_Sirin_small() {

  TString files = "/net/user/sodrowski/IC79NeutrinoLevelData/exp/L6_up1.1_downv4.2_data_v13_runs_x*_grl2.5.root";
  TChain *tree = new TChain("CutTree");

  tree->Add(files);
  tree->GetEntries();

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");

  double livetimeTotal = 27315619.1; // 316 days
  TString livetimeTotalStr = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(livetimeTotal);
  tree->SetAlias("livetimeTotal", livetimeTotalStr);

  tree->SetAlias("mueEn",                 "MuE_FINAL_MUON_En");
  tree->SetAlias("muexEn",                "MuEX_FINAL_MUON_En");
  tree->SetAlias("NChan",                 "NCh_FINAL_MUON");
  tree->SetAlias("mZd",                   "MPE_FINAL_MUON_Zd");
  tree->SetAlias("mAd",                   "MPE_FINAL_MUON_Ad");
  tree->SetAlias("mpbSigmaDeg",           "MPE_Pb_FINAL_MUON_SigmaDeg");

  TString tmin = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(55347.3); //in MJD
  TString tmax = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(55694.1); //in MJD

  tree->SetAlias("tmin", tmin);
  tree->SetAlias("tmax", tmax);

  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;


}

TTree* LoadTree_IC79_GoodRuns_BDT_33d_small() {

  TString files = LOADTREE_DIR+"/final_up.v6.GoodRunList_33d.root";

  TChain *tree = new TChain("L3CutsTree");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  double livetimeTotal = 2858652.0; //seconds
  TString livetimeTotalStr = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(livetimeTotal);
  double tmin = 54971.154549; //in MJD
  double tmax = 55347.284392;

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree



  tree->SetAlias("livetimeTotal", livetimeTotalStr);
  tree->SetAlias("tmin", TStringify(tmin));
  tree->SetAlias("tmax", TStringify(tmax));
  //--------------------//

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");


  tree->SetAlias("mueEn",                 "MuE_FINAL_MUON_En");
  tree->SetAlias("muexEn",                "MuEX_FINAL_MUON_En");
  tree->SetAlias("NChan",                 "NCh_FINAL_MUON");
  tree->SetAlias("mZd",                   "MPE_FINAL_MUON_Zd");
  tree->SetAlias("mAd",                   "MPE_FINAL_MUON_Ad");
  tree->SetAlias("mpbSigmaDeg",           "MPE_Pb_FINAL_MUON_SigmaDeg");


  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;
}

TTree* LoadTree_IC79_GoodRuns_BDT_316d_small() {

  TString files = LOADTREE_DIR+"/final.v19.GoodRunList_316d_RunID_new.root";

  TChain *tree = new TChain("L3CutsTree");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  //double livetimeTotal = 28068353.1;//325 days seconds

  double livetimeTotal = 27315619.1;//316 days seconds

  TString livetimeTotalStr = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(livetimeTotal);
  TString tmin = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(55347.3); //in MJD
  TString tmax = " NCh_FINAL_MUON/NCh_FINAL_MUON * " + TStringify(55694.1); //in MJD

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree



  tree->SetAlias("livetimeTotal", livetimeTotalStr);
  tree->SetAlias("tmin", tmin);
  tree->SetAlias("tmax", tmax);
  //--------------------//

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");


  tree->SetAlias("mueEn",                 "MuE_FINAL_MUON_En");
  tree->SetAlias("muexEn",                "MuEX_FINAL_MUON_En");
  tree->SetAlias("NChan",                 "NCh_FINAL_MUON");
  tree->SetAlias("mZd",                   "MPE_FINAL_MUON_Zd");
  tree->SetAlias("mAd",                   "MPE_FINAL_MUON_Ad");
  tree->SetAlias("mpbSigmaDeg",           "MPE_Pb_FINAL_MUON_SigmaDeg");


  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;
}


