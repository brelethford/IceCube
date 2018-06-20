/*
 * compile using:
 * g++ -Wall EffAreaIC86-II-III.cpp `root-config --cflags --glibs` -o EffAreaIC86-II-III
*/
#include <iostream>
#include <sstream> 

#include <TString.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCut.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

TH1D* FastEffectiveArea(TChain *tr, TCut cut="", double zenMinDeg=90, double zenMaxDeg=180, int EBins=40, double LogEMin=1, double LogEMax=9,TString hname="hFastEffArea") {

    TH1D *hFastEffArea = new TH1D(hname,"Effective Area",EBins, LogEMin, LogEMax);

    cout << zenMinDeg << " " << zenMaxDeg << endl;

    double zenMinRad = zenMinDeg*TMath::DegToRad();
    double zenMaxRad = zenMaxDeg*TMath::DegToRad();

    double solidAngleFactor = 1./(2 * TMath::Pi() * (TMath::Cos(zenMinRad)-TMath::Cos(zenMaxRad)));

    //I3MCWeightDict.PowerLawIndex
    
    stringstream EBinsPerDecade;
    EBinsPerDecade << EBins/(LogEMax-LogEMin);
    tr->SetAlias("mcLogEBin", "int(log10(MCPrimary1.energy)*"+TString(EBinsPerDecade.str())+")");
    tr->SetAlias("mcEMin", "pow(10., mcLogEBin/"+TString(EBinsPerDecade.str())+")");
    tr->SetAlias("mcEMax", "pow(10., (1+mcLogEBin)/"+TString(EBinsPerDecade.str())+")");		 

    stringstream s;
    s << "MCPrimary1.zenith>" << zenMinRad;
    TCut angleMinCut = TCut(TString(s.str()));
    s.str("");
    s.clear();
    s << "MCPrimary1.zenith<" << zenMaxRad;    
    TCut angleMaxCut = TCut(TString(s.str()));

    s.str("");
    s.clear();    
    s << solidAngleFactor*1e-4;
    
    tr->Draw("log10(MCPrimary1.energy) >> "+ TString(hFastEffArea->GetName()),angleMinCut*angleMaxCut*TCut(TString(s.str())+"*I3MCWeightDict.OneWeight/(NFiles.NFiles*I3MCWeightDict.NEvents*(mcEMax-mcEMin))"));

    return hFastEffArea;
}

int main(){
    TString LOADTREE_ANADIR="/data/IceCube/IC862moreyeartars/";
    vector<TString> in_file;
    
    in_file.push_back("IC86_2012_NuGen_11069_Downgoing_6089files.root");
    in_file.push_back("IC86_2012_NuGen_11069_Upgoing_6089files.root");
    in_file.push_back("IC86_2012_NuGen_11070_Downgoing_4829files.root");
    in_file.push_back("IC86_2012_NuGen_11070_Upgoing_4829files.root");
    
    TChain *tr  =new TChain("MasterTree");
    TChain *trF =new TChain("NFiles");

    for (unsigned int i=0;i<in_file.size();i++){
        tr->Add(LOADTREE_ANADIR+in_file[i]);
        trF->Add(LOADTREE_ANADIR+in_file[i].ReplaceAll(".root","_NfilesTree.root"));
    }
    tr->AddFriend(trF);
    
    tr->GetEntries();

    vector <double*> Zdeg;
    double Ztmp0[3]={0, 60,2};
    Zdeg.push_back(Ztmp0);
    double Ztmp1[3]={60, 85,42};
    Zdeg.push_back(Ztmp1);
    double Ztmp2[3]={85, 120,8};
    Zdeg.push_back(Ztmp2);
    double Ztmp3[3]={120, 180,4};
    Zdeg.push_back(Ztmp3);
    
    vector <TH1D*> hFastEffArea;
    
    TCanvas *can=new TCanvas();
    
    TLegend *leg=new TLegend(0.1,0.8,0.3,1.);
    stringstream s;
    for (unsigned int i=0; i < Zdeg.size(); i++){
        s.str("");
        s.clear();    
        s << "hFastEffArea_" << i;
        cout << s.str() << endl;
        hFastEffArea.push_back(FastEffectiveArea(tr, "", Zdeg[i][0], Zdeg[i][1], 40, 1, 9, s.str()));
        hFastEffArea[i]->SetLineColor(Zdeg[i][2]);
        hFastEffArea[i]->SetStats(0);
    }
    for (unsigned int i=0; i < Zdeg.size(); i++){
        if (i==0) {
            hFastEffArea[i]->Draw();
            hFastEffArea[i]->GetXaxis()->SetTitle("log_{10}( E_{#nu} [GeV]) ");
            hFastEffArea[i]->GetYaxis()->SetTitle("Effective Area [m^{2}]");
        }
        else {
            hFastEffArea[i]->Draw("same");
        }
        can->Update();
        can->Modified();
        s.str("");
        s.clear();    
        s << Zdeg[i][0] -90. << " to " <<  Zdeg[i][1]-90.;
        leg->AddEntry(hFastEffArea[i],TString(s.str()), "l");
    }
    
    can->SetLogy();
    can->SetGridx();
    can->SetGridy();
    leg->Draw();
    can->SaveAs("EffAr.root");
    return 0;
}