Float_t ddeclination;
TString outname;
Float_t dRA;
gROOT->SetBatch();

void runIC86II_III_IV_AllSkyTimeDep_discoE2(Float_t declination, Float_t RA){
    gROOT->SetBatch();
    ddeclination=declination;
    dRA=RA;
    std::ostringstream ss;
    ss << ddeclination;
    std::ostringstream ss2;
    ss2 << dRA;
    outname="IC86II_III_IV_AllSkyTimeDep_discoE2_dec"+ss.str()+"_RA_"+ss2.str()+".root";
    cout << "dec is "<< ddeclination << endl;
    gROOT->ProcessLine(".x IC86II_III_IV_allSkyTimeDep_disco.C"); //everything is done and saved in here
}
