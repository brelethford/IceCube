{
    gROOT->Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C");
    initialize_ran1(-55);

    gROOT->ProcessLine(".L ArkTime.C+");
    bool OPT_USEREALDATA = false;

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC86_II_III_IV_lc.C");
    I3Ark arkIC86II_III_IV;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic86_II_III_IV_tdep_lc.C(arkIC86II_III_IV, OPT_USEREALDATA,\"SplineMPE\")");

    int monLev = 0;

    NewLlhGausTime newllh86II_to_IV;
    newllh86II_to_IV.SetUseEnergy(true);
    //newllh86II_to_IV.SetOptimizeAngleDeg(10.);
    newllh86II_to_IV.SetOptimizeTolerance(0.01);
    newllh86II_to_IV.SetMonitorLevel(0);
    newllh86II_to_IV.SetEMaxRatioWarnOnlyOnce(1);
    newllh86II_to_IV.close_ = 10.;
    newllh86II_to_IV.JimsTerm_ = true;
    newllh86II_to_IV.SpectralPenalty_ = false;
    //newllh86II_to_IV.ndof_ = 3.;
    newllh86II_to_IV.ndof_ = 4.;
    newllh86II_to_IV.SetLivetime(arkIC86II_III_IV.livetime/86400.);
    newllh86II_to_IV.SetLocalCoordBkgProb(arkIC86II_III_IV.lcBkgProb);


    EquatorialDeg testSearch(343.491,16.148);
    double spectralIndex = -2;
    double tmean = (arkIC86II_III_IV.tmax + arkIC86II_III_IV.tmin)/2.;

    TimePdf * tPdf = new GaussianTimePdf(arkIC86II_III_IV.tmin, arkIC86II_III_IV.tmax, tmean, 1e-5, 1.);
    arkIC86II_III_IV.SetPointSource(testSearch, PowerLawFlux(1.,spectralIndex), tPdf);
    arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( arkIC86II_III_IV.lcBkgProb.nbAz );
    newllh86II_to_IV.SetTimeBounds(tPdf);

        // plot ene pdf
    int nbinsEn = 440;
    int nbinsZn = 2500;
    double val;
    arkIC86II_III_IV.decBkgProb.FixToBase();
    //TH2D * hePdf = new TH2D("","",nbinsEn,2,8,nbinsZn,-1,1); // This matches the space bkg pdf
    //TH2D * hePdf = new TH2D("","",nbinsZn,-1,1,nbinsEn,2,8);
    TH2D * hePdf = new TH2D("","",nbinsZn,-1,1,nbinsEn,2,8);
    EquatorialDeg coo;
    I3Event ev;
    I3EventParameters pa;
    for (int iEn=1;iEn<=nbinsEn;iEn++) {
        //pa.energyValue = hePdf->GetXaxis()->GetBinCenter(iEn);
        pa.energyValue = hePdf->GetYaxis()->GetBinCenter(iEn);
        for (int iDec=1;iDec<=nbinsZn;iDec++) {
            //pa.recoZenithDeg = hePdf->GetYaxis()->GetBinCenter(iDec);
            //pa.recoZenithDeg = 57.3*acos(hePdf->GetYaxis()->GetBinCenter(iDec));
            pa.recoZenithDeg = 57.3*acos(hePdf->GetXaxis()->GetBinCenter(iDec));
            //pa.recoZenithDeg = 90.+57.3*asin(hePdf->GetYaxis()->GetBinCenter(iDec));
            //cout << pa <<  " " << flush;
            //cout << "pa.recoZenithDeg " << pa.recoZenithDeg << " pa.energyValue " << pa.energyValue << endl;
            ev.SetParams(pa);
            val = arkIC86II_III_IV.eProb->GetEnergyProbBkg(ev);
            //hePdf->SetBinContent(iEn,iDec,val); // this matches local coord bkg in zenith
            hePdf->SetBinContent(iDec,iEn,val);
        }
    }  
    hePdf->GetYaxis()->SetTitle("log_{10} (Energy Proxy)");
    //hePdf->GetYaxis()->SetTitle("Zenith (#circ)");
    hePdf->GetXaxis()->SetTitle("cos(Zenith)");
    hePdf->GetZaxis()->SetTitle("B^{Energy}_{i}");
    hePdf->GetZaxis()->SetTitleOffset(1);
    TCanvas * c = new TCanvas("c","c");
    c->SetRightMargin(0.125);
    hePdf->GetZaxis()->SetTitleOffset(0.5);
    hePdf->Draw("colz");

    
    
    newllh86II_to_IV.SetAnalysis(arkIC86II_III_IV.psData, testSearch);
    arkIC86II_III_IV.psData->GenerateDataSet_with_nSrcEvents(0);
    
    
    
    
    /*
    TH1D hTestStatistic;
    hTestStatistic.SetBins(100,0,1);
    hTestStatistic.SetTitle(";2 ln #lambda;trials");
    hTestStatistic.SetBit(TH1::kCanRebin);
    for (int i=0;i<1000;i++) {
        if (i%100==0) cout << i << endl;
        arkIC86II_III_IV.psData->GenerateDataSet_with_nSrcEvents(0);
        newllh86II_to_IV.MaximizeLlh();
        hTestStatistic.Fill(newllh86II_to_IV.GetTestStatistic() * 2.);
    }
    hTestStatistic.Draw();
    */
    newllh86II_to_IV.MaximizeLlh();
    cout <<
    newllh86II_to_IV.GetPar(0) << " " <<
    newllh86II_to_IV.GetPar(1) << " " <<
    newllh86II_to_IV.GetPar(2) << " " <<
    newllh86II_to_IV.GetPar(3) << " " << endl;
  return 1; // signal correct finish of script
}
