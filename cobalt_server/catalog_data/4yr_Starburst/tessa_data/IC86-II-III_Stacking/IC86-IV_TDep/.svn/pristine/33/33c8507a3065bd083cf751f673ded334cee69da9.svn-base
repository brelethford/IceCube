{
    gROOT->Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C");
    initialize_ran1(-55);

    gROOT->ProcessLine(".L ArkTime.C+");
    bool OPT_USEREALDATA = false;

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE_lc.C");
    I3Ark arkIC79_IV_MESE;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC79_86_I_to_IV_MESE_tdep_lc.C(arkIC79_IV_MESE, OPT_USEREALDATA,\"SplineMPE\")");

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
    newllh86II_to_IV.SetLivetime(arkIC79_IV_MESE.livetime/86400.);
    newllh86II_to_IV.SetLocalCoordBkgProb(arkIC79_IV_MESE.lcBkgProb);


    EquatorialDeg testSearch(343.491,16.148);
    double spectralIndex = -2;
    double tmean = (arkIC79_IV_MESE.tmax + arkIC79_IV_MESE.tmin)/2.;

    TimePdf * tPdf = new GaussianTimePdf(arkIC79_IV_MESE.tmin, arkIC79_IV_MESE.tmax, tmean, 1e-5, 1.);
    arkIC79_IV_MESE.SetPointSource(testSearch, PowerLawFlux(1.,spectralIndex), tPdf);
    arkIC79_IV_MESE.psData.GetSource().SetTimeAzBins( arkIC79_IV_MESE.lcBkgProb.nbAz );
    newllh86II_to_IV.SetTimeBounds(tPdf);

    newllh86II_to_IV.SetAnalysis(arkIC79_IV_MESE.psData, testSearch);
    arkIC79_IV_MESE.psData->GenerateDataSet_with_nSrcEvents(0);
    
    TH1D hTestStatistic;
    hTestStatistic.SetBins(100,0,1);
    hTestStatistic.SetTitle(";2 ln #lambda;trials");
    hTestStatistic.SetBit(TH1::kCanRebin);
    for (int i=0;i<1000;i++) {
        if (i%100==0) cout << i << endl;
        arkIC79_IV_MESE.psData->GenerateDataSet_with_nSrcEvents(0);
        newllh86II_to_IV.MaximizeLlh();
        hTestStatistic.Fill(newllh86II_to_IV.GetTestStatistic() * 2.);
    }
    hTestStatistic.Draw();
    
    newllh86II_to_IV.MaximizeLlh();
    cout <<
    newllh86II_to_IV.GetPar(0) << " " <<
    newllh86II_to_IV.GetPar(1) << " " <<
    newllh86II_to_IV.GetPar(2) << " " <<
    newllh86II_to_IV.GetPar(3) << " " << endl;
  return 1; // signal correct finish of script
}
