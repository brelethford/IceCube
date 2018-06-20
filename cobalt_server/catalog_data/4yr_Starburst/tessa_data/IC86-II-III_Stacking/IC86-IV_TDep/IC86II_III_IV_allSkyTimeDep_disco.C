{

    gROOT->Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C");
    gROOT->SetBatch();
    cout << "disco for dec " <<  ddeclination << " and RA "<< dRA << " will be stored in " << outname << endl;

    initialize_ran1(-55);

    bool OPT_USEREALDATA = false;

    gROOT->ProcessLine(".L ArkTime.C+");
    gROOT->ProcessLine(".L SetDisco.C");
    gROOT->ProcessLine(".L DetectionStudy.C");

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
    newllh86II_to_IV.ndof_ = 3.;
    newllh86II_to_IV.SetLivetime(arkIC86II_III_IV.livetime/86400.);
    newllh86II_to_IV.SetLocalCoordBkgProb(arkIC86II_III_IV.lcBkgProb);

    EquatorialDeg testSearch(dRA, ddeclination); //To macth Mike's previous discovery plots at delta = 16 deg
    cout << "EquatorialDeg dec: " << ddeclination << endl;
    double spectralIndex = -2;
    double tmean = (arkIC86II_III_IV.tmax + arkIC86II_III_IV.tmin)/2.;

    TimePdf * tPdf = new GaussianTimePdf(arkIC86II_III_IV.tmin, arkIC86II_III_IV.tmax, tmean, 1e-5, 1.);
    arkIC86II_III_IV.SetPointSource(testSearch, PowerLawFlux(1.,spectralIndex), tPdf);
    arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( arkIC86II_III_IV.lcBkgProb.nbAz );
    newllh86II_to_IV.SetTimeBounds(tPdf);
    newllh86II_to_IV.SetAnalysis(arkIC86II_III_IV.psData, testSearch);
  
    DiscoveryPotential disco;
  
    // parameters:  loops, optMedianUpperLimit,  significance,  power);
    SetDisco(disco, 50, false, 2.87e-7, 0.5);
    // SetDisco(disco, 15, true, 0.5, 0.9);
  
    cout << "setting sigmas" << endl;
    double n_disco[12], ldsigmas[12],fscale[12];
    double dsigmas[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 1., 3., 5., 10, 1e2};
    for (int i = 0; i < 12; i++) {
        TimePdf * tPdf2 = new GaussianTimePdf(arkIC86II_III_IV.tmin, arkIC86II_III_IV.tmax, tmean,dsigmas[i],1.);
        arkIC86II_III_IV.SetPointSource(testSearch, PowerLawFlux(1.,spectralIndex), tPdf2);
        
	if (dsigmas[i] < 8e-2) {
            arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( dsigmas[i] * 360.);
            cout << "short itme bin fix" << endl;
        }
	else arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( arkIC86II_III_IV.lcBkgProb.nbAz );
        
        newllh86II_to_IV.SetTimeBounds(tPdf2);
  
        n_disco[i] = DetectionStudy_d(arkIC86II_III_IV, newllh86II_to_IV, disco,fscale[i]);
        ldsigmas[i] = log10(dsigmas[i]);
        cout << "dsigmas " << dsigmas[i] << " n_disco " << n_disco[i] << endl;
    }
    cout << "sigmas loop done" << endl;
    
    TGraph * g = new TGraph(12,ldsigmas,n_disco);
    g->SetNameTitle("disco","disco");
    
    TGraph * g2 = new TGraph(12,ldsigmas,fscale);
    g2->SetNameTitle("fscale","fscale");

    SetDisco(disco, 50, true, 0.5, 0.9);
    cout << "setting sigmas" << endl;
    for (int i = 0; i < 12; i++) {
        TimePdf * tPdf2 = new GaussianTimePdf(arkIC86II_III_IV.tmin, arkIC86II_III_IV.tmax, tmean,dsigmas[i],1.);
        arkIC86II_III_IV.SetPointSource(testSearch, PowerLawFlux(1.,spectralIndex), tPdf2);
        if (dsigmas[i] < 8e-2) {
            arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( dsigmas[i] * 360.);
            cout << "short itme bin fix" << endl;
        }
        else arkIC86II_III_IV.psData.GetSource().SetTimeAzBins( arkIC86II_III_IV.lcBkgProb.nbAz );

        newllh86II_to_IV.SetTimeBounds(tPdf2);

        n_disco[i] = DetectionStudy_d(arkIC86II_III_IV, newllh86II_to_IV, disco,fscale[i]);
        ldsigmas[i] = log10(dsigmas[i]);
        cout << "dsigmas " << dsigmas[i] << " n_disco " << n_disco[i] << endl;
    }
    cout << "sigmas loop done" << endl;
    TGraph * s = new TGraph(12,ldsigmas,n_disco);
    s->SetNameTitle("sens","sens");

    TGraph * s2 = new TGraph(12,ldsigmas,fscale);
    s2->SetNameTitle("fscale","fscale");


    TFile *outf=new TFile(outname,"recreate");
    g->Write();
    g2->Write();
    s->Write();
    s2->Write();

    outf->Close();
    return 1; // signal correct finish of script
}
