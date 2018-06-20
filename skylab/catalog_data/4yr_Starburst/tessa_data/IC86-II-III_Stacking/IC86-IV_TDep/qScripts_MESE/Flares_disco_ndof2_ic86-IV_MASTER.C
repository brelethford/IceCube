{
    gROOT->SetBatch();
    
    gROOT->Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C");

    string blocksFile =  "MASTER_BLOCKFILE_MASTER";
    double ra = MASTER_RA_MASTER;
    double declination = MASTER_DECLINATION_MASTER;
    TString outfile="MASTER_OUT_FILE";
    TString siteatt="MASTER_SITEATT";

    cout << "Opening blockfile " << blocksFile.c_str() << endl;
    cout << "With coordinates " << ra << " " << declination << endl;

    initialize_ran1(-156); // seed has to be a *NEGATIVE* integer

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/ArkTime.C+");

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SimpleMultiBlockAnalysis.C");
    
    bool OPT_USEREALDATA = false;

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+");
    loadSplines_IC79_86_I_to_IV_SplineMPE_MESE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE"+siteatt+".C");

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
    loadSplines_IC86_II_III_IV_SplineMPE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC86_II_III_IV"+siteatt+".C");

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_SplineMPE.C+");
    loadSplines_IC86_SplineMPE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC86"+siteatt+".C");

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79.C+");
    loadSplines_IC79();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_Final"+siteatt+".C");

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC59.C+");
    loadSplines_IC59();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC59_Final"+siteatt+".C");

    I3Ark arkIC79_IV_MESE;
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC79_to_86IV_MESE_tdep_at.C(arkIC79_IV_MESE, OPT_USEREALDATA,\"SplineMPE\")");

    I3Ark arkIC86II_III_IV;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic86_II_III_IV_tdep"+siteatt+".C(arkIC86II_III_IV, OPT_USEREALDATA,\"SplineMPE\")");
    I3Ark arkIC86I;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic86_BDT_tdep"+siteatt+".C(arkIC86I, OPT_USEREALDATA,\"SplineMPE\")");
    I3Ark arkIC79;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic79_Final_tdep"+siteatt+".C(arkIC79, OPT_USEREALDATA)");
    I3Ark arkIC59;
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic59_BDT_tdep"+siteatt+".C(arkIC59, OPT_USEREALDATA)");

    MultiArk mark;
    mark.AddArk(arkIC79_IV_MESE);
    mark.AddArk(arkIC86II_III_IV);
    mark.AddArk(arkIC86I);
    mark.AddArk(arkIC79);
    mark.AddArk(arkIC59);

    int nSrcEvents = 0;
    EquatorialDeg testSearch(ra, declination);

    NewLlhBlockTime newllh79_IV_MESE;
    NewLlhBlockTime newllh86II_III_IV;
    NewLlhBlockTime newllh86I;
    NewLlhBlockTime newllh79;
    NewLlhBlockTime newllh59;

    MultiBlockAnalysisFn maf;
    NewLlhBlockTime_ParTranslator pt;

    bool useE = true;
    int monLevel = 0;
    double laglimit = 0.5;

    TimePdf * tPdf;

    newllh79_IV_MESE.SetUseEnergy(true);
    newllh79_IV_MESE.SetOptimizeTolerance(0.01);
    newllh79_IV_MESE.SetMonitorLevel(monLevel);
    newllh79_IV_MESE.SetEMaxRatioWarnOnlyOnce(1);
    newllh79_IV_MESE.JimsTerm_ = false;
    newllh79_IV_MESE.SpectralPenalty_ = false;
    newllh79_IV_MESE.ndof_ = 2.;
    newllh79_IV_MESE.laglimit_ = laglimit;
    newllh79_IV_MESE.SetLivetime(arkIC79_IV_MESE.livetime/86400.);
    newllh79_IV_MESE.SetBlocks(blocksFile,0.);
    newllh79_IV_MESE.SetAnalysisSet(arkIC79_IV_MESE.psData);

    newllh86II_III_IV.SetUseEnergy(true);
    newllh86II_III_IV.SetOptimizeTolerance(0.01);
    newllh86II_III_IV.SetMonitorLevel(monLevel);
    newllh86II_III_IV.SetEMaxRatioWarnOnlyOnce(1);
    newllh86II_III_IV.JimsTerm_ = false;
    newllh86II_III_IV.SpectralPenalty_ = false;
    newllh86II_III_IV.ndof_ = 2.;
    newllh86II_III_IV.laglimit_ = laglimit;
    newllh86II_III_IV.SetLivetime(arkIC86II_III_IV.livetime/86400.);
    newllh86II_III_IV.SetBlocks(blocksFile,0.);
    newllh86II_III_IV.SetAnalysisSet(arkIC86II_III_IV.psData);
    
    newllh86I.SetUseEnergy(true);
    newllh86I.SetOptimizeTolerance(0.01);
    newllh86I.SetMonitorLevel(monLevel);
    newllh86I.SetEMaxRatioWarnOnlyOnce(1);
    newllh86I.JimsTerm_ = false;
    newllh86I.SpectralPenalty_ = false;
    newllh86I.ndof_ = 2.;
    newllh86I.laglimit_ = laglimit;
    newllh86I.SetLivetime(arkIC86I.livetime/86400.);
    newllh86I.SetBlocks(blocksFile,0.);
    newllh86I.SetAnalysisSet(arkIC86I.psData);

    newllh79.SetUseEnergy(true);
    newllh79.SetOptimizeTolerance(0.01);
    newllh79.SetMonitorLevel(monLevel);
    newllh79.SetEMaxRatioWarnOnlyOnce(1);
    newllh79.JimsTerm_ = false;
    newllh79.SpectralPenalty_ = false;
    newllh79.ndof_ = 2.;
    newllh79.laglimit_ = laglimit;
    newllh79.SetLivetime(arkIC79.livetime/86400.);
    newllh79.SetBlocks(blocksFile,0.);
    newllh79.SetAnalysisSet(arkIC79.psData);

    newllh59.SetUseEnergy(true);
    newllh59.SetOptimizeTolerance(0.01);
    newllh59.SetMonitorLevel(monLevel);
    newllh59.SetEMaxRatioWarnOnlyOnce(1);
    newllh59.JimsTerm_ = false;
    newllh59.SpectralPenalty_ = false;
    newllh59.ndof_ = 2.;
    newllh59.laglimit_ = laglimit;
    newllh59.SetLivetime(arkIC59.livetime/86400.);
    newllh59.SetBlocks(blocksFile,0.);
    newllh59.SetAnalysisSet(arkIC59.psData);

    newllh79_IV_MESE.SetTimeBounds(arkIC79_IV_MESE.tmin, arkIC79_IV_MESE.tmax);
    cout << "newllh79_IV_MESE.SetTimeBounds " << arkIC79_IV_MESE.tmin << " "<<  arkIC79_IV_MESE.tmax << endl;
    newllh86II_III_IV.SetTimeBounds(arkIC86II_III_IV.tmin, arkIC86II_III_IV.tmax);
    cout << "newllh86II_III_IV.SetTimeBounds " << arkIC86II_III_IV.tmin << " "<<  arkIC86II_III_IV.tmax << endl;
    newllh86I.SetTimeBounds(arkIC86I.tmin, arkIC86I.tmax);
    cout << "newllh86I.SetTimeBounds " << arkIC86I.tmin << " "<<  arkIC86I.tmax << endl;
    newllh79.SetTimeBounds(arkIC79.tmin, arkIC79.tmax);
    cout << "newllh79.SetTimeBounds " << arkIC79.tmin << " "<<  arkIC79.tmax << endl;
    newllh59.SetTimeBounds(arkIC59.tmin, arkIC59.tmax);
    cout << "newllh59.SetTimeBounds " << arkIC59.tmin << " "<<  arkIC59.tmax << endl;

    maf.AddAnalysisFn(&newllh79_IV_MESE);
    maf.AddAnalysisFn(&newllh86II_III_IV);
    maf.AddAnalysisFn(&newllh86I);
    maf.AddAnalysisFn(&newllh79);
    maf.AddAnalysisFn(&newllh59);

    maf.blocksFile = blocksFile;

    maf.SetSearchCoord( testSearch );

    tPdf = new BlockTimePdf1();
    tPdf->SetBlockLevels(blocksFile,0.);

    mark.SetPointSource( testSearch, PowerLawFlux(1.,-2.), tPdf);
    pt.SetRange(1,4,31);
    pt.SetUpBlocks(blocksFile);
    MultiAnalysisSet* mas=dynamic_cast<MultiAnalysisSet*>(mark.psData);
    pt.SetTranslator(mas);
    maf.SetParTranslator(&pt);

    mark.psData->GenerateDataSet_with_nSrcEvents(0);
    
    DiscoveryPotential disco;

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SetDisco.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/MultiBlockDetectionStudy.C");
    // parameters:  loops, optMedianUpperLimit,  significance,  power);
    SetDisco(disco, 40, false, 2.87e-7, 0.5);
    double thMax = ( GetHighestBlock(blocksFile) + GetSecondHighestBlock(blocksFile) ) /2.;
    cout << "Max thershold " << thMax << endl;

    double n_disco[9], thTest[9], thPlot[9], ndays[9];

    for (int i = 0 ; i < 9; i++) {
        thTest[i] = thMax * i / 8.;
        thPlot[i] = thTest[i] * 1e6;
        ndays[i] = BlockTimeAboveThresh(blocksFile, thTest[i],arkIC59.tmin,arkIC86I.tmax);
        cout << "Testing Threshold of: " << thTest[i] << endl;
        TimePdf * tPdf2 = new BlockTimePdf1();
        tPdf2->SetBlockLevels(blocksFile,thTest[i]);
        mark.SetPointSource( testSearch, PowerLawFlux(1.,-2.), tPdf2);
        n_disco[i] = MultiBlockDetectionStudy_d(mark, maf, mas, disco);
    }
    TCanvas * canComp = new TCanvas("canComp","canComp",1200,500);
    canComp->Divide(2,1);
    
    canComp->cd(1);
    
    TGraph * g = new TGraph(9,thTest,n_disco);
    g->SetName("gDisco");
    g->GetXaxis()->SetTitle("Threshold for Emission (Photon cm^{-2} s^{-1})");
    g->GetYaxis()->SetTitle("N Events for Discovery in 50% of trials");
    g->SetLineWidth(3);
    g->SetLineColor(4);
    g->Draw("AL");
    leg = new TLegend(0.25,0.8,0.95,0.99);
    leg->AddEntry(g,"P=0.5 5#sigma discovery, Fit Threshold","l");
    leg->Draw();
    
    canComp->cd(2);

    TGraph * gd = new TGraph(9,ndays,n_disco);
    gd->SetName("gDiscoDays");
    gd->GetXaxis()->SetTitle("Time above threshold for emission (days)");
    
    gd->GetYaxis()->SetTitle("N Events for Discovery in 50% of trials");
    gd->SetLineWidth(3);
    gd->SetLineColor(4);
    gd->Draw("AL");
    canComp->SetLogx();
    cout << "saving file "<< outfile << endl;
    canComp->SaveAs(outfile+".png");
    canComp->SaveAs(outfile+".root");

}
