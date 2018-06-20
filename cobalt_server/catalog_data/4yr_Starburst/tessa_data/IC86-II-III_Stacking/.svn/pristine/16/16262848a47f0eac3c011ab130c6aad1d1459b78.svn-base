{
    gROOT->Macro("$LAB_MAIN_DIR/llh/loadlibs.C");
    gInterpreter->AddIncludePath( gSystem->Getenv("LAB_MAIN_DIR") );

    if (!LOADSUCCESS) { return 0.; } // stop, and signal that script failed.
        
    initialize_ran1(-55);  // seed has to be a *NEGATIVE* integer

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/ArkExt.C+");

    bool OPT_USEREALDATA = true;
        
    //gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/mpfSigmaDegRescaledIC59.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC59.C+");
    loadSplines_IC59();
    //gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/mpfSigmaDegRescaledIC79Sirin.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC79.C+");
    loadSplines_IC79();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC86_SplineMPE.C+");
    loadSplines_IC86_SplineMPE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
    loadSplines_IC86_II_III_IV_SplineMPE();

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC59_Final.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC79_Final.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC86.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC86_II_III_IV.C");
    
    vector<I3Ark*> ark;
    ark.clear();
    TString RecoName = "SplineMPE";
    I3Ark *ark59=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_IC59.C(*ark59, OPT_USEREALDATA)");
    ark.push_back(ark59);
    
    
    I3Ark *ark79=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic79_Sirin.C(*ark79, OPT_USEREALDATA)");
    ark.push_back(ark79);
    
    I3Ark *ark86=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic86_BDT.C(*ark86, OPT_USEREALDATA, RecoName)");
    ark.push_back(ark86);
    
    I3Ark *ark862more=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic86_II_III_IV.C(*ark862more, OPT_USEREALDATA, RecoName)");
    ark.push_back(ark862more);
    

    MultiArk *mark=new MultiArk();
    for(int arki=0;arki<ark.size();arki++){
        mark->AddArk(*ark[arki]);
    }

    MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark->psData);

    int monLev = 0;
        
    // Specific to Milagro6:
    printf("Setting up Stacking Analysis Sources (6 Milagro SNR Sources)...\n");
    printf("Spectra are fit on the fly from Halzen, Kappes, O'Murchada (2008)\n");
        
        
    bool OPT_LOADSIGNAL = true;
    

    vector<I3Event> sourceEvents;
    sourceEvents.clear();
        
    // Set source sizes, used both for signal generation and spatial pdf
    vector<double> srcSigmas;
    srcSigmas.clear();
    for (int i=0; i<6; i++) { // Treat most as point sources
        srcSigmas.push_back(0.0); // Exceptions follow...
    }
    // Exceptions: ("Extent Diameter" taken as sigma to 2D circular Gaussian)
    srcSigmas[0] = 0.64; // MGRO J2019+37
    srcSigmas[1] = 1.3; //MGRO 1908+06 
    srcSigmas[2] = 1.5; // MGRO J2031+41
    srcSigmas[3] = 1.0; // MGRO J2043+36

    vector<EquatorialDeg> srcLocations;
    srcLocations.clear();
    vector<char*> srcNames;
    srcNames.clear();
    vector<double> weights;
    weights.clear();
    
    

    // Each Ark has its own event loader.  Each ark needs its own msg...
    //  (I.e. its own list of candidate source events)

    vector<EventLoaderExt*> evLoader;
    evLoader.clear();

    vector<I3MultiSignalGenerator*> msg;
    msg.clear();

    for(int arki=0;arki<ark.size();arki++){ 
        evLoader.push_back(&(ark[arki]->evLoader));
        evLoader[arki]->SetNUpSamples(20); // Default number of times each signal event is moved
        I3MultiSignalGenerator *msg_tmp=new I3MultiSignalGenerator();
        msg.push_back(msg_tmp);
    }
    // Source 0   r.a.,  dec.
    //EquatorialDeg srcLocation14(305.22,36.83);  // Fermi position
    //EquatorialDeg srcLocation0(305.03,36.72);  // MGRO J2019+37, Original
    //EquatorialDeg srcLocation0(304.68,36.70);  // MGRO J2019+37, Updated
    srcLocations.push_back(EquatorialDeg(304.85,36.80));
    srcNames.push_back("MGRO J2019+37");
    // Source 1
    //EquatorialDeg srcLocation1(286.89,6.03);  // Fermi position
    //EquatorialDeg srcLocation1(287.27,6.18);  // MGRO J1908+06, Original
    srcLocations.push_back(EquatorialDeg(286.68,5.83));
    srcNames.push_back("MGRO J1908+06");
    // Source 2
    //???? MGRO J2034+41 on Jim's pageark59
    //EquatorialDeg srcLocation2(308.06,41.38);  // TeV 2032+41 (assoc with J2031+42)
    //EquatorialDeg srcLocation2(308.04,41.57);  // MGRO J2031+41, Original
    srcLocations.push_back(EquatorialDeg(307.93,40.67));
    srcNames.push_back("MGRO J2031+42");
    // Source 3
    srcLocations.push_back(EquatorialDeg(310.98,36.3));
    srcNames.push_back("MGRO J2043+36");
    // Source 4
    srcLocations.push_back(EquatorialDeg(307.75,36.52));
    srcNames.push_back("MGRO J2032+37");
    // Source 5, added to v2 of Francis's paper, in the text it states that
    //  it has a flux about 2.5 times that of J2019+37 (next-brightest)
    srcLocations.push_back(EquatorialDeg(283.12,0.51));
    srcNames.push_back("MGRO J1852+01");
        
    double params[3];
    vector<double> enhancementFactors;
    enhancementFactors.clear();
    gROOT->ProcessLine(".L FitMilagroSpectra.C");

    for(int srcIndex=0;srcIndex<srcNames.size();srcIndex++){ 
        // For the theoretical spectra, the differences in flux are built
        //  directly into the flux fits of Halzen, Kappes, O'Murchada's
        //  simulated spectral data directly.
        // Just give 1.0 for each source here.
        enhancementFactors.push_back(1.);

        FitMilagroSpectra(srcNames[srcIndex],params);
        weights.push_back(pow(10,params[0]));
        if (OPT_LOADSIGNAL) {
            TF1 *fCutoff=new TF1(TString("fCutoff")+TString(srcNames[srcIndex]),"pow(10,[0]+[1]*log10(x)-pow(10,log10(x))/pow(10,[2])*log10(exp(1)))",1.e3,1.e6);
            fCutoff->SetParameters(params);
            FormulaFlux formflux(*fCutoff);
            FluxBase *myFluxPtr = &formflux;
            for(int arki=0;arki<ark.size();arki++){ 
                evLoader[arki]->LoadSourcePdf_Kent(srcLocations[srcIndex], srcSigmas[srcIndex]);
                evLoader[arki]->LoadSourceEvents(sourceEvents);
                I3PointGenerator i3point(sourceEvents, *myFluxPtr, srcLocations[srcIndex], ark[arki]->livetime);
                msg[arki]->AddSignal(i3point,enhancementFactors[srcIndex]);
            }
        }
    }
    // ... finished adding sources
    vector<SourceModule*> srcMods;
    
    for(int arki=0;arki<ark.size();arki++){ 
        srcMods.push_back(msg[arki]);
    }
    // CALCULATE STACKING WEIGHTS VS GAMMA TABLE
    // Set theoretical bias, i.e. if you have reason to believe one source emits
    //  a higher flux, add that enhancement factor here and the stacking search
    //  will use that enhancement factor to construct the signal spatial pdfs
    //  e.g., if you believe one source emits twice as much flux as another, the
    //  srcWeightsArray will scale that source pdf up by a factor of 2 if set here.
    //  (It will not inject more flux, such as setting the enhancement factor
    //   of the I3PointGenerator)
        
    vector<double> pdfEnhancement;
    pdfEnhancement.clear();
    
    for (int i=0; i<6; i++) {
        // For equal weights
        pdfEnhancement.push_back(1./6.); // Identical, i.e. no src enhancments
        // Use norms from fits to theoretical spectra
        //double weightSum = weights[0]+weights[1]+weights[2]+weights[3]+weights[4]+weights[5];
        //pdfEnhancement.push_back(weights[i]/weightSum);
        //cout << pdfEnhancement[i] << " " << weights[i] << endl;
    }
    
    // At least make sure these are the same length as number of sources
    assert ( srcSigmas.size() == srcLocations.size() ); // else misconfigured
    assert ( enhancementFactors.size() == srcLocations.size() ); // else misconfigured
    assert ( pdfEnhancement.size() == srcLocations.size() ); // else misconfigured
    
    cout << "\nCalculating Stacking Weight vs Gamma Tables...\n ";
    
    double gammaMin = -4;
    double gammaMax = -1;
    int nGammaBins  = 30;

    vector<vector<vector<double> > > srcWeightsArray;

    for (int arki=0; arki<ark.size(); arki++) {     
        vector<vector<double> > srcWeightsArray_tmp;
        int vSize = srcLocations.size();
        srcWeightsArray_tmp.resize( vSize );
        double srcWeight=0, testGamma=0;
        for (int srcIndex=0; srcIndex < srcLocations.size(); srcIndex++) {   
            evLoader[arki]->LoadSourceEvents(sourceEvents, srcLocations[srcIndex]);
            for ( int gammaIndex=0; gammaIndex<nGammaBins; gammaIndex++) {
                double testGamma = (gammaIndex+0.5)*(gammaMax-gammaMin)/nGammaBins + gammaMin;
                PowerLawFlux testFlux(1,testGamma);  //  1 GeV^-1 cm^-2 s^-1;  index = testGamma
                FluxBase *myFluxPtr = &testFlux;
                I3PointGenerator i3point(sourceEvents, *myFluxPtr, srcLocations[srcIndex], ark[arki]->livetime);
                srcWeight = i3point->GetMeanSrcNev();
                srcWeightsArray_tmp[srcIndex].push_back(srcWeight*pdfEnhancement[srcIndex]);
            }
            
        }
    
        // Normalize each row for sanity's sake (stacking method does not require this)
        cout << endl;
        cout << "Normalized table of Src Weights for different gamma:" << endl;
        cout << "x-axis = Source Index, y-axis = gamma index (soft to hard)"<<endl;
        cout << endl;
        for ( int gammaIndex=0; gammaIndex<nGammaBins; gammaIndex++) {
            printf("%02d:  ",gammaIndex);
            double sumRow = 0;
            // Find norm constant
            for (int srcIndex=0; srcIndex<vSize; srcIndex++){
                sumRow += srcWeightsArray_tmp[srcIndex][gammaIndex];
            }
            // And normalize
            for (int srcIndex=0; srcIndex<vSize; srcIndex++){
                srcWeightsArray_tmp[srcIndex][gammaIndex]/=sumRow;
                printf("%.4f ",srcWeightsArray_tmp[srcIndex][gammaIndex]);
            }
            cout << endl;
        }
        srcWeightsArray.push_back(srcWeightsArray_tmp);
    }
    
    vector<NewLlhStack*> llhs;
    MultiAnalysisFn maf;
    for (int arki=0; arki<ark.size(); arki++) { 
        NewLlhStack *tmp_llh=new NewLlhStack();
        tmp_llh->SetUseEnergy(true);
        tmp_llh->SetOptimizeTolerance(0.00);
        tmp_llh->SetMonitorLevel(monLev);
        tmp_llh->SetAnalysisSet(ark[arki]->psData);
        tmp_llh->SetSourceCoords(srcLocations);
        tmp_llh->SetStackedWeightTable(srcWeightsArray[arki]);
        tmp_llh->SetSourceSigmas(srcSigmas);
        llhs.push_back(tmp_llh);
        maf.AddAnalysisFn(tmp_llh);
        ark[arki]->SetSource(srcMods[arki]);//msg[arki]
    }

    NewLlhEnergy_ParTranslator pt;
    pt.SetRange(1,4,31);
    pt.SetTranslator(mas);
    maf.SetParTranslator(&pt);

    vector<MinuitParDef> pdv;

    pdv.push_back( MinuitParDef("nSrc",2,0.1, 0.,100.) );
    
    pdv.push_back( MinuitParDef("gamma",2.5,0.5, 1., 4.) );
    maf.SetParDefs(pdv);

    maf.MaximizeLlh();

    cout<<"GetEstProb:  "<<maf.GetEstProb()<<endl;
    cout<<"GetTestStat: "<<maf.GetTestStatistic()<<endl;
    cout<<"Get Ns:      "<<maf.GetPar(0)<<endl;
    cout<<"Get gamma:   "<<maf.GetPar(1)<<endl;

     gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/macro_genHotSpot_new.C");

}
