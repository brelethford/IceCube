{
    gROOT->Macro("$LAB_MAIN_DIR/llh/loadlibs.C");
    gInterpreter->AddIncludePath( gSystem->Getenv("LAB_MAIN_DIR") );

    if (!LOADSUCCESS) { return 0.; } // stop, and signal that script failed.
        
    initialize_ran1(-55);  // seed has to be a *NEGATIVE* integer

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/ArkExt.C+");

    bool OPT_USEREALDATA = true;
        
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC40.C+");
    loadSplines_IC40();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC59.C+");
    loadSplines_IC59();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79.C+");
    loadSplines_IC79();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_SplineMPE.C+");
    loadSplines_IC86_SplineMPE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
    loadSplines_IC86_II_III_IV_SplineMPE();
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+");
    loadSplines_IC79_86_I_to_IV_SplineMPE_MESE();

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC40_CutA6_Fix_final.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC59_Final.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC79_Final.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC86.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC86_II_III_IV.C");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/TreeLoader_IC79_86_I_to_IV_MESE.C");    

    vector<I3Ark*> ark;
    ark.clear();

    I3Ark *ark40=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic40.C(*ark40, OPT_USEREALDATA)");
    ark.push_back(ark40);

    I3Ark *ark59=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_IC59.C(*ark59, OPT_USEREALDATA)");
    ark.push_back(ark59);
    
    
    I3Ark *ark79=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic79_Sirin.C(*ark79, OPT_USEREALDATA)");
    ark.push_back(ark79);
    
    I3Ark *ark86=new I3Ark();
    TString RecoName = "SplineMPE";
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic86_BDT.C(*ark86, OPT_USEREALDATA, RecoName)");
    ark.push_back(ark86);
    
    I3Ark *ark862more=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic86_II_III_IV.C(*ark862more, OPT_USEREALDATA, RecoName)");
    ark.push_back(ark862more);
    
    I3Ark *arkMESE=new I3Ark();
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_IC79_86_I_to_IV_MESE.C(*arkMESE, OPT_USEREALDATA, RecoName)");
    ark.push_back(arkMESE);

    MultiArk *mark=new MultiArk();
    for(int arki=0;arki<ark.size();arki++){
        mark->AddArk(*ark[arki]);
    }

    MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark->psData);

    int monLev = 0;
        
    printf("Setting up Stacking Analysis Sources ...\n");
        
        
    bool OPT_LOADSIGNAL = true;
    

    vector<I3Event> sourceEvents;
    sourceEvents.clear();
        
    // Set source sizes, used both for signal generation and spatial pdf
    vector<double> srcSigmas;
    srcSigmas.clear();
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

// -----------------Reading Files from Script -------------------------------------
    string single_line;

    TString s;

    ifstream cloud_data_count("HAWC_sources.csv");

    TObjArray *tx;

    float Dec=0, Ra=0, Sig=0, extension, TotalWeight=0.0;
    string detector, name;
    while (1){       

        if (!getline(cloud_data_count,single_line )){cout<< "Breaking now!" <<endl;break;};

        s=TString(single_line);

        cout << s << endl;

        tx=s.Tokenize(" ");
        detector=((TObjString *)(tx->At(0)))->String();
        name=((TObjString *)(tx->At(1)))->String();
        Dec=((TObjString *)(tx->At(2)))->String()->Atof();      
        Ra=((TObjString *)(tx->At(3)))->String()->Atof();      
        Sig=((TObjString *)(tx->At(4)))->String()->Atof();      
        extension=((TObjString *)(tx->At(5)))->String()->Atof(); 

        cout << detector << " " << name ;
        cout << " This is a float " << Dec << " RA: "<< Ra << " Sig: " << Sig << " ext: "<< extension << endl;
	srcNames.push_back(malloc(strlen(name.c_str())) + 1);
	strcpy(srcNames.back(),name.c_str());
	weights.push_back(1.0);//(Sig);
	srcSigmas.push_back(extension);
	TotalWeight+=1.0;
	srcLocations.push_back(EquatorialDeg(Ra,Dec));
    }

    for(int srcIndex=0;srcIndex<srcNames.size();srcIndex++){
        cout<<srcNames[srcIndex]<<" "<<weights[srcIndex]<<" "<< srcSigmas[srcIndex]<<" "<<srcLocations[srcIndex].GetRa() << " " << srcLocations[srcIndex].GetDec() <<" "<<TotalWeight<<endl;
    }

    double params[3];
    vector<double> enhancementFactors;
    enhancementFactors.clear();
    double spectralIndex=-2.0;

    for(int srcIndex=0;srcIndex<srcNames.size();srcIndex++){ 
        // For the theoretical spectra, the differences in flux are built
        //  directly into the flux fits of Halzen, Kappes, O'Murchada's
        //  simulated spectral data directly.
        // Just give 1.0 for each source here.
        enhancementFactors.push_back(weights[srcIndex]/TotalWeight);
	cout<<srcNames[srcIndex]<<" "<<weights[srcIndex]/TotalWeight<<endl;
        weights.push_back(pow(10,params[0]));
        if (OPT_LOADSIGNAL) {
	    PowerLawFlux pflux(1,spectralIndex);
            for(int arki=0;arki<ark.size();arki++){ 
                evLoader[arki]->LoadSourcePdf_Kent(srcLocations[srcIndex], srcSigmas[srcIndex]);
                evLoader[arki]->LoadSourceEvents(sourceEvents);
                I3PointGenerator i3point(sourceEvents, pflux, srcLocations[srcIndex], ark[arki]->livetime);
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
    
   for (int i=0; i<srcNames.size(); i++) {
        // For equal weights
        pdfEnhancement.push_back(weights[i]/TotalWeight); // Identical, i.e. no src enhancments
        //cout << pdfEnhancement[i] << " " << weights[i] << endl;
    }
 
    // At least make sure these are the same length as number  f s urces
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


}
