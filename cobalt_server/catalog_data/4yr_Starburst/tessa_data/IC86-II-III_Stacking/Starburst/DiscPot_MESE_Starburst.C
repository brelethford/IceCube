{
    gROOT->Macro("$LAB_MAIN_DIR/llh/loadlibs.C");
    gInterpreter->AddIncludePath( gSystem->Getenv("LAB_MAIN_DIR") );

    if (!LOADSUCCESS) { return 0.; } // stop, and signal that script failed.
        
    initialize_ran1(-55);  // seed has to be a *NEGATIVE* integer

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/ArkExt.C+");

    bool OPT_USEREALDATA = false;
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/mpfSigmaDegRescaledIC40.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/mpfSigmaDegRescaledIC59.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/mpfSigmaDegRescaledIC79Sirin.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC86_SplineMPE.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/RescaledSigma_IC79_86_I_to_IV_SplineMPE.C+");
    
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
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/load_ark_ic59_BDT.C(*ark59, OPT_USEREALDATA)");
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

    int monLev = 5;
        
    printf("Setting up Stacking Analysis Sources \n");
    printf("Spectra are fit on the fly from Halzen, Kappes, O'Murchada (2008)\n");
        
        
    bool OPT_LOADSIGNAL = true;
    

    vector<I3Event> sourceEvents;
    sourceEvents.clear();
    vector<EquatorialDeg> srcLocations;
    srcLocations.clear();
    vector<char*> srcNames;
    srcNames.clear();
    vector<double> weights;
    weights.clear();
    vector<double> srcSigmas;
    srcSigmas.clear(); 
     
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
    int nStarburstListExpected = 127; // Just a precaution // Just a precaution
    int nStarburstList = 0;

    char* filename = "20mJy_sample_181208.csv";

    // Can configure analysis to use Flux at 60um as weight, defaults to false
    bool OPT_USE_THEORY_WEIGHT_SOURCE = 1;  // Source Simulation Weight
    bool OPT_USE_THEORY_WEIGHT_SEARCH = 1;  // Search Hypothesis Weight

    printf("Weight source simulation by 60um Flux: %i\n",
         OPT_USE_THEORY_WEIGHT_SOURCE);
    printf("Weight search hypothesis by 60um Flux: %i\n",
         OPT_USE_THEORY_WEIGHT_SEARCH);

    //FILE *fp = fopen(filename,"r");
    ifstream in;
    in.open(filename);

    char buffer[100];
    double raDeg, decDeg, Flux60um;
    double dummy=0;

// Theoretical enhancements for I3MultiSignalGenerator signal generation
    vector<double> enhancementFactors;
    enhancementFactors.clear(); 
 
    // Search enhancements for likelihood pdfs
    vector<double> pdfEnhancements;
    pdfEnhancements.clear();

    // Read the Starburst catalogue, store all info in vectors
    int srcIndex=0;
    //while( fscanf(fp,"%s %lf %lf\n", buffer, &raDeg, &decDeg) == 3) 
    // Read in header
    for (int i=0; i<21; i++) {
      in >> buffer;
    }
    while (1){
      in >> buffer >> raDeg >> decDeg >> dummy >> dummy >> dummy >> dummy >> dummy >> Flux60um >> dummy >> dummy >> dummy >>
 dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
      if ( !in.good() ) break;
      printf("%03d %14s %6.2f %6.2f %5.2f\n",srcIndex,buffer,raDeg,decDeg,Flux60um);
 
      srcLocations.push_back( EquatorialDeg(raDeg, decDeg) );
      srcSigmas.push_back(0.0); // all are point-like
 
      // These vectors are normalized a bit later
      if (OPT_USE_THEORY_WEIGHT_SOURCE) { enhancementFactors.push_back(Flux60um); }
      else { enhancementFactors.push_back(1.); }
      if (OPT_USE_THEORY_WEIGHT_SEARCH) { pdfEnhancements.push_back(Flux60um);}
      else { pdfEnhancements.push_back(1.); }
 
      ++srcIndex;
      cout << srcIndex << " " << flush;
    }
    nStarburstList = srcIndex;
 
    // At least make sure vector lengths are same size
    assert ( nStarburstList == nStarburstListExpected );
    assert ( nStarburstList == srcLocations.size());
    assert ( srcSigmas.size() == srcLocations.size() ); // else misconfigured
    assert ( enhancementFactors.size() == srcLocations.size() ); // else misconfigured
    assert ( pdfEnhancements.size() == srcLocations.size() ); // else misconfigured
 
    cout << nStarburstList << " total Starburst sources read.\n";
    //fclose(fp);
    in.close();
    cout<<"TEST 1"<<endl;
    cout<<srcSigmas.size()<<endl;
   

 // Normalize the enhancement factors (easier to interpret final flux)
    double sumEF=0;
    double sumPE=0;
    for (int i=0; i<srcSigmas.size(); i++){
      sumEF += enhancementFactors[i];
      sumPE += pdfEnhancements[i];
    }
    for (int i=0; i<srcSigmas.size(); i++){
      enhancementFactors[i]/=sumEF;
      pdfEnhancements[i]/=sumPE;
    }


    for(int srcIndex=0;srcIndex<srcSigmas.size();srcIndex++){
       cout<<srcIndex<<" "<< srcSigmas[srcIndex]<<" "<<srcLocations[srcIndex].GetRa() << " " << srcLocations[srcIndex].GetDec() <<" "<<pdfEnhancements[srcIndex]<<endl;
    }
   
    cout<<"TEST 2 "<<endl;
    cout<<srcSigmas.size()<<endl;
    double spectralIndex=-2.0;

    for(int srcIndex=0;srcIndex<srcSigmas.size();srcIndex++){ 
        // For the theoretical spectra, the differences in flux are built
        //  directly into the flux fits of Halzen, Kappes, O'Murchada's
        //  simulated spectral data directly.
        // Just give 1.0 for each source here.
        //enhancementFactors.push_back(weights[srcIndex]/TotalWeight);
        cout<<"Error weight test : "<<srcIndex<<endl;
        if (OPT_LOADSIGNAL) {
    //        TF1 *fCutoff=new TF1(TString("fCutoff")+TString(srcNames[srcIndex]),"pow(10,[0]+[1]*log10(x)-pow(10,log10(x))/pow(10,[2])*log10(exp(1)))",1.e3,1.e6);
     //       fCutoff->SetParameters(params);
    //        FormulaFlux formflux(*fCutoff);
              PowerLawFlux pflux(1,spectralIndex);
   //           FluxBase *myFluxPtr = &formflux;
            for(int arki=0;arki<ark.size();arki++){ 
                evLoader[arki]->LoadSourcePdf_Kent(srcLocations[srcIndex], srcSigmas[srcIndex]);
                evLoader[arki]->LoadSourceEvents(sourceEvents);
                I3PointGenerator i3point(sourceEvents, pflux, srcLocations[srcIndex], ark[arki]->livetime);
                msg[arki]->AddSignal(i3point,enhancementFactors[srcIndex]);
            }
        }
    }
    
    cout<<"TEST 3"<<endl;
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
    cout<<"TEST 4 "<<endl; 
    //pdfEnhancement.clear();
    double TotalWeight=0.0;    
    
    // At least make sure these are the same length as number of sources
    assert ( srcSigmas.size() == srcLocations.size() ); // else misconfigured
    assert ( enhancementFactors.size() == srcLocations.size() ); // else misconfigured
    assert ( pdfEnhancements.size() == srcLocations.size() ); // else misconfigured
    
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
                srcWeightsArray_tmp[srcIndex].push_back(srcWeight*pdfEnhancements[srcIndex]);
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

    // Disco (Discovery Potential and Sensitivity Estimator)
    DiscoveryPotential disco;
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/SetDisco.C");
    // parameters:  loops, optMedianUpperLimit,  significance,  power);
    //SetDisco(disco, 15, false, 2.87e-7, 0.5); // Discovery Potential
    SetDisco(disco, 30, true, 0.5, 0.9); // Sensitivity 

    gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-II-III_Stacking/MultiDetectionStudy.C");
     
    TCanvas *c = MultiDetectionStudy(*mark, maf, mas, disco); //


    c->SaveAs("Discpotwithextensions.png");
    c->SaveAs("Discpotwithextensions.root");
    c->SaveAs("Discpot.C");


}
