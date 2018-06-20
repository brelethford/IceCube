
class SimpleMultiBlockAnalysis {
 public:
  TH1D hTestStatistic;
  TH1D hLogProb;
  TH1D hPar;
  TH2D hGammaNs;
  TH2D hThreshNs;
  TH1D hLogPInt;
  TCanvas* can;
  int ndof_;

  SimpleAnalysis() : can(NULL) , ndof_(2) { }
  ~SimpleAnalysis() { }

  void Usage();
  void SetNDoF(int ndof) { ndof_ = ndof; }
  void Execute(Ark& ark, AnalysisFn& llh, int nTrials, int nSrcEvents);
  bool Write(char *filename, char* fileoption);
};


void SimpleMultiBlockAnalysis::Usage() {
  cout << "SimpleAnalysis sa;\n";
  cout << "sa.SetNDoF(int ndof);   // default degrees of freedom is 2\n";
  cout << "sa.Execute(Ark& ark, Analysis& llh,\n"
          "           int nTrials, int nSrcEvents = 0)\n";
}


void SimpleMultiBlockAnalysis::Execute(MultiArk& ark, MultiBlockAnalysisFn& llh, 
			     string blocksFile, int nTrials, int nSrcEvents)
{

// I think the easiest way to do this is to set everything before
//  llh.SetAnalysisSet(ark.psData);
//  llh.SetSearchCoord(ark.mySrcLocation);

  // find out if energy is being used (changes # of degrees of freedom)
  //if (ndof_<1 || ndof_>2) {
  //  cout << "ndof==" << ndof_ << " not allowed right now.  Setting to 1.\n";
  //  ndof_ = 1;
  //}

  hTestStatistic.Reset();
  hLogProb.Reset();
  hPar.Reset();
  hGammaNs.Reset();
  hLogPInt.Reset();
  hThreshNs.Reset();

  hTestStatistic.SetBins(100,0,1);
  hTestStatistic.SetTitle(";2 ln #lambda;trials");
  hLogProb.SetBins(100,-1,1);
  hLogProb.SetTitle(";log_{10} Prob;trials");
  hPar.SetBins(100,-1,1);
  hPar.SetTitle(";nSrcBest; trials");
  hGammaNs.SetBins(50,-1,1,30,1,4);
  hGammaNs.SetTitle(";nSrcBest;gammaBest");
  
  hThreshNs.SetBins(20,0,20,50,0.,2.5e-5);
  hThreshNs.SetTitle(";nSrcBest;threshold_Best");

  // Let Root figure out how to extend the bin range automatically
  hTestStatistic.SetBit(TH1::kCanRebin);
  hLogProb.SetBit(TH1::kCanRebin);
  hPar.SetBit(TH1::kCanRebin);
  hGammaNs.SetBit(TH1::kCanRebin);

  //
  // MAIN LOOP
  //

  CountMonitor countMon(10.,nTrials);
  cout << "going for " << nTrials << endl;
  
  //double laglimit = 0.5;
  double threshMax = 
     ( GetHighestBlock(blocksFile) + GetSecondHighestBlock(blocksFile) ) /2.;
  
  //MinuitParDef pd = MinuitParDef("nSrc",2,0.1, 0,100);
  
  //double guessThresh;  
  //double lagGuess;//SearchForLag(llh, laglimit);
  //maf.nPar = 4;
      
  vector<MinuitParDef> pdv;// (4,pd);
  //pdv[0] = 
  //pdv.push_back( MinuitParDef("nSrc",2,0.1, 0,100) );
  //pdv[1] = MinuitParDef("gamma",2.5,0.5, 1., 4.);
  //pdv.push_back( MinuitParDef("gamma",2.5,0.5, 1., 4.) );
  //pdv[2] = MinuitParDef("lag", lagGuess, laglimit/10., -1.0*laglimit, laglimit);
  //pdv.push_back( MinuitParDef("lag", lagGuess, laglimit/10., -1.0*laglimit, laglimit) );
  //pdv[3] = MinuitParDef("thresh",guessThresh, guessThresh/10., 0., threshMax);
  //pdv.push_back( MinuitParDef("thresh",guessThresh, guessThresh/10., 0., threshMax) );
  //llh.SetParDefs(pdv);
  
  for (int i=0; i<nTrials; ++i) {
  
    //pdv.clear();
    //pdv.push_back( MinuitParDef("nSrc",2.5,1., 0.,100) );
    //pdv.push_back( MinuitParDef("gamma",2.5,0.5, 1., 4.) );
  
    ark.psData->GenerateDataSet_with_nSrcEvents(nSrcEvents);

    cout << "got dataset" << endl;

    //llh.PrepareAnalysis();

    //double lagGuess = llh.SearchForLag(laglimit);
    //pdv.push_back( MinuitParDef("lag", lagGuess, laglimit/10., -1.0*laglimit, laglimit) );

    //llh.SearchBlockSpace(blocksFile, lagGuess, guessThresh, threshMax);
    //pdv.push_back( MinuitParDef("thresh",guessThresh, threshMax/10., 0., threshMax) );
    //cout << "Guess Threshold: " << guessThresh << endl;
    //llh.SetParDefs(pdv);

    llh.MaximizeLlh();

    hTestStatistic.Fill(llh.GetTestStatistic() * 2.);
    hLogProb.Fill(log10(llh.GetEstProb()));
    hPar.Fill(llh.GetPar(0));
    if (ndof_>1) {
      hGammaNs.Fill(llh.GetPar(0),llh.GetPar(1));
      hThreshNs.Fill(llh.GetPar(0),llh.GetPar(3));
    }
    countMon.UpdateCount();
  } //*/

  //
  // PLOT OUTPUT
  //

  // only make new canvas if it didn't already exist:
  //
    can = dynamic_cast<TCanvas*> 
      ( gROOT->GetListOfCanvases()->FindObject("canSimpleAnalysis") );
  if ( !can ) {
    can =  new TCanvas("canSimpleAnalysis","canSimpleAnalysis",20,20,950,650);
    can->Divide(3,2,0.005,0.005);
//    for (int i=0; i<6; ++i) {
//      can->cd(i);
//      gPad->SetRightMargin(0.01);
//      gPad->SetTopMargin(0.01);
//    }
  }

  can->cd(1);

  TF1 *fchi2 = new TF1("fchi2","[0] * (1./sqrt(2*TMath::Pi())) "
		       "* exp(-x/2)/sqrt(x)",0,50);
  double xRange = hTestStatistic.GetXaxis()->GetXmax() -
    hTestStatistic.GetXaxis()->GetXmin();
  double binsPerUnit = hTestStatistic.GetNbinsX() / xRange;
  double normFactor = (nTrials/2.) / binsPerUnit; // to match histogram plot
  // divide by 2 because we only count excess tail contribution to chi2 
  fchi2->SetParameter(0,normFactor);
  fchi2->SetLineColor(kBlue);

  hTestStatistic.Draw();
  if (ndof_==1) { fchi2->Draw("same"); }

  if (ndof_==2) {
    TF1 *fchi2_2d = new TF1("fchi2_2d","[0] * (1./2.) * exp(-x/2)",0,50);
    fchi2_2d->SetParameter(0,normFactor);
    fchi2_2d->SetLineColor(kRed);
    fchi2_2d->Draw("same"); 
  }
  
  if (ndof_== 3) {
    TF1 *fchi2_3d = new TF1("fchi2_2d","[0] * [1] * ( sqrt([1]*x) * exp(-[1]*x/2) )",0,50);
    fchi2_3d->SetParameter(0,normFactor);
    fchi2_3d->SetParameter(1,2); // correction for chi = 2 * log lambda
    fchi2_3d->SetLineColor(kRed);
    fchi2_3d->Draw("same"); 
  }

  
  gPad->SetLogy(1);


  can->cd(2);
  hLogProb.Draw();
  cout << "Median Log Prob of Max Like: " << histMedian(&hLogProb) << endl;

  can->cd(3);
  gPad->SetLogy(1);
  hLogPInt = hLogProb;  // first make binning the same
  AscendingCumulate(&hLogProb, &hLogPInt);
  hLogPInt.Draw("e");
  double entries = hLogProb.GetEntries();
  TLine *line = new TLine(log10(1./entries),1.,0.,entries);
  line->Draw();

  can->cd(4);
  hPar.Draw();
  gPad->SetLogy(1);

  can->cd(5);
  if (ndof_>2) {
    hThreshNs.Draw("colz");
  }
}

bool SimpleMultiBlockAnalysis::Write(char*filename, char* fileoption) {
  TFile *fileOutput = new TFile(filename, fileoption);
  if (fileOutput->IsZombie()) {
    cout << "Try using 'recreate' option?\n";
    return false; // no file saved
  }
  
  cout << "Writing Histograms  to: " << filename << endl;

  hTestStatistic.Write("hTestStatistic");
  hLogProb.Write("hLogProb");
  hLogPInt.Write("hLogPInt");
  hPar.Write("hPar");
  hGammaNs.Write("hGammaNs");
  can->Write("can");

  if(ndof_ == 2)fchi2_2d->Write("fchi2_2d");
  if(ndof_ == 1)fchi2->Write("fchi2_2d");
  if(ndof_ == 3)fchi2_3d->Write("fchi2_3d");

  fileOutput->Close();
  return true;  // file saved (apparently!)                                                                                                                                                                 
}
