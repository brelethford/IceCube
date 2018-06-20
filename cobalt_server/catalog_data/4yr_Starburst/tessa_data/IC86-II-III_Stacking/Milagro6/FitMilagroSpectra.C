void FitMilagroSpectra(char* sourceName, double *params)
{
  cout<<"test 1 "<<endl;
  char* filename = gSystem->ExpandPathName("MGRO_1908+06.txt");
  ifstream in;
  in.open(filename);
  if (!in.good()) {
    cout << "Bad Milagro stacking spectrum txt file!\n";
    exit(1);
  }

  const int nSamples = 123;

  double xbuf=0, ybuf=0;
  double x[nSamples+1], y[nSamples+1];
  int i=0;

  //char *sourceName = "MGRO_1908+06";
  double sourceScale = 1.;

  cout<<"test 2 "<<endl;
  // Used GraphClick to extract J1908+06, all others scaled to this
  // Used Halzen, Kappes, O'Murchada (2008)
  if (sourceName == "MGRO J2019+37") {
    sourceScale = 1.06059;
  }
  if (sourceName == "MGRO J1908+06") {
    sourceScale = 1.; // This is source from txt file
  }
  if (sourceName == "MGRO J2031+42") {
    sourceScale = 0.72272;
  }
  if (sourceName == "MGRO J2043+36") {
    sourceScale = 0.51467;
  }
  if (sourceName == "MGRO J2032+37") {
    sourceScale = 0.39460;
  }
  // Taken from text of paper (not in plot) "about 2.5 times more than MGRO J2019+37"
  if (sourceName == "MGRO J1852+01") {
    sourceScale = 2.65148;
  }

  cout<<"test 3 "<<endl;
  while (1){
    in >> xbuf >> ybuf;
    if ( !in.good() ) break;
    // Convert E^2*Phi to Phi, and TeV to GeV
    y[i] = log10( ybuf*sourceScale/(xbuf**2.)/1.e3 );
    x[i] = log10( xbuf*1.e3 ); // Convert TeV to GeV// impose lower cutoffs
    //cout << x[i] << " " << y[i] << endl;
    i++;
  }

  cout<<"test 4 "<<endl;
  TGraph *gSpec = new TGraph(nSamples, x, y);

  //TF1 *fCutoff = new TF1("fCutoff","[0]*pow(x,[1])*exp(-x/[2])");
  // Fit in log, much more stable!
  TF1 *fCutoffLog = new TF1("fCutoffLog","[0]+[1]*x-pow(10,x)/pow(10,[2])*log10(exp(1))",2,5);
  fCutoffLog->SetParameters(-5, -2.1, 5);
  gSpec->Fit(fCutoffLog);
  fCutoffLog->GetParameters(params);

  cout<<"test 5 "<<endl;
  // Now in linear, to pass to PS code, cast to use same params
  ///TF1 *fCutoff = new TF1("fCutoff","[0]*pow(x,[1])*exp(-x/[2])");
  //TF1 *fCutoff = new TF1("fCutoff","pow(10,[0]+[1]*log10(x)-pow(10,log10(x))/pow(10,[2])*log10(exp(1)))",1.e3,1.e6);
  //fCutoff->SetParameters(params);

  //gSpec->Draw("AP");

  //TCanvas *can = new TCanvas("can","can",600,500);
  //can->SetLogy(1);
  //can->SetLogx(1);
  //fCutoff->Draw();

}
