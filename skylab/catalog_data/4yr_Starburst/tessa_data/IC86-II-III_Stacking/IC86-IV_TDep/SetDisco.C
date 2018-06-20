
// Default macro, gives user instructions

void SetDisco() {
  cout << "Usage:\n";
  cout << ".x SetDisco.C(DiscoveryPotential& disco,\n";
  cout << " int loops, bool optMedianUpperLimit = false,\n";
  cout << "              double significance = 2.87e-7, ";
  cout << " double power = 0.5)\n\n";
  cout << " Loops: 5 is fast test;  20 is typical;  100 is high accuracy\n\n";
  cout << " optMedianUpperLimit:\n";
  cout << "  false ==> Use signifance and power given";
  cout << " (e.g. 5 sigma disc. pot.)\n";
  cout << "  true  ==> significance will be calculated from median p-value of bkg trials\n";
  cout << "            and power is set to 0.9\n";
}


void SetDisco(DiscoveryPotential& disco, 
	      int loops,  bool optMedianUpperLimit = false,
	      double significance = 2.87e-7, // one-sided 5-sigma 
	      double power = 0.5)
{
  cout << "ark_SetDisco:\n";
  
  disco.monitor_ = true;

  disco.method_ = 2;
  disco.loops_ = loops;
  cout << "!! disco set with:\n   Loops: " << disco.loops_ << endl;

  if (optMedianUpperLimit) {
    int nBkgTrials = loops*50;  // scale number of trials to loops
    disco.SetForMedianUpperLimit(nBkgTrials);
    cout << "   Setting for Median Upper Limit with nBkgTrials: ";
    cout << nBkgTrials << endl;
    cout << "   Significance: T.B.D.";
  }
  else {
    disco.SetDetectionSignificance(significance);
    disco.SetDetectionPower(power);
    cout << "   Significance: " << disco.GetDetectionSignificance();
  }  

  cout << "   Power: " << disco.GetDetectionPower() << endl;
}
