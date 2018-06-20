double mpfSigmaDegRescaledIC59(double mpfSigmaDegOld, double mmueEn)
{
  //IC59 MPE sigma correction
  double x = log10(mmueEn);
  
  return mpfSigmaDegOld*(31.9051 - 24.5591*pow(x,1.)+7.19661*pow(x,2.)-0.908188*pow(x,3.)+0.043108*pow(x,4.));
}
