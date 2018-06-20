double mpfSigmaDegRescaledIC79Sirin(double mpfSigmaDegOld, double nch)
{
  //IC79 MPE sigma correction
  
  return mpfSigmaDegOld * (5.620-7.448*log10(nch)+3.519*pow(log10(nch), 2.0)+0.252*pow(log10(nch), 3.0)-0.509*pow(log10(nch), 4.0)+0.095*pow(log10(nch), 5.0));

}
