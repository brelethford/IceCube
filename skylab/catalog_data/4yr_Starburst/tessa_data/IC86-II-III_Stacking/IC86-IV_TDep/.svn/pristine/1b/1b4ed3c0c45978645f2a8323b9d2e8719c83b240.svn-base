double mpfSigmaDegRescaledIC79V2(double mpfSigmaDegOld, double mmueEn)
{
  //IC79 MPE sigma correction
  double x = log10(mmueEn);
  
  return mpfSigmaDegOld * (72.793 -69.4918 * pow(x,1.) + 26.17 * pow(x,2.) -4.78732 * pow(x,3.) + 0.429997 * pow(x,4.) -0.0151473 * pow(x, 5.));
}
