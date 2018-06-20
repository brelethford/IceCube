double RescaledSigma_IC862plus3_SplineMPE(double Sigma, double Energy, double Zenith)
{
  double x = log10(Energy);
  
  //MuEXUnc: return Sigma * (1.179820 - 0.348988 * pow(x,1.) + 0.0684238 * pow(x,2.));
  //return Sigma * (2.912332 - 1.500813 * pow(x,1.) + 0.265407 * pow(x,2.)); // CR (old)
  //return Sigma * (5.127137 - 2.301534 * pow(x,1.) + 0.3385327 * pow(x,2.)); // Pb (old)
  //return Sigma * (4.469995 - 5.20666386 * pow(x,1.) + 2.7868725 * pow(x,2.) - 0.6573967 * pow(x,3.) + 0.05888417897 * pow(x,4.)); // CR (new)
  //return Sigma * (15.851093 - 12.4445132 * pow(x,1.) + 3.7039592 * pow(x,2.) - 0.453491478 * pow(x,3.) + 0.0206292269 * pow(x,4.)); // Pb (new)
  
if (Zenith > 85.)  
{
  return Sigma * pow(10.,(1.176 - 0.834 * pow(x,1.) + 0.241 * pow(x,2.) - 0.026 * pow(x,3.) +0.001 * pow(x,4.))); // Pb (NCh)
}
else
{
  return Sigma * pow(10.,(1.637 - 1.265 * pow(x,1.) + 0.371 * pow(x,2.) - 0.042 * pow(x,3.) +0.002 * pow(x,4.))); // Pb (NCh)
}  
  
// Pb (real)
  //return Sigma * (4.09787 - 2.589478 * pow(x,1.) + 0.525268 * pow(x,2.) + 0.14046896 * pow(x,3.) - 0.013224799 * pow(x,4.)); // Pb (NCh)
}
