#include "TFile.h"
#include "TSpline.h"

TSpline5 *spline_IC79;

void loadSplines_IC79(){
    TFile *f=new TFile("$LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SplineMPEPullCorrectionIC79.root");
    spline_IC79=(TSpline5*)f->Get("Spline5");
    f->Close();
}


double RescaledSigma_IC79(double Sigma, double Energy)
{
    double x = log10(Energy);
    if (x<spline_IC79->GetXmin()) x=spline_IC79->GetXmin();
    else if (x>spline_IC79->GetXmax()) x=spline_IC79->GetXmax();
    return Sigma * spline_IC79->Eval(x)/ 1.1774;

}  

