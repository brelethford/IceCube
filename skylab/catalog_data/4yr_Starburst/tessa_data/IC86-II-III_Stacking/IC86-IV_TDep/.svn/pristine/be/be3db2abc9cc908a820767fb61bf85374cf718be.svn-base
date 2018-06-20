#include "TFile.h"
#include "TSpline.h"

TSpline5 *spline_IC40;

void loadSplines_IC40(){
    TFile *f=new TFile("$LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SplineMPEPullCorrectionIC40.root");
    spline_IC40=(TSpline5*)f->Get("Spline5");
    f->Close();
}


double RescaledSigma_IC40(double Sigma, double Energy)
{
    double x = log10(Energy);
    if (x<spline_IC40->GetXmin()) x=spline_IC40->GetXmin();
    else if (x>spline_IC40->GetXmax()) x=spline_IC40->GetXmax();
    return Sigma * spline_IC40->Eval(x)/ 1.1774;

}  

