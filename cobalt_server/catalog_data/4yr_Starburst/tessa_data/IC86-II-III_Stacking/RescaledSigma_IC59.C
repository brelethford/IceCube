#include "TFile.h"
#include "TSpline.h"

TSpline5 *spline_IC59;

void loadSplines_IC59(){
    TFile *f=new TFile("$LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SplineMPEPullCorrectionIC59.root");
    spline_IC59=(TSpline5*)f->Get("Spline5");
    f->Close();
}


double RescaledSigma_IC59(double Sigma, double Energy)
{
    double x = log10(Energy);
    if (x<spline_IC59->GetXmin()) x=spline_IC59->GetXmin();
    else if (x>spline_IC59->GetXmax()) x=spline_IC59->GetXmax();    
    return Sigma * spline_IC59->Eval(x)/ 1.1774;
}  

