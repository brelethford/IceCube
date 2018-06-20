#include "TFile.h"
#include "TSpline.h"

TSpline5 *spline_IC86_SplineMPE;

void loadSplines_IC86_SplineMPE(){
    TFile *f=new TFile("$LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SplineMPEPullCorrectionIC86_I.root");
    spline_IC86_SplineMPE=(TSpline5*)f->Get("Spline5");
    f->Close();
}


double RescaledSigma_IC86_SplineMPE(double Sigma, double Energy)
{
    double x = log10(Energy);
    if (x<spline_IC86_SplineMPE->GetXmin()) x=spline_IC86_SplineMPE->GetXmin();
    else if (x>spline_IC86_SplineMPE->GetXmax()) x=spline_IC86_SplineMPE->GetXmax();    
    return Sigma * spline_IC86_SplineMPE->Eval(x)/ 1.1774;

}  

