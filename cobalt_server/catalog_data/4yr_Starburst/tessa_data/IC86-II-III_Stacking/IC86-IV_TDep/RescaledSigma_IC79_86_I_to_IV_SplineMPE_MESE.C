#include "TFile.h"
#include "TSpline.h"

TSpline5 *spline_IC79_86_I_to_IV_SplineMPE_MESE;

void loadSplines_IC79_86_I_to_IV_SplineMPE_MESE(){
    TFile *f=new TFile("$LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/SplineMPEPullCorrectionMESE.root");
    spline_IC79_86_I_to_IV_SplineMPE_MESE=(TSpline5*)f->Get("Spline5");
    f->Close();
}


double RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(double Sigma, double Energy)
{
    double x = log10(Energy);
    if (x<spline_IC79_86_I_to_IV_SplineMPE_MESE->GetXmin()) x=spline_IC79_86_I_to_IV_SplineMPE_MESE->GetXmin();
    else if (x>spline_IC79_86_I_to_IV_SplineMPE_MESE->GetXmax()) x=spline_IC79_86_I_to_IV_SplineMPE_MESE->GetXmax();
    return Sigma * spline_IC79_86_I_to_IV_SplineMPE_MESE->Eval(x)/ 1.1774;
}  

