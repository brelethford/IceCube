#include "TFile.h"
#include "TSpline.h"

TSpline5 *upspline;
TSpline5 *downspline;

void loadSplines_IC86_II_III_IV_SplineMPE(){
    TFile *fup=new TFile("SplineMPEPullCorrectionIC86_II_III_IV_upgoing.root");
    upspline=(TSpline5*)fup->Get("Spline5");
    fup->Close();
    TFile *fdown=new TFile("SplineMPEPullCorrectionIC86_II_III_IV_downgoing.root");
    downspline=(TSpline5*)fdown->Get("Spline5");
    fdown->Close();
}

double RescaledSigma_IC86_II_III_IV_SplineMPE(double Sigma, double Energy, double Zenith) {
    double x = log10(Energy);
    
    if (Zenith > 85.){
        if (x<upspline->GetXmin()) x=upspline->GetXmin();
        else if (x>upspline->GetXmax()) x=upspline->GetXmax();
        return Sigma * upspline->Eval(x)/ 1.1774;
    } else {
        if (x<downspline->GetXmin()) x=downspline->GetXmin();
        else if (x>downspline->GetXmax()) x=downspline->GetXmax();
        return Sigma * downspline->Eval(x)/ 1.1774;
    }  
}
