#!/usr/bin/env python
import ROOT,sys
ROOT.PyConfig.IgnoreCommandLineOptions = True
gROOT = ROOT.gROOT
TCanvas = ROOT.TCanvas
TH2D = ROOT.TH2D
TH1D = ROOT.TH1D
TGraph =ROOT.TGraph
TFile = ROOT.TFile

import argparse,os
from math import sqrt,log10

import numpy as np
import SimpleMultiAnalysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab for MESE',prefix_chars='@')
    parser.add_argument('@ra', required=True,  type=float, help='RA of the source in degrees')
    parser.add_argument('@dec',required=True,  type=float, help='DEC of the source in degrees')
    parser.add_argument('@o', required=True,  type=str, help='output dir, if set to "auto" it will be IC86II_IV_MESE_PSsample_allSkyTimeDep_disco_RA_DEC')
    parser.add_argument('@plot', action='store_true' )
    args = parser.parse_args()
    
    source_ra=args.__dict__["ra"]
    source_dec=args.__dict__["dec"]
    outdir=args.__dict__["o"]
    if not args.__dict__["plot"]:
        gROOT.SetBatch()

    print "search coordinates",source_ra,source_dec
    
    if outdir=="auto":
        outdir="IC86II_IV_MESE_PSsample_allSkyTimeDep_TS_RA_"+str(source_ra)+"_DEC_"+str(source_dec)
    try:
        os.makedirs(outdir)
    except Exception as e:
        if e.errno == 17:
            print outdir+" alredy exists"
        else:
            sys.exit(e)

    ROOT.OPT_USEREALDATA=False
    
    gROOT.Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C")
    
    ROOT.initialize_ran1(-156)
    NbinsZen=1
    NbinsAz=10
    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/ArkTime.C+")
    
    #for MESE    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+")
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE_lc.C")
    gROOT.ProcessLine("I3Ark MESEark;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC86_II_to_IV_MESE_tdep_lc_folded.C(MESEark, OPT_USEREALDATA,\"SplineMPE\","+str(NbinsZen)+","+str(NbinsAz)+")")
    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_II_III_IV_SplineMPE.C+")
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC86_II_III_IV_lc.C")
    gROOT.ProcessLine("I3Ark psark;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic86_II_III_IV_tdep_lc.C(psark, OPT_USEREALDATA,\"SplineMPE\")")
    """
    c_dataAZ_MESE=TCanvas()
    c_dataAZ_MESE.SetTitle("dataAZ_MESE")
    dataAZ_MESE=ROOT.MESEark.lcBkgProb.dataAZ
    dataAZ_MESE.SetTitle("MESE local coord raw data;azimuth folded;-cos(zenith)")
    dataAZ_MESE.Draw("col")
    c_dataAZ_MESE.SaveAs(outdir+"/dataAZ_MESE.eps")

    c_dataAZ_PS=TCanvas()
    c_dataAZ_PS.SetTitle("dataAZ_PS")
    dataAZ_PS=ROOT.psark.lcBkgProb.dataAZ
    dataAZ_PS.SetTitle("PS local coord raw data;azimuth folded;-cos(zenith)")
    dataAZ_PS.Draw("col")
    c_dataAZ_PS.SaveAs(outdir+"/dataAZ_PS.eps")

    naz=500
    nze=500
    pdfAZ_MESE_llh=TH2D("pdfAZ_MESE_llh","MESE local coord pdf;azimuth folded;-cos(zenith)",naz,0,60,nze,-1,0)
    pdfAZ_PS_llh=TH2D("pdfAZ_PS_llh","PS local coord pdf;azimuth folded;-cos(zenith)",naz,0,60,nze,-1,0)

    for a in range(1,naz):
        for z in range(1,nze):
            pdfAZ_MESE_llh.SetBinContent(a, z, ROOT.MESEark.lcBkgProb.BackgroundLCProbA_mcozZ_folded(pdfAZ_MESE_llh.GetXaxis().GetBinCenter(a),pdfAZ_MESE_llh.GetYaxis().GetBinCenter(z)))
            pdfAZ_PS_llh.SetBinContent(a, z, ROOT.psark.lcBkgProb.BackgroundLCProbA_mcozZ(pdfAZ_PS_llh.GetXaxis().GetBinCenter(a),pdfAZ_PS_llh.GetYaxis().GetBinCenter(z)))
    
    c_pdfAZ_MESE_llh=TCanvas()   
    c_pdfAZ_MESE_llh.SetTitle("dataAZ_MESE_llh")
    pdfAZ_MESE_llh.Draw("col")
    c_pdfAZ_MESE_llh.SaveAs(outdir+"/c_pdfAZ_MESE_llh.eps")
    
    c_pdfAZ_PS_llh=TCanvas()   
    c_pdfAZ_PS_llh.SetTitle("dataAZ_PS_llh")
    pdfAZ_PS_llh.Draw("col")
    c_pdfAZ_PS_llh.SaveAs(outdir+"/c_pdfAZ_PS_llh.eps")

    c=TCanvas()   
    c.Divide(1,NbinsZen)
    h=[]
    for i in range(1,NbinsZen+1):
        c.cd(i)
        h.append(pdfAZ_MESE_llh.ProjectionX("px"+str(i),nze*(i-1)/NbinsZen,nze*i/NbinsZen,""))
        h[-1].Scale(NbinsZen/float(nze))
        h[-1].Draw()
    if NbinsZen==1:
        print "drawing:"
        ROOT.MESEark.lcBkgProb.pdfAZ_1D.Draw("pl same")
    c.SaveAs(outdir+"/pdfAZ_MESE_slices.eps")
    """
    
    newllh86II_to_IV=ROOT.NewLlhGausTime()
    newllh86II_to_IV.SetUseEnergy(True);
    #newllh86II_to_IV.SetOptimizeAngleDeg(10.)
    newllh86II_to_IV.SetOptimizeTolerance(0.001)
    newllh86II_to_IV.SetMonitorLevel(0)
    newllh86II_to_IV.SetEMaxRatioWarnOnlyOnce(1)
    newllh86II_to_IV.close_ = 10.
    newllh86II_to_IV.JimsTerm_ = True
    newllh86II_to_IV.SpectralPenalty_ = False
    newllh86II_to_IV.ndof_ = 3.
    newllh86II_to_IV.SetLivetime(ROOT.MESEark.livetime/86400.)
    newllh86II_to_IV.SetLocalCoordBkgProb(ROOT.MESEark.lcBkgProb,True) #true is there for the folded option
    newllh86II_to_IV.SetAnalysisSet(ROOT.MESEark.psData)
    
    
    newllh86PSII_to_IV=ROOT.NewLlhGausTime()
    newllh86PSII_to_IV.SetUseEnergy(True);
    #newllh86PSII_to_IV.SetOptimizeAngleDeg(10.)
    newllh86PSII_to_IV.SetOptimizeTolerance(0.01)
    newllh86PSII_to_IV.SetMonitorLevel(0)
    newllh86PSII_to_IV.SetEMaxRatioWarnOnlyOnce(1)
    newllh86PSII_to_IV.close_ = 10.
    newllh86PSII_to_IV.JimsTerm_ = True
    newllh86PSII_to_IV.SpectralPenalty_ = False
    newllh86PSII_to_IV.ndof_ = 3.
    newllh86PSII_to_IV.SetLivetime(ROOT.psark.livetime/86400.)
    newllh86PSII_to_IV.SetLocalCoordBkgProb(ROOT.psark.lcBkgProb) 
    newllh86PSII_to_IV.SetAnalysisSet(ROOT.psark.psData)
    
    #ROOT.psark.psData.GetSource().SetTimeAzBins( ROOT.psark.lcBkgProb.nbAz)#times 6 because of the folding
    gROOT.ProcessLine("MultiArk mark;")
    ROOT.mark.AddArk(ROOT.MESEark)
    ROOT.mark.AddArk(ROOT.psark)          

    maf=ROOT.MultiGaussAnalysisFn()
    maf.AddAnalysisFn(newllh86II_to_IV)
    maf.AddAnalysisFn(newllh86PSII_to_IV)
    maf.SetTimeBounds(ROOT.MESEark.tmin,ROOT.MESEark.tmax )
    maf.SetTimeBounds(ROOT.psark.tmin,ROOT.psark.tmax )

    testSearch=ROOT.EquatorialDeg(source_ra, source_dec)
    tmean = (ROOT.MESEark.tmax + ROOT.MESEark.tmin)/2.
    tmean = (ROOT.psark.tmax + ROOT.psark.tmin)/2.
    #tPdf = ROOT.GaussianTimePdf(ROOT.MESEark.tmin, ROOT.MESEark.tmax, tmean, 1., 1.)
    tPdf = ROOT.GaussianTimePdf(ROOT.psark.tmin, ROOT.psark.tmax, tmean, 1e-5, 1.)
    newllh86II_to_IV.SetTimeBounds(tPdf)
    newllh86PSII_to_IV.SetTimeBounds(tPdf)
    ROOT.mark.SetPointSource(testSearch, ROOT.PowerLawFlux(1,-2), tPdf)
    maf.SetSearchCoord(testSearch)
    
    newllh86PSII_to_IV.SetAnalysis(ROOT.psark.psData, testSearch);
    newllh86II_to_IV.SetAnalysis(ROOT.MESEark.psData, testSearch);

    pt=ROOT.NewLlhGausTime_ParTranslator()
    pt.SetRange(1,4,31) #gamma_min, gamma_max, nBins
    gROOT.ProcessLine("MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark.psData);")

    pt.SetTranslator(ROOT.mas)

    maf.SetParTranslator(pt)
    
    sa=SimpleMultiAnalysis.SimpleMultiAnalysis()

    sa.Execute(ROOT.mark, maf, 50000, 0)

    sa.Write(outdir+"/TS.root","RECREATE")

    if args.__dict__["plot"]:
        raw_input("Press Enter to continue...")
    
    

    