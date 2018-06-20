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

import SimpleMultiAnalysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab for MESE',prefix_chars='@')
    parser.add_argument('@ra', required=True,  type=float, help='RA of the source in degrees')
    parser.add_argument('@dec',required=True,  type=float, help='DEC of the source in degrees')
    parser.add_argument('@o', required=True,  type=str, help='output dir, if set to "auto" it will be IC86II_IV_MESEonly_allSkyTimeDep_disco_RA_DEC')
    parser.add_argument('@plot', action='store_true' )
    args = parser.parse_args()
    
    source_ra=args.__dict__["ra"]
    source_dec=args.__dict__["dec"]
    outdir=args.__dict__["o"]
    if not args.__dict__["plot"]:
        gROOT.SetBatch()

    print "search coordinates",source_ra,source_dec
    
    if outdir=="auto":
        outdir="IC86II_IV_MESEonly_allSkyTimeDep_TS_RA_"+str(source_ra)+"_DEC_"+str(source_dec)
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
    #gROOT.ProcessLine(".L SimpleMultiAnalysis.C")

    #for MESE    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+")
    ROOT.loadSplines_IC79_86_I_to_IV_SplineMPE_MESE()
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE_lc.C")
    gROOT.ProcessLine("I3Ark tark;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC86_II_to_IV_MESE_tdep_lc_folded.C(tark, OPT_USEREALDATA,\"SplineMPE\","+str(NbinsZen)+","+str(NbinsAz)+")")

    c_dataAZ_MESE=TCanvas()
    c_dataAZ_MESE.SetTitle("dataAZ_MESE")
    dataAZ_MESE=ROOT.tark.lcBkgProb.dataAZ
    dataAZ_MESE.SetTitle("MESE local coord raw data;azimuth folded;-cos(zenith)")
    dataAZ_MESE.Draw("col")
    c_dataAZ_MESE.SaveAs(outdir+"/dataAZ_MESE.eps")

    """    if not args.__dict__["plot"]:
    c_pdfAZ_MESE=TCanvas()    
    c_pdfAZ_MESE.SetTitle("pdfAZ_MESE")
    pdfAZ_MESE=ROOT.tark.lcBkgProb.pdfAZ
    pdfAZ_MESE.Draw("col")
    """
    naz=500
    nze=500
    pdfAZ_MESE_llh=TH2D("pdfAZ_MESE_llh","MESE local coord pdf;azimuth folded;-cos(zenith)",naz,0,60,nze,-1,0)

    for a in range(1,naz):
        for z in range(1,nze):
            pdfAZ_MESE_llh.SetBinContent(a, z, ROOT.tark.lcBkgProb.BackgroundLCProbA_mcozZ_folded(pdfAZ_MESE_llh.GetXaxis().GetBinCenter(a),pdfAZ_MESE_llh.GetYaxis().GetBinCenter(z)))
    
    c_pdfAZ_MESE_llh=TCanvas()   
    c_pdfAZ_MESE_llh.SetTitle("dataAZ_MESE_llh")
    pdfAZ_MESE_llh.Draw("col")
    c_pdfAZ_MESE_llh.SaveAs(outdir+"/c_pdfAZ_MESE_llh.eps")

    c=TCanvas()   
    c.Divide(1,NbinsZen)
    h=[]
    for i in range(1,NbinsZen+1):
        c.cd(i)
        h.append(pdfAZ_MESE_llh.ProjectionX("px"+str(i),nze*(i-1)/NbinsZen,nze*i/NbinsZen,""))
        h[-1].Scale(NbinsZen/float(nze))
        h[-1].Draw()
    if NbinsZen==1:
        ROOT.tark.lcBkgProb.pdfAZ_1D.Draw("pl same")

    c.SaveAs(outdir+"/pdfAZ_MESE_slices.eps")
    raw_input("Press Enter to continue...")

    newllh86II_to_IV=ROOT.NewLlhGausTime()
    newllh86II_to_IV.SetUseEnergy(True);
    #newllh86II_to_IV.SetOptimizeAngleDeg(10.)
    newllh86II_to_IV.SetOptimizeTolerance(0.01)
    newllh86II_to_IV.SetMonitorLevel(0)
    newllh86II_to_IV.SetEMaxRatioWarnOnlyOnce(1)
    newllh86II_to_IV.close_ = 10.
    newllh86II_to_IV.JimsTerm_ = True
    newllh86II_to_IV.SpectralPenalty_ = False
    newllh86II_to_IV.ndof_ = 3.
    newllh86II_to_IV.SetLivetime(ROOT.tark.livetime/86400.)
    newllh86II_to_IV.SetLocalCoordBkgProb(ROOT.tark.lcBkgProb,True) #true is there for the folded option

    testSearch=ROOT.EquatorialDeg(source_ra, source_dec) #if you want to macth Mike's previous discovery plots use declination = 16 deg

    spectralIndex = -2.
    tmean = (ROOT.tark.tmax + ROOT.tark.tmin)/2.

    tPdf = ROOT.GaussianTimePdf(ROOT.tark.tmin, ROOT.tark.tmax, tmean, 1e-5, 1.)
    ROOT.tark.SetPointSource(testSearch, ROOT.PowerLawFlux(1.,spectralIndex), tPdf)
    ROOT.tark.psData.GetSource().SetTimeAzBins( ROOT.tark.lcBkgProb.nbAz *6)#times 6 because of the folding
    newllh86II_to_IV.SetTimeBounds(tPdf)
    newllh86II_to_IV.SetAnalysis(ROOT.tark.psData, testSearch)

    sa=SimpleMultiAnalysis.SimpleMultiAnalysis()

    sa.Execute(ROOT.tark, newllh86II_to_IV, 1000, 0)

    sa.Write(outdir+"/TS1000.root","RECREATE")

    if args.__dict__["plot"]:
        raw_input("Press Enter to continue...")
    

    