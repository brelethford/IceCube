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
from ctypes import c_double
from math import sqrt,log10

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
        outdir="IC86II_IV_MESEonlyNonAzResetSigleSampleScript_allSkyTimeDep_disco_RA_"+str(source_ra)+"_DEC_"+str(source_dec)
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
    gROOT.ProcessLine(".L SetDisco.C")
    gROOT.ProcessLine(".L DetectionStudy.C")

    #for MESE    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+")
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE_lc.C")
    gROOT.ProcessLine("I3Ark tark;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC86_II_to_IV_MESE_tdep_lc_folded.C(tark, OPT_USEREALDATA,\"SplineMPE\","+str(NbinsZen)+","+str(NbinsAz)+")")

    c_dataAZ_MESE=TCanvas()
    c_dataAZ_MESE.SetTitle("dataAZ_MESE")
    dataAZ_MESE=ROOT.tark.lcBkgProb.dataAZ
    dataAZ_MESE.SetTitle("MESE local coord raw data;azimuth folded;-cos(zenith)")
    dataAZ_MESE.Draw("col")
    c_dataAZ_MESE.SaveAs(outdir+"/dataAZ_MESE.eps")

    naz=500
    nze=500
    pdfAZ_MESE_llh=TH2D("pdfAZ_MESE_llh","MESE local coord pdf;azimuth folded;-cos(zenith)",naz,0,60,nze,-1,0)

    for a in range(1,naz):
        for z in range(1,nze):
            pdfAZ_MESE_llh.SetBinContent(a, z, ROOT.tark.lcBkgProb.BackgroundLCProbA_mcozZ(pdfAZ_MESE_llh.GetXaxis().GetBinCenter(a),pdfAZ_MESE_llh.GetYaxis().GetBinCenter(z)))
    
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
        h[-1].Draw()
    c.SaveAs(outdir+"/pdfAZ_MESE_slices.eps")
    
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


    disco=ROOT.DiscoveryPotential()
  
    #parameters:  loops, optMedianUpperLimit,  significance,  power);
    ROOT.SetDisco(disco, 50, False, 2.87e-7, 0.5)
  
    dsigmas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 1., 3., 5., 10, 1e2]
    fscale =[]
    n_disco=[]
    c_fscale=c_double(0.)
    gn_disco = TGraph(len(dsigmas))
    gn_disco.SetNameTitle("disco","Discovery pot.;log_{10}(sigma [days]);n_{src}")
    
    gs_disco = TGraph(len(dsigmas))
    gs_disco.SetNameTitle("fscale","Discovery pot.;log_{10}(sigma [days]);fscale")

    for i in range(len(dsigmas)):
        print "processing sigma", dsigmas[i],"for disco"
        tPdf2 = ROOT.GaussianTimePdf(ROOT.tark.tmin, ROOT.tark.tmax, tmean,dsigmas[i],1.)
        ROOT.tark.SetPointSource(testSearch, ROOT.PowerLawFlux(1.,spectralIndex), tPdf2)
        nAzReset=ROOT.tark.lcBkgProb.nbAz *6#times 6 because of the folding
        """
        if dsigmas[i] < 8e-2:
            nAzReset=int(360./sqrt(dsigmas[i]))
            print "setting TimeAzBins to",nAzReset
        """
        ROOT.tark.psData.GetSource().SetTimeAzBins( nAzReset) 
        newllh86II_to_IV.SetTimeBounds(tPdf2)
        n_disco.append( ROOT.DetectionStudy_d(ROOT.tark, newllh86II_to_IV, disco,c_fscale))
        fscale.append(c_fscale.value)
        print "sigma",dsigmas[i],"n_disco",n_disco[i],"fscale",fscale[i]
        gn_disco.SetPoint(i,log10(dsigmas[i]),n_disco[i])
        gs_disco.SetPoint(i,log10(dsigmas[i]),fscale[i])
   
   
    ROOT.SetDisco(disco, 50, True, 0.5, 0.9)
  
    fscale_sens =[]
    n_sens=[]

    gn_sens = TGraph(len(dsigmas))
    gn_sens.SetNameTitle("sens","Discovery pot.;log_{10}(sigma [days]);n_{src}")
    
    gs_sens = TGraph(len(dsigmas))
    gs_sens.SetNameTitle("fscale_sens","Discovery pot.;log_{10}(sigma [days]);fscale")

    for i in range(len(dsigmas)):
        print "processing sigma", dsigmas[i],"for sens"
        tPdf2 = ROOT.GaussianTimePdf(ROOT.tark.tmin, ROOT.tark.tmax, tmean,dsigmas[i],1.)
        ROOT.tark.SetPointSource(testSearch, ROOT.PowerLawFlux(1.,spectralIndex), tPdf2)
        nAzReset=ROOT.tark.lcBkgProb.nbAz *6#times 6 because of the folding
        """
        if dsigmas[i] < 8e-2:
            nAzReset=int(360./sqrt(dsigmas[i]))
            print "setting TimeAzBins to",nAzReset
        """
        ROOT.tark.psData.GetSource().SetTimeAzBins( nAzReset) 

        newllh86II_to_IV.SetTimeBounds(tPdf2)
        n_sens.append( ROOT.DetectionStudy_d(ROOT.tark, newllh86II_to_IV, disco,c_fscale))
        fscale_sens.append(c_fscale.value)
        print "sigma",dsigmas[i],"n_sens",n_sens[i],"fscale",fscale_sens[i]
        gn_sens.SetPoint(i,log10(dsigmas[i]),n_sens[i])
        gs_sens.SetPoint(i,log10(dsigmas[i]),fscale_sens[i])
   
    fgr=TFile(outdir+"/disco_sens.root","recreate")
    gn_disco.Write()
    gs_disco.Write()
    gn_sens.Write()
    gs_sens.Write()
    fgr.Close()
    
    outf=open(outdir+"/disco_sens.dat","w")
    outf.write("sigma[days] ")
    for i in range(len(dsigmas)):
        outf.write(str(dsigmas[i])+" ")
    outf.write("\n")
    outf.write("disco_ns ")
    for i in range(len(dsigmas)):
        outf.write(str(n_disco[i])+" ")
    outf.write("\n")
    outf.write("disco_fs ")
    for i in range(len(dsigmas)):
        outf.write(str(fscale[i])+" ")
    outf.write("\n")
    outf.write("sens_ns ")
    for i in range(len(dsigmas)):
        outf.write(str(n_sens[i])+" ")
    outf.write("\n")
    outf.write("sens_fs ")
    for i in range(len(dsigmas)):
        outf.write(str(fscale_sens[i])+" ")
    outf.write("\n")
    outf.close()
    if args.__dict__["plot"]:
        raw_input("Press Enter to continue...")
    

    