#!/usr/bin/env python
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
gROOT = ROOT.gROOT
TH2D  =ROOT.TH2D
TFile =ROOT.TFile
TTree =ROOT.TTree
TCanvas=ROOT.TCanvas

if __name__ == "__main__":
    ROOT.OPT_USEREALDATA=True
    gROOT.Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C")
    NbinsZen=1
    NbinsAz=12
    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/ArkTime.C+")

    #for MESE    
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+")
    ROOT.loadSplines_IC79_86_I_to_IV_SplineMPE_MESE()
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC79_86_I_to_IV_MESE_lc.C")
    gROOT.ProcessLine("I3Ark MESEark")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_IC86_II_to_IV_MESE_tdep_lc_1D_splines.C(MESEark, OPT_USEREALDATA,\"SplineMPE\","+str(NbinsZen)+","+str(NbinsAz)+")")

    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/RescaledSigma_IC86_II_III_IV_SplineMPE.C+")
    ROOT.loadSplines_IC86_II_III_IV_SplineMPE()
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/TreeLoader_IC86_II_III_IV_lc_MESEremoved_DOUBLESremoved.C")
    gROOT.ProcessLine("I3Ark psark")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/macro_llh/IC86-IV_TDep/load_ark_ic86_II_III_IV_tdep_lc_MESEremoved_DOUBLESremoved.C(psark, OPT_USEREALDATA,\"SplineMPE\")")

    newllh86II_to_IV=ROOT.NewLlhGausTime()
    newllh86II_to_IV.SetUseEnergy(True)
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
    newllh86PSII_to_IV.SetUseEnergy(True)
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

    gROOT.ProcessLine("MultiArk mark")
    ROOT.mark.AddArk(ROOT.MESEark)
    ROOT.mark.AddArk(ROOT.psark)          

    maf=ROOT.MultiGaussAnalysisFn()
    maf.AddAnalysisFn(newllh86II_to_IV)
    maf.AddAnalysisFn(newllh86PSII_to_IV)
    maf.SetTimeBounds(ROOT.MESEark.tmin,ROOT.MESEark.tmax )
    maf.SetTimeBounds(ROOT.psark.tmin,ROOT.psark.tmax )

    testSearch=ROOT.EquatorialDeg(0., 0.)
    tmean = (ROOT.MESEark.tmax + ROOT.MESEark.tmin)/2.
    tmean = (ROOT.psark.tmax + ROOT.psark.tmin)/2.
    tPdf = ROOT.GaussianTimePdf(ROOT.psark.tmin, ROOT.psark.tmax, tmean, 1e-5, 1.)
    newllh86II_to_IV.SetTimeBounds(tPdf)
    newllh86PSII_to_IV.SetTimeBounds(tPdf)
    ROOT.mark.SetPointSource(testSearch, ROOT.PowerLawFlux(1,-2), tPdf)
    maf.SetSearchCoord(testSearch)
    
    pt=ROOT.NewLlhGausTime_ParTranslator()
    pt.SetRange(1,4,31) #gamma_min, gamma_max, nBins
    gROOT.ProcessLine("MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark.psData)")
    pt.SetTranslator(ROOT.mas)

    maf.SetParTranslator(pt)

    nAzResetPS=ROOT.psark.lcBkgProb.nbAz
    nAzResetMESE=ROOT.MESEark.lcBkgProb.nbAz
    ROOT.psark.psData.GetSource().SetTimeAzBins( nAzResetPS) 
    ROOT.MESEark.psData.GetSource().SetTimeAzBins( nAzResetMESE) 

    resultThreshold=2.5
    ra_range=[0, 360]
    dec_range=[-85, 85]
    coarseBinsPerDeg = 2
    fineBinsPerDeg = 10

    nBinsCoarseRa = int((ra_range[1]-ra_range[0])*coarseBinsPerDeg)
    nBinsCoarseDec = int((dec_range[1]-dec_range[0])*coarseBinsPerDeg)

    nBinsFineRa = int((ra_range[1]-ra_range[0])*fineBinsPerDeg)
    nBinsFineDec = int((dec_range[1]-dec_range[0])*fineBinsPerDeg)

    hAllSkyCoarse=TH2D("hAllSkyCoarse","hAllSkyCoarse",nBinsCoarseRa,ra_range[0],ra_range[1],nBinsCoarseDec,dec_range[0],dec_range[1])
    hAllSkyFine=TH2D("hAllSkyFine","hAllSkyFine",nBinsFineRa,ra_range[0],ra_range[1],nBinsFineDec,dec_range[0],dec_range[1])
    print "Coarse Grid: RA Bins: " , nBinsCoarseRa,"   Dec Bins: " , nBinsCoarseDec

    print  "Fine Grid: RA Bins: " , nBinsFineRa, "   Dec Bins: " , nBinsFineDec 

    print  "Result Threshold for fine-grid follow-up: " , resultThreshold
    

    resultMax = -9999
    decDegMax = 0
    raDegMax = 0
    nsBest = 0
    gammaBest = 0
    meanBest = 0
    sigmaBest = 0
    llhBest = 0

    canAllSkyCoarse =  TCanvas("canAllSkyCoarse","canAllSkyCoarse",800,400)
    canAllSkyFine = TCanvas("canAllSkyFine","canAllSkyFine",20,40,800,400)
    canAllSkyCoarse.cd()
    canAllSkyCoarse.SetLogz()
    for iDec in range(1,nBinsCoarseDec):
        decDeg = hAllSkyCoarse.GetYaxis().GetBinCenter(iDec)
        for iRa in range(1,nBinsCoarseRa):
            raDeg = hAllSkyCoarse.GetXaxis().GetBinCenter(iRa)
            print "point dec ra " , decDeg, raDeg
            testSearch=ROOT.EquatorialDeg(raDeg, decDeg)
            maf.SetSearchCoord(testSearch)
            maf.MaximizeLlh()

            result = -np.log10(maf.GetEstProb())
            hAllSkyCoarse.SetBinContent(iRa, iDec, result)  
            hAllSkyCoarse.Draw("colz")
            canAllSkyCoarse.Update()

            if result > 1:
                print result, " (" , raDeg , ":" , decDeg , ") " , maf.Get_logLambdaBest() , " " , maf.GetPar(2) , " " , maf.GetPar(3)
            if result>resultMax:
                resultMax = result
                decDegMax = decDeg
                raDegMax = raDeg
                nsBest = maf.GetPar(0)
                gammaBest = maf.GetPar(1)
                meanBest = maf.GetPar(2)
                sigmaBest = maf.GetPar(3)
                llhBest = maf.Get_logLambdaBest()
    
    print "Coarse Grid Hottest Spot:"
    print "   Ra: " , raDegMax , " , Dec: " , decDegMax 
    print "   logLambda =  " , llhBest , "      "
    print "   -log10(p) =  " , resultMax 
    print "          ns = " , nsBest 
    print "       gamma = " , gammaBest 
    canAllSkyFine.cd()
    canAllSkyFine.SetLogz()

    resultMax = -9999.
    for iDec in range(1,nBinsFineDec):
        decDeg = hAllSkyFine.GetYaxis().GetBinCenter(iDec)
        for iRa in range(1,nBinsFineRa):
            raDeg = hAllSkyFine.GetXaxis().GetBinCenter(iRa)
            print "point dec ra " , decDeg, raDeg
            coarseResult = hAllSkyCoarse.GetBinContent( hAllSkyCoarse.FindBin(raDeg,decDeg) )
            
            if coarseResult > resultThreshold:
                testSearch=ROOT.EquatorialDeg(raDeg, decDeg)
                maf.SetSearchCoord(testSearch)
                maf.MaximizeLlh()

                result = -np.log10(maf.GetEstProb())
                resultLLH = 2*maf.GetTestStatistic()
                if result>resultMax:
                    resultMax = result
                    decDegMax = decDeg
                    raDegMax = raDeg
                    nsBest = maf.GetPar(0)
                    gammaBest = maf.GetPar(1)
                    meanBest = maf.GetPar(2)
                    sigmaBest = maf.GetPar(3)
                    llhBest = maf.Get_logLambdaBest()
            else:
                result = coarseResult
            hAllSkyFine.SetBinContent(iRa, iDec, result)
            hAllSkyFine.Draw("colz")
            canAllSkyFine.Update()

    print "Fine Grid Hottest Spot:"
    print "   Ra: " , raDegMax , " , Dec: " , decDegMax 
    print "   logLambda =  " , llhBest , "      "
    print "   -log10(p) =  " , resultMax 
    print "          ns = " , nsBest 
    print "       gamma = " , gammaBest 
    print "  meanBest   = " , meanBest 
    print "  sigmaBest  = " , pow(10,sigmaBest) 

    fileOutput = TFile("IC86II-IV_AllSkyunblinded_ndof3_MESE_doubles_removed.root", "RECREATE")
    hAllSkyCoarse.Write("hAllSkyCoarse")
    hAllSkyFine.Write("hAllSkyFine")
    fileOutput.Close()
    raw_input("Press Enter to continue...")
