#!/usr/bin/env python
import sys,os,getopt
from ROOT import TFile,TH2D,TCanvas,gStyle,TF1,TGraphAsymmErrors,TChain,TH1D,TGraph
from math import log10,pi,sqrt,pow,acos,sin,cos
import numpy as np
import ROOT
import ctypes
TSpline5=ROOT.TSpline5

def CircGaus(X,par):
    x=X[0]
    sigma=par[0]
    norm=par[1]
    return x*norm*np.e**(-x**2/(2*sigma**2))

def SpaceAngleRad(zen1, azi1, zen2, azi2):
    return acos( sin(zen1)*sin(zen2)*cos(azi1-azi2)+cos(zen1)*cos(zen2))

def pSigma(err1,err2):
    return sqrt(pow(err1,2)+pow(err2,2))/sqrt(2)

def MESECut(FirstLoss_XYDist,FirstLoss_ZDist,ene,LineFitzenith,LineFitazimuth,SMPEzenith,SMPEazimuth):
    if not SMPEzenith<np.radians(85.):
        return False
    #if not min(FirstLoss_XYDist, FirstLoss_ZDist) > -81.*log10(ene) + 426.:
    #    return False
    if not acos(sin(LineFitzenith)*sin(SMPEzenith)*cos(LineFitazimuth - SMPEazimuth)+cos(LineFitzenith)*cos(SMPEzenith))<np.radians(10.**1.62):
        return False
    return True

def weight(oneW,E,nEv,maxMCE):
    if maxMCE==7:
        nfiles=3999
    elif maxMCE==9:
        nfiles=4999
    else:
        print "unknown sample",maxMCE
        exit()
    return oneW*pow(E,-2.)/(nEv*nfiles)

def weighted_percentile(percents, data, weights=None):
      percents=100*percents
      if weights is None:
        return np.percentile(data, percents)
      ind=np.argsort(data)
      d=data[ind]
      w=weights[ind]
      p=1.*w.cumsum()/w.sum()*100
      y=np.interp(percents, p, d)
      return y

def main(argv):

    LOADTREE_ANADIR= "/net/user/tcarver/Data/PS/IC86_2012/";

    tr = TChain("MasterTree")
    tr.Add("/net/user/tcarver/Data/PS/MESE/Nugen/MergedNugenMillipede_3999.root")
    tr.GetEntries()
    
    use11070=True
    
    tr1 = TChain("MasterTree")
    tr1.Add("/net/user/tcarver/Data/PS/MESE/Nugen/MESE_Merged_11070_Nfiles_4999.root")
    tr1.GetEntries()
    
    MCzenith=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.zenith",MCzenith)
    
    MCazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.azimuth",MCazimuth)

    eneMC=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.energy",eneMC)
    
    oneW=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.OneWeight",oneW)
    tr1.SetBranchAddress("I3MCWeightDict.OneWeight",oneW)

    nEv=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.NEvents",nEv)
    tr1.SetBranchAddress("I3MCWeightDict.NEvents",nEv)
 
    maxMCE=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.MaxEnergyLog",maxMCE)
    tr1.SetBranchAddress("I3MCWeightDict.MaxEnergyLog",maxMCE)


    tr1.SetBranchAddress("I3MCWeightDict.PrimaryNeutrinoZenith",MCzenith)
    
    tr1.SetBranchAddress("I3MCWeightDict.PrimaryNeutrinoAzimuth",MCazimuth)
    
    tr1.SetBranchAddress("I3MCWeightDict.PrimaryNeutrinoEnergy",eneMC)

    SMPEzenith=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPE.zenith",SMPEzenith)

    SMPEazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPE.azimuth",SMPEazimuth)
    
    tr1.SetBranchAddress("SplineMPE.zenith",SMPEzenith)

    tr1.SetBranchAddress("SplineMPE.azimuth",SMPEazimuth)

    LineFitzenith=ctypes.c_double(0)
    tr.SetBranchAddress("LineFit.zenith",LineFitzenith)

    LineFitazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("LineFit.azimuth",LineFitazimuth)

    ene=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEMuEXDifferential.energy",ene)
    
    tr1.SetBranchAddress("SplineMPEMuEXDifferential.energy",ene)

    err1=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEParaboloidFitParams.err1",err1)

    err2=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEParaboloidFitParams.err2",err2)
    
    tr1.SetBranchAddress("SplineMPEParaboloidFitParams.err1",err1)

    tr1.SetBranchAddress("SplineMPEParaboloidFitParams.err2",err2)
    
    FirstLoss_XYDist=ctypes.c_double(0)
    tr.SetBranchAddress("Millipede_FirstLoss_XYDist.value",FirstLoss_XYDist)
    
    FirstLoss_ZDist=ctypes.c_double(0)
    tr.SetBranchAddress("Millipede_FirstLoss_ZDist.value",FirstLoss_ZDist)

    h=TH2D("pull","pull",50,2,9,200,0.01,50)
    h.SetStats(0)
    h.GetXaxis().SetTitle("log_{10}(MuEx_energy [GeV])")
    h.GetYaxis().SetTitle("Pull: #Delta#Psi/#sigma_{parab.}")
    g=TGraphAsymmErrors(50)
    medians=[[] for i in range(0,50)]
    mediansW=[[] for i in range(0,50)]
    go=TGraph(50)
    go.SetMarkerStyle(22)
    go.SetMarkerColor(4)

    hEneW1=TH1D("hEneW1","MCPrimary energy; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2} ",100,1,9)
    hEneW2=TH1D("hEneW2","MCPrimary energy; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2} ",100,1,9)
    hEneW2.SetLineColor(2)
    hREneW1=TH1D("hREneW1","Reco energy SplineMPEMuEXDiff; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2}",100,1,9)
    hREneW2=TH1D("hREneW2","Reco energy SplineMPEMuEXDiff; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2}",100,1,9)
    hREneW2.SetLineColor(2)
    
    h_dAngle1=TH1D("h_dAngle1","h_dAngle1",100,0,0.2)
    h_Sigma1=TH1D("h_Sigma1","h_Sigma1",100,0,0.2)
    h_dAngle2=TH1D("h_dAngle2","h_dAngle2",100,0,0.2)
    h_dAngle2.SetLineColor(2)
    h_Sigma2=TH1D("h_Sigma2","h_Sigma2",100,0,0.2)
    h_Sigma2.SetLineColor(2)
    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not MESECut(FirstLoss_XYDist.value,FirstLoss_ZDist.value,ene.value,LineFitzenith.value,LineFitazimuth.value,SMPEzenith.value,SMPEazimuth.value): 
            continue
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        if dAngle>np.radians(10.):continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)        
        hEneW1.Fill(log10(eneMC.value),w)
        hREneW1.Fill(log10(ene.value),w)
        h_dAngle1.Fill(dAngle,w)
        h_Sigma1.Fill(pSigma(err1.value,err2.value),w)
        h.Fill(log10(ene.value),dAngle/pSigma(err1.value,err2.value),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/pSigma(err1.value,err2.value))
            mediansW[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
    if use11070:
        for i in range(0,tr1.GetEntries()):
            tr1.GetEntry(i)
            dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
            if dAngle>np.radians(10.):continue
            w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
            hEneW2.Fill(log10(eneMC.value),w)
            hREneW2.Fill(log10(ene.value),w)
            h_dAngle2.Fill(dAngle,w)
            h_Sigma2.Fill(pSigma(err1.value,err2.value),w)
            h.Fill(log10(ene.value),dAngle/pSigma(err1.value,err2.value),w)
            if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
                medians[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/pSigma(err1.value,err2.value))
                mediansW[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
            
    npoints=0
    for i in range(1,h.GetNbinsX()+1):
        if len(medians[i-1])>3:
            g.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),weighted_percentile(0.5,
                                                                               np.asarray(medians[i-1]),
                                                                               np.asarray(mediansW[i-1])))
            g.SetPointError(npoints,0.,0.,weighted_percentile(0.5-0.341,
                                                              np.asarray(medians[i-1]),                                                                               np.asarray(mediansW[i-1])),
                                          weighted_percentile(0.5+0.341,
                                                              np.asarray(medians[i-1]),                                                                               np.asarray(mediansW[i-1])))                          
            npoints+=1
            
    c0=TCanvas()
    hEneW1.Draw()
    if use11070:
        hEneW2.Draw("same")
    c00=TCanvas()
    hREneW1.Draw()
    if use11070:
        hREneW2.Draw("same")
    c01=TCanvas()
    h_dAngle1.Draw()
    if use11070:
        h_dAngle2.Draw("same")
    c02=TCanvas()
    h_Sigma1.Draw()
    if use11070:
        h_Sigma2.Draw("same")
        
    g.Set(npoints)
    g.SetMarkerStyle(20)
    gStyle.SetPalette(1)
    gStyle.SetOptLogz(1)
    spline= TSpline5("spl", TGraph(g))
    fout=TFile("SplineMPEPullCorrectionMESE_out.root","RECREATE")
    spline.Write()
    fout.Close()

    
    ROOT.gROOT.ProcessLine(".L RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE.C+")
    try:
        ROOT.loadSplines()
    except:
        pass
    npoints=0
    for i in range(1,h.GetNbinsX()+1):
        if len(medians[i-1])>3:
            go.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),ROOT.RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(1.,pow(10,h.GetXaxis().GetBinCenter(i))))
            npoints+=1
    go.Set(npoints)

    del medians[:]
    del medians
    del mediansW[:]
    del mediansW
    
    c=TCanvas()
    c.SetLogy()
    c.SetGridx()
    c.SetGridy()
    h.SetMinimum(1)
    h.Draw("colz")
    g.Draw("pl same")
    go.Draw("pl same")
    spline.Draw("lp same")
    h2=TH2D("pull - corrected","pull - corrected",50,2,9,200,0.01,50)
    h2.GetXaxis().SetTitle("log_{10}(E [GeV])")
    h2.GetYaxis().SetTitle("Pull_{corr}: #Delta#Psi/f_{corr}(#sigma_{parab.})")

    g2=TGraphAsymmErrors(50)
    medians2=[[] for i in range(0,50)]
    medians2W=[[] for i in range(0,50)]
    g2.SetMarkerStyle(20)
    
    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not MESECut(FirstLoss_XYDist.value,FirstLoss_ZDist.value,ene.value,LineFitzenith.value,LineFitazimuth.value,SMPEzenith.value,SMPEazimuth.value): 
            continue
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        if dAngle>np.radians(10.):continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
        h2.Fill(log10(ene.value),dAngle/(pSigma(err1.value,err2.value)*spline.Eval(log10(ene.value))),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians2[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/(pSigma(err1.value,err2.value)*spline.Eval(log10(ene.value))))
            medians2W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
    if use11070:
        for i in range(0,tr1.GetEntries()):
            tr1.GetEntry(i)
            w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
            dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
            if dAngle>np.radians(10.):continue
            h2.Fill(log10(ene.value),dAngle/(pSigma(err1.value,err2.value)*spline.Eval(log10(ene.value))),w)
            if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
                medians2[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/(pSigma(err1.value,err2.value)*spline.Eval(log10(ene.value))))
                medians2W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)

    npoints=0
    for i in range(1,h2.GetNbinsX()+1):
        if len(medians2[i-1])>3:
            g2.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),weighted_percentile(0.5,
                                                                               np.asarray(medians2[i-1]),
                                                                               np.asarray(medians2W[i-1])))
            g2.SetPointError(npoints,0.,0.,weighted_percentile(0.5-0.341,
                                                              np.asarray(medians2[i-1]),                                                                               np.asarray(medians2W[i-1])),
                                          weighted_percentile(0.5+0.341,
                                                              np.asarray(medians2[i-1]),                                                                               np.asarray(medians2W[i-1])))                          
            npoints+=1

    del medians2[:]
    del medians2
    del medians2W[:]
    del medians2W
        
    g2.Set(npoints)
    go.Set(npoints)
    c2=TCanvas()
    c2.SetLogy()
    c2.SetGridx()
    c2.SetGridy()
    h2.SetStats(0)
    h2.SetMinimum(1)
    h2.Draw("colz")
    g2.Draw("pl same")

    h3=TH2D("pull - corrected xcheck","pull - corrected xcheck",50,2,9,200,0.01,50)
    h3.GetXaxis().SetTitle("log_{10}(E [GeV])")
    h3.GetYaxis().SetTitle("Pull_{corr}: #Delta#Psi/f_{corr}(#sigma_{parab.})")
    g3=TGraphAsymmErrors(50)
    medians3=[[] for i in range(0,50)]
    medians3W=[[] for i in range(0,50)]
    g3.SetMarkerStyle(20)
    
    SigmaBins=[(0.5,0.9),(0.9,1.3),(1.3,1.7)]
    h_pdf_corr=[]
    for b in range(len(SigmaBins)):
        h_pdf_corr.append(TH1D("h_pdf_corr"+str(b),"h_pdf_corr"+str(b)+";#Delta#Psi [deg]",100,0,3))
        h_pdf_corr[-1].SetLineColor(b+1)

    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not MESECut(FirstLoss_XYDist.value,FirstLoss_ZDist.value,ene.value,LineFitzenith.value,LineFitazimuth.value,SMPEzenith.value,SMPEazimuth.value): continue
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        if dAngle>np.radians(10.):continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
        pull=dAngle/ROOT.RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(pSigma(err1.value,err2.value),ene.value)
        h3.Fill(log10(ene.value),pull,w)
        s_cor=ROOT.RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(pSigma(err1.value,err2.value),ene.value)
        for b in range(len(SigmaBins)):
            if SigmaBins[b][0]<=np.degrees(s_cor) and SigmaBins[b][1]>np.degrees(s_cor):
                h_pdf_corr[b].Fill(np.degrees(dAngle),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians3[h.GetXaxis().FindBin(log10(ene.value))-1].append(pull)
            medians3W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
        #if i>200:
            #break
    if use11070:
        for i in range(0,tr1.GetEntries()):
            tr1.GetEntry(i)
            dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
            if dAngle>np.radians(10.):continue
            w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
            pull=dAngle/ROOT.RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(pSigma(err1.value,err2.value),ene.value)
            h3.Fill(log10(ene.value),pull,w)
            s_cor=ROOT.RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE(pSigma(err1.value,err2.value),ene.value)
            for b in range(len(SigmaBins)):
                if SigmaBins[b][0]<=np.degrees(s_cor) and SigmaBins[b][1]>np.degrees(s_cor):
                    h_pdf_corr[b].Fill(np.degrees(dAngle),w)
                                          
            if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
                medians3[h.GetXaxis().FindBin(log10(ene.value))-1].append(pull)
                medians3W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
        
    npoints=0
    for i in range(1,h2.GetNbinsX()+1):
        if len(medians3[i-1])>3:
            g3.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),weighted_percentile(0.5,
                                                                               np.asarray(medians3[i-1]),
                                                                               np.asarray(medians3W[i-1])))
            g3.SetPointError(npoints,0.,0.,weighted_percentile(0.5-0.341,
                                                              np.asarray(medians3[i-1]),                                                                               np.asarray(medians3W[i-1])),
                                          weighted_percentile(0.5+0.341,
                                                              np.asarray(medians3[i-1]),                                                                               np.asarray(medians3W[i-1])))                                      
            npoints+=1

    c30=TCanvas()
    leg=ROOT.TLegend(0.7,0.7,1.,1.)
    ffunc=[]
    for b in range(len(SigmaBins)):
        if b==0:
            h_pdf_corr[b].Draw()
        else:
            h_pdf_corr[b].Draw("same")
        ffunc.append(TF1("ffunc"+str(b),CircGaus,0.001,3.,2))
        ffunc[-1].SetLineColor(b+1)
        ffunc[-1].SetParameter(0,1.)
        ffunc[-1].SetParameter(1,1.)

        h_pdf_corr[b].Fit(ffunc[-1],"R")
        leg.AddEntry(h_pdf_corr[b], "range: "+str(SigmaBins[b][0])+" to "+str(SigmaBins[b][1])+" #sigma_{fit}="+'{:.2f}'.format(ffunc[-1].GetParameter(0)),"l")
    leg.Draw()
        

    g3.Set(npoints)
    c3=TCanvas()
    c3.SetLogy()
    c3.SetGridx()
    c3.SetGridy()
    h3.SetStats(0)
    h3.SetMinimum(1)
    h3.Draw("colz")
    g3.Draw("pl same")
    c3.SaveAs("RescaledSigma_IC79_86_I_to_IV_SplineMPE_MESE_corrcted.png")
    raw_input("Press Enter to continue...")


#~-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])
