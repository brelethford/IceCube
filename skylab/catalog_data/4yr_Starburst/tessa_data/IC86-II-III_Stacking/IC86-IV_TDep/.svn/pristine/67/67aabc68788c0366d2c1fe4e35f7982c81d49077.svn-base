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

def pSigma(err1,err2,med,stat):
    if stat==0:
        return sqrt(pow(err1,2)+pow(err2,2))/sqrt(2)
    return med
    

def hem(zen):
    if zen>=np.radians(85.):
        return True
    return False

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
  

def weight(oneW,E,nEv,maxMCE):
    if maxMCE==7:
        nfiles=6089
    elif maxMCE==9:
        nfiles=4829
    else:
        print "unknown sample"
        exit()
    return oneW*pow(E,-2.)/(nEv*nfiles)

def main(argv):

    LOADTREE_ANADIR= "/net/user/tcarver/Data/PS/IC86_2012/";

    tr = TChain("MasterTree")
    tr.Add( LOADTREE_ANADIR+"Merged_11069_upgoing_6089filesIC86.2012.root"  )
    tr.Add( LOADTREE_ANADIR+"Merged_11069_downgoing_6089filesIC86.2012.root"  )
    tr.Add( LOADTREE_ANADIR+"Merged_11070_upgoing_4829filesIC86.2012.root"  )
    tr.Add( LOADTREE_ANADIR+"Merged_11070_downgoing_4829filesIC86.2012.root"  )
    tr.GetEntries()
    
    MCzenith=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.zenith",MCzenith)
    
    MCazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.azimuth",MCazimuth)

    eneMC=ctypes.c_double(0)
    tr.SetBranchAddress("MCPrimary1.energy",eneMC)

    SMPEzenith=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPE.zenith",SMPEzenith)

    SMPEazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPE.azimuth",SMPEazimuth)
    
    ene=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEMuEXDifferential.energy",ene)
    
    err1=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEParaboloidFitParams.err1",err1)

    err2=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEParaboloidFitParams.err2",err2)
    
    parstat=ctypes.c_int(0)
    tr.SetBranchAddress("SplineMPEParaboloidFitParams.status",parstat)
    
    sigma=ctypes.c_double(0)
    tr.SetBranchAddress("SplineMPEBootstrapVectStats.median",sigma)
    
    oneW=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.OneWeight",oneW)
    
    nEv=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.NEvents",nEv)
    
    maxMCE=ctypes.c_double(0)
    tr.SetBranchAddress("I3MCWeightDict.MaxEnergyLog",maxMCE)
    

    h=TH2D("pull","pull",50,1,9,100,0.01,20)
    h.SetStats(0)
    h.GetXaxis().SetTitle("log_{10}(MuEx_energy [GeV])")
    h.GetYaxis().SetTitle("Pull: #Delta#Psi/#sigma_{parab.}")
    g=TGraphAsymmErrors(50)
    medians=[[] for i in range(0,50)]
    mediansW=[[] for i in range(0,50)]

    go=TGraph(50)
    go.SetMarkerStyle(22)
    go.SetMarkerColor(4)
    
    hEneW=TH1D("hEneW","MCPrimary energy; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2} ",100,1,9.5)
    hREneW=TH1D("hREneW","Reco energy SplineMPEMuEXDiff; log_{10}(E [GeV]); #frac{d#Phi}{dE} [GeV^{-1} cm^{-2} s^{-1} sr^{-1} ( #frac{E}{GeV} )^{-2}",100,1,9)

    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not hem(SMPEzenith.value): 
            continue
        #if pSigma(err1.value,err2.value,sigma.value,parstat.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
        hEneW.Fill(log10(eneMC.value),w)
        hREneW.Fill(log10(ene.value),w)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        h.Fill(log10(ene.value),dAngle/pSigma(err1.value,err2.value,sigma.value,parstat.value),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/pSigma(err1.value,err2.value,sigma.value,parstat.value))
            mediansW[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
        #if i>200:
            #break
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
    hEneW.Draw()

    g.Set(npoints)
    g.SetMarkerStyle(20)
    gStyle.SetPalette(1)
    gStyle.SetOptLogz(1)
    spline= TSpline5("spl", TGraph(g))
    filename="SplineforCorr_"
    if hem(0.):
        filename=filename+"down.root"
    else:
        filename=filename+"up.root"
    fout=TFile(filename,"RECREATE")
    spline.Write()
    fout.Close()

    print "loading RescaledSigma_IC86_II_III_IV_SplineMPE.C"
    ROOT.gROOT.ProcessLine(".L RescaledSigma_IC86_II_III_IV_SplineMPE.C+")
    ROOT.loadSplines()
    npoints=0
    for i in range(1,h.GetNbinsX()+1):
        if len(medians[i-1])>3:
            if hem(0.):
                plotZen=0.
            else:
                plotZen=100.
            go.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),ROOT.RescaledSigma_IC86_II_III_IV_SplineMPE(1.,pow(10,h.GetXaxis().GetBinCenter(i)),plotZen))
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

    
    h2=TH2D("pull - corrected","pull - corrected",50,1,9,100,0.01,20)
    h2.GetXaxis().SetTitle("log_{10}(E [GeV])")
    h2.GetYaxis().SetTitle("Pull_{corr}: #Delta#Psi/f_{corr}(#sigma_{parab.})")
    
    
    g2=TGraphAsymmErrors(50)
    medians2=[[] for i in range(0,50)]
    medians2W=[[] for i in range(0,50)]
    g2.SetMarkerStyle(20)
    
    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not hem(SMPEzenith.value): continue
        #if pSigma(err1.value,err2.value,sigma.value,parstat.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        h2.Fill(log10(ene.value),dAngle/(pSigma(err1.value,err2.value,sigma.value,parstat.value)*spline.Eval(log10(ene.value))),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians2[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/(pSigma(err1.value,err2.value,sigma.value,parstat.value)*spline.Eval(log10(ene.value))))
            medians2W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
        #if i>200:
            #break

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
    c2=TCanvas()
    c2.SetLogy()
    c2.SetGridx()
    c2.SetGridy()
    h2.SetStats(0)
    h2.SetMinimum(1)
    h2.Draw("colz")
    g2.Draw("pl same")

    h3=TH2D("pull - corrected xcheck","pull - corrected xcheck",50,1,9,100,0.01,20)
    h3.GetXaxis().SetTitle("log_{10}(E [GeV])")
    h3.GetYaxis().SetTitle("Pull_{corr}: #Delta#Psi/f_{corr}(#sigma_{parab.})")
    g3=TGraphAsymmErrors(50)
    medians3=[[] for i in range(0,50)]
    medians3W=[[] for i in range(0,50)]
    g3.SetMarkerStyle(20)
    SigmaBins=[(0.5,0.6,1),(0.6,0.7,2),(0.7,0.9,4),(0.9,1.3,1),(1.3,1.7,2)]
    h_pdf_corr=[]
    for b in range(len(SigmaBins)):
        h_pdf_corr.append(TH1D("h_pdf_corr"+str(b),"h_pdf_corr"+str(b)+";#Delta#Psi [deg]",50,0,3))
        h_pdf_corr[-1].SetLineColor(SigmaBins[b][2])
        h_pdf_corr[-1].SetLineWidth(2)

    for i in range(0,tr.GetEntries()):
        tr.GetEntry(i)
        if not hem(SMPEzenith.value): continue
        #if pSigma(err1.value,err2.value,sigma.value,parstat.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value,maxMCE.value)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , SMPEzenith.value, SMPEazimuth.value)
        pull=dAngle/ROOT.RescaledSigma_IC86_II_III_IV_SplineMPE(pSigma(err1.value,err2.value,sigma.value,parstat.value),ene.value,np.degrees(SMPEzenith.value))
        s_cor=ROOT.RescaledSigma_IC86_II_III_IV_SplineMPE(pSigma(err1.value,err2.value,sigma.value,parstat.value),ene.value,np.degrees(SMPEzenith.value))
        h3.Fill(log10(ene.value),pull,w)
        for b in range(len(SigmaBins)):
            if SigmaBins[b][0]<=np.degrees(s_cor) and SigmaBins[b][1]>np.degrees(s_cor):
                h_pdf_corr[b].Fill(np.degrees(dAngle),w)
        
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians3[h.GetXaxis().FindBin(log10(ene.value))-1].append(pull)
            medians3W[h.GetXaxis().FindBin(log10(ene.value))-1].append(w)
        #if i>200:
            #break
    
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
    for b in range(3):
        h_pdf_corr[b].Sumw2()
        if b==0:
            h_pdf_corr[b].Draw("HIST")
        else:
            h_pdf_corr[b].Draw("same HIST")
        
        ffunc.append(TF1("ffunc"+str(b),CircGaus,0.001,3.,2))
        ffunc[-1].SetLineColor(SigmaBins[b][2])
        ffunc[-1].SetParameter(0,1.)
        ffunc[-1].SetParameter(1,1.)

        h_pdf_corr[b].Fit(ffunc[-1],"R")
        ffunc[-1].Draw("same")
        leg.AddEntry(h_pdf_corr[b], "range: "+str(SigmaBins[b][0])+" to "+str(SigmaBins[b][1])+" #sigma_{fit}="+'{:.2f}'.format(ffunc[-1].GetParameter(0)),"l")
    leg.Draw()
    
    c31=TCanvas()
    leg1=ROOT.TLegend(0.7,0.7,1.,1.)
    ffunc=[]
    for b in range(3,5):
        h_pdf_corr[b].Sumw2()
        if b==0:
            h_pdf_corr[b].Draw("HIST")
        else:
            h_pdf_corr[b].Draw("same HIST")
        ffunc.append(TF1("ffunc"+str(b),CircGaus,0.001,3.,2))
        ffunc[-1].SetLineColor(SigmaBins[b][2])
        ffunc[-1].SetParameter(0,1.)
        ffunc[-1].SetParameter(1,1.)

        h_pdf_corr[b].Fit(ffunc[-1],"R")
        ffunc[-1].Draw("same")
        leg1.AddEntry(h_pdf_corr[b], "range: "+str(SigmaBins[b][0])+" to "+str(SigmaBins[b][1])+" #sigma_{fit}="+'{:.2f}'.format(ffunc[-1].GetParameter(0)),"l")
    leg1.Draw()
       
    g3.Set(npoints)
    c3=TCanvas()
    c3.SetLogy()
    c3.SetGridx()
    c3.SetGridy()
    h3.SetStats(0)
    h3.SetMinimum(1)
    h3.Draw("colz")
    g3.Draw("pl same")
    c3.SaveAs("RescaledSigma_IC86_II_III_IV_SplineMPE_corrected.png")
    raw_input("Press Enter to continue...")


#~-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])
