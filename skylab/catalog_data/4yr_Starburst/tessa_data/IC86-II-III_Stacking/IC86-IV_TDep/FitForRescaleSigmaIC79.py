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
  

def weight(oneW,E,nEv):
    nfiles=9981
    return oneW*pow(E,-2.)/(nEv*nfiles)

def main(argv):

    LOADTREE_ANADIR = "/net/user/aguilar/work/IceCube/scripts/IC79/";
    
    tr = TChain("L3CutsTree")
    tr.Add(LOADTREE_ANADIR+ "final.v19.nugen_numu_6308_9981f_RunID_new.root")

    tr.GetEntries()
    
    MCzenith=ctypes.c_double(0)
    tr.SetBranchAddress("MC_Zr",MCzenith)
    
    MCazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("MC_Ar",MCazimuth)

    eneMC=ctypes.c_double(0)
    tr.SetBranchAddress("MC_En",eneMC)

    SMPEzenith=ctypes.c_double(0)
    tr.SetBranchAddress("MPE_FINAL_MUON_Zd",SMPEzenith)

    SMPEazimuth=ctypes.c_double(0)
    tr.SetBranchAddress("MPE_FINAL_MUON_Ad",SMPEazimuth)
    
    ene=ctypes.c_double(0)
    tr.SetBranchAddress("MuE_FINAL_MUON_En",ene)
    
    sigma=ctypes.c_double(0)
    tr.SetBranchAddress("MPE_Pb_FINAL_MUON_SigmaDeg",sigma)
            
    oneW=ctypes.c_double(0)
    tr.SetBranchAddress("OW",oneW)
    
    nEv=ctypes.c_double(0)
    tr.SetBranchAddress("NEvents",nEv)
    

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
        #if not hem(np.radians(SMPEzenith.value)): 
        #    continue
        #if np.radians(sigma.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value)
        hEneW.Fill(log10(eneMC.value),w)
        hREneW.Fill(log10(ene.value),w)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , np.radians(SMPEzenith.value), np.radians(SMPEazimuth.value))
        h.Fill(log10(ene.value),dAngle/np.radians(sigma.value),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/np.radians(sigma.value))
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
    filename="SplineforCorr_IC79.root"
    #if hem(0.):
    #    filename=filename+"down.root"
    #else:
    #    filename=filename+"up.root"
    fout=TFile(filename,"RECREATE")
    spline.Write()
    fout.Close()

    #print "loading mpfSigmaDegRescaledIC79V2.C"
    #ROOT.gROOT.ProcessLine(".L mpfSigmaDegRescaledIC79V2.C+")
    print "loading RescaledSigma_IC79.C"
    ROOT.gROOT.ProcessLine(".L RescaledSigma_IC79.C+")
    try:
        ROOT.loadSplines()
    except:
        pass
    npoints=0
    for i in range(1,h.GetNbinsX()+1):
        if len(medians[i-1])>3:
            #if hem(0.):
                #plotZen=0.
            #else:
                #plotZen=100.
            #go.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),ROOT.mpfSigmaDegRescaledIC79V2(1.,pow(10,h.GetXaxis().GetBinCenter(i))))
            go.SetPoint(npoints,h.GetXaxis().GetBinCenter(i),ROOT.RescaledSigma_IC79(1.,pow(10,h.GetXaxis().GetBinCenter(i))))
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
        #if not hem(np.radians(SMPEzenith.value)): continue
        #if np.radians(sigma.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , np.radians(SMPEzenith.value), np.radians(SMPEazimuth.value))
        h2.Fill(log10(ene.value),dAngle/(np.radians(sigma.value)*spline.Eval(log10(ene.value))),w)
        if h.GetXaxis().FindBin(log10(ene.value))>0 and h.GetXaxis().FindBin(log10(ene.value))<51:
            medians2[h.GetXaxis().FindBin(log10(ene.value))-1].append(dAngle/(np.radians(sigma.value)*spline.Eval(log10(ene.value))))
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
        #if not hem(np.radians(SMPEzenith.value)): continue
        #if np.radians(sigma.value)>np.radians(5.): 
        #    continue
        w=weight(oneW.value,eneMC.value,nEv.value)
        dAngle=SpaceAngleRad(MCzenith.value,MCazimuth.value , np.radians(SMPEzenith.value), np.radians(SMPEazimuth.value))
        #pull=dAngle/ROOT.mpfSigmaDegRescaledIC79V2(np.radians(sigma.value),ene.value)
        #s_cor=ROOT.mpfSigmaDegRescaledIC79V2(np.radians(sigma.value),ene.value)
        pull=dAngle/ROOT.RescaledSigma_IC79(np.radians(sigma.value),ene.value)
        s_cor=ROOT.RescaledSigma_IC79(np.radians(sigma.value),ene.value)
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
    c3.SaveAs("mpfSigmaDegRescaledIC79V2_corrected.png")
    raw_input("Press Enter to continue...")


#~-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])
