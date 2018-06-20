#!/usr/bin/env python
"""
this script creates the list of times used for timescrambling
it is based on the previous script by Juanan (psLab_RHEL_6.0_amd64/macro_llh/ic79/PrintOutAllEvTimesIC79.C)
it samples 1E7 times from a histogram of the times of the original (real data) events
I take the same original (real data) events as Jake in (/net/user/jfeintzeig/2012/PointSources/scripts/BDTs/CalcTotalLiveTime.py)
"""

import glob, tables
import numpy as np
from sys import exit
import ROOT
TH1D   =ROOT.TH1D
TFile  =ROOT.TFile 
TChain =ROOT.TChain
ROOT.gROOT.SetBatch()
import argparse

def main(infile_list,outfile_name):

    excludedruns=[]
    excludebadfiles=["Level3_IC79_data_Run00117788_Part0000000x.root",
                     "Level3_IC79_data_Run00116606_Part0000000x.root",
                     "Level3_IC79_data_Run00117407_Part0000000x.root",
                     "Level3_IC79_data_Run00117997_Part0000000x.root",
                     "Level3_IC79_data_Run00117059_Part0000007x.root",
                     "Level3_IC79_data_Run00117476_Part0000000x.root",
                     "Level3_IC79_data_Run00117275_Part0000010x.root",
                     "Level3_IC79_data_Run00116085_Part0000000x.root",
                     "Level3_IC79_data_Run00116077_Part0000011x.root",
                     "Level3_IC79_data_Run00116918_Part0000000x.root",
                     "Level3_IC79_data_Run00118012_Part0000000x.root"]

    tmin = 55347.3
    tmax = 55694.1
    nbins = (tmax-tmin)*(24.*15.) #4min bins?
    
    minhisto = TH1D("minhisto","minhisto",int(nbins),tmin,tmax)
    minhisto1= TH1D("minhisto1","minhisto1",int(nbins),tmin,tmax)

    runs={}
    runlist=[]
    for infile in infile_list.split(","):
        runlist.append(file(infile))

    for rr in runlist:
        rr.readline();rr.readline()
        for line in rr.readlines():
            ll=line.split()
            if ll[1]=="1" and ll[2]=="1":
                ll_date=ll[7].split("/")
                runs[ll[0]]=(ll_date[4],ll_date[7][:2],ll_date[7][2:])

    nruns=0
    l3files=[]
    for run in runs.keys():
        nruns+=1
        if run in excludedruns: 
            print "skipping",run
            continue
        print "filling times for run", run, " ", nruns,"/",len(runs)
        # get date
        year,month,day=runs[run]
        l3files=l3files+glob.glob('/data/exp/IceCube/2010/filtered/level3-mu/exp/run_xxxxx?/Level3_IC79_data_Run00%s_Part*.root'% (run))

    ch = TChain("L3CutsTree");

    for item in l3files:
        exclude=False
        for exf in excludebadfiles:
            if exf in item:
                exclude=True
        if exclude: continue
        ch.Add(item)
    print ch.GetEntries()
    #ch.Draw("(MJDay+(MJSec+MJNanoSec*1e-9)/86400.) >>minhisto1("+str(int(nbins))+","+str(tmin)+","+str(tmax)+")" )
    ch.Draw("MJDay >>minhisto1("+str(int(nbins))+","+str(tmin)+","+str(tmax)+")" )
    print ROOT.minhisto1.GetEntries()
    for i in range(ROOT.minhisto1.GetNbinsX()):
        if ROOT.minhisto1.GetBinContent(i)>0.:
            minhisto.SetBinContent(i,1.)
            
    fth=TFile(outfile_name.replace(".txt",".root"),"RECREATE")
    ROOT.minhisto1.Write()
    minhisto.Write()
    fth.Close()

    outfile = open(outfile_name, 'w')
    for i in range(int(1E6)):
        rndtime=minhisto.GetRandom()
        outfile.write(str("%.11f" % rndtime)+"\n")
    outfile.close()


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a list of random times for time scrambling')  
    parser.add_argument('-i',metavar='input file list',    type=str,  default="", help='the runlists, coma separated, in doublequotes')
    parser.add_argument('-o',metavar='output file',   type=str,  default="", help='HugeListOf ....txt')
    args = parser.parse_args()
    infile=args.__dict__["i"]
    outfile=args.__dict__["o"]

    main(infile,outfile)
