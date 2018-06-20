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
from ROOT import TH1D,TFile
import argparse

def main(infile_list,outfile_name,tmin,tmax):

    excludedruns=[]

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
    for run in runs.keys():
        nruns+=1
        if run in excludedruns: 
            print "skipping",run
            continue
        # get date
        year,month,day=runs[run]
        print "filling times for run", run, month,day , nruns,"/",len(runs)
        l3files=glob.glob('/data/ana/Muon/level3/exp/%s/%s%s/Run00%s/Level3_IC86.201?_data_Run00%s_Merged.hdf5' % (year,month,day,run,run))
        l3files=glob.glob('/data/ana/Muon/level3/exp/%s/%s%s/Level3_IC86.201?_data_Run00%s_Part*.hd5' % (year,month,day,run))      
        if len(l3files)==0:
            l3files=glob.glob('/data/ana/Muon/level3/exp/%s/%s%s/Run00%s/Level3_IC86.201?_data_Run00%s_Subrun*.hdf5' % (year,month,day,run,run))    
        print l3files
        for item in l3files:
            try:
                f=tables.openFile(item)
                time_start_mjd_day=f.root.I3EventHeader.col('time_start_mjd_day')
                time_start_mjd_sec=f.root.I3EventHeader.col('time_start_mjd_sec')
                time_start_mjd_ns =f.root.I3EventHeader.col('time_start_mjd_ns')
                time = time_start_mjd_day+(time_start_mjd_sec+time_start_mjd_ns*1e-9)/86400.;
                bins = map(lambda x: minhisto1.Fill(x), time)
                mask = map(lambda b: b if minhisto.GetBinContent(b)==0 else -1, bins )
                
                nonz_bins=np.ma.masked_equal(mask,-1)
                b_to_set=set(nonz_bins.compressed())
                
                for b in b_to_set:
                    minhisto.SetBinContent(int(b),1.)
                f.close()
            except Exception, e:
                    print "problem with file ", item, " pr: ",e
    fth=TFile(outfile_name.replace(".txt",".root"),"RECREATE")
    minhisto1.Write()
    minhisto.Write()
    fth.Close()

    outfile = open(outfile_name, 'w')
    for i in range(int(1E7)):
        rndtime=minhisto.GetRandom()
        outfile.write(str("%.11f" % rndtime)+"\n")
    outfile.close()


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a list of random times for time scrambling')  
    parser.add_argument('-i',metavar='input file list',    type=str,  default="", help='the runlists, coma separated, in doublequotes')
    parser.add_argument('-o',metavar='output file',   type=str,  default="", help='HugeListOf ....txt')
    parser.add_argument('-t',metavar='start time',   type=float,  default=-1., help='start time mjd')
    parser.add_argument('-T',metavar='end time',   type=float,  default=-1., help='end time mjd')
    args = parser.parse_args()
    infile=args.__dict__["i"]
    outfile=args.__dict__["o"]
    tmin=args.__dict__["t"]
    tmax=args.__dict__["T"]

    main(infile,outfile,tmin,tmax)
    #python CreateHugeListOf_HESE_TimesNewMethod.py -t 55694.99085730 -T 57160.0410 -i "IC86I_HESE_GoodRunInfo.txt,IC86II_a_HESE_GoodRunInfo.txt,IC86II_b_HESE_GoodRunInfo.txt" -o HugeListOfTimes_IC86_I_II_III_HESE.txt
    
