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

def main():
    excludedruns=[120218,    120229,    120230,    120231,    120232,    121107,    121108,    121110,    121307,    121704,    121712,    121747,    121812,    121852,    122014,    122015,    122040,    122069,    122078,    122643,    122809,    122961,    123153,    123297,    123304,    123389,    123392,    123395,    123399,    123418,    123434,    123455,    123555,    123558,    123561,    123563,    123580,    123581,    123590,    123593,    123595,    123598,    123627,    123639,    123648,    123651,    123652,    123660,    123678,    123707,    123713,    123730,    123738,    123770,    123790,    123810,    123850,    123851,    123860,    123861,    123864,    123924,    123951,    123952,    123954,    124004,    124005,    124023,    124024,    124042,    124043,    124091,    124092,    124126,    124127,    124157,    124158,    124188,    124189,    124231,    124232,    124254,    124255,    124261,    124272,    124290,    124291,    124292,    124385,    124609,    124906,    124916,    124958,    124959,    125128,    125199,    125255,    125377,    125435,    125630,    125632,    125633,    125635,    125636,    125637,    125641,    125648,    125649,    125685,    125693,    125705,    125706,    125708,    125709,    125821,    125822,    125824,    125825,    125826,    125890,    125891,    125933,    126093]

    tmin = 56062.4207
    tmax = 57160.0410
    nbins = (tmax-tmin)*(24.*15.) #4min bins?
    
    minhisto = TH1D("minhisto","minhisto",int(nbins),tmin,tmax)
    minhisto1= TH1D("minhisto1","minhisto1",int(nbins),tmin,tmax)

    runs={}
    runlist=[]
    runlist.append(file('IC86_2012_GoodRunInfo.txt'))
    runlist.append(file('IC86_2013_GoodRunInfo.txt'))
    runlist.append(file('IC86_2014_GoodRunInfo.txt'))
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
        print "filling times for run", run, " ", nruns,"/",len(runs)
        # get date
        year,month,day=runs[run]
        l3files=glob.glob('/data/ana/Muon/level3/exp/%s/%s%s/Run00%s/Level3_IC86.201?_data_Run00%s_Merged.hdf5' % (year,month,day,run,run))
        for item in l3files:
            print l3files
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
    fth=TFile("HugeListOfICII_to_IV_TimesNewMethod.root","RECREATE")
    minhisto1.Write()
    minhisto.Write()
    fth.Close()

    outfile = open('HugeListOfICII_to_IV_TimesNewMethod.txt', 'w')
    for i in range(int(1E7)):
        rndtime=minhisto.GetRandom()
        outfile.write(str("%.11f" % rndtime)+"\n")
    outfile.close()


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
