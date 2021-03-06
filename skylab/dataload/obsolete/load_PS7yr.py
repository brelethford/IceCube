# -*-coding:utf8-*-
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import tables
import os
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite, astro, neutrinoflux, NewNuFlux
r"""Load different years of data in and load it correctly to the LLH classes."""

#Livetimes (in days)
livetime_IC40  = 375.539
livetime_IC59  = 348.138
livetime_IC79  = 315.506
livetime_IC79_sirin  = 316.15
livetime_IC86I = 332.61
livetime_IC86II = 330.38
livetime_IC86III = 359.95
livetime_IC86IV = 367.21
livetime3yr = livetime_IC86II+livetime_IC86III+livetime_IC86IV
#getting the sample from coenders, so using his livetime too
livetimeMESE3yr = 988.54
livetimeMESEfollowup = 358.402+368.381

hem = np.sin(np.radians(-5.))

#Loading Zone#
projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/AGN_Core_Sample/'
filename_pickle = datafolder+'pickle/'

#Weighting
honda = NewNuFlux.makeFlux('honda2006')
honda.knee_reweighting_model = 'gaisserH3a_elbert' 
flux = honda.getFlux

#For Skylab I only need OneWeight (I think?) I'll keep the function for posterity, but let's just keep oneweight:

    ### Pull Correction ###
#A polynomial exists to make sure dpsi and sigma have the same median. This is different for each year.
#The polynomial was done with incorrect assumptions from IC40 - IC86I- to fix this, we multiply by the following factor:

#Now these are all handled courtesy of their respective prepdata folder
#pullCorrect = 1./1.1774
#
### for IC40:
##https://wiki.icecube.wisc.edu/index.php/IC-40_PS_Cut_Variables
#def RescaledSigma_IC40_SplineMPE( Sigma, Energy):
#  x = np.log10(Energy);
#  return Sigma  * pullCorrect * (1.189207*(4.974823 -1.967809*x +0.2706778*pow(x,2.)))
#  
### for IC59:
##https://wiki.icecube.wisc.edu/index.php/IC59_PS_Cut_Variables#Paraboloid_Sigma_.28corrected.29
#def RescaledSigma_IC59_SplineMPE( Sigma, Energy):
# x = np.log10(Energy);
# return Sigma * pullCorrect * (31.9051 - 24.5591*pow(x,1.)+7.19661*pow(x,2.)-0.908188*pow(x,3.)+0.043108*pow(x,4.))
#
##for IC79:
##https://wiki.icecube.wisc.edu/index.php/IC79_Point_Source_Analysis/Reconstruction
#
#def RescaledSigma_IC79_SplineMPE(Sigma, nch):
#    x = np.log10(nch)
#    return Sigma * pullCorrect * (5.620-7.448*x+3.519*pow(x, 2.0)+0.252*pow(x, 3.0) -0.509*pow(x, 4.0)+0.095*pow(x, 5.0))
#
#    #Currently, the plan is to use a dpsi/sigma of 1, because the function which causes a median of 1 is already found. But, we may need to change this later. The function comes from , and is: https://wiki.icecube.wisc.edu/index.php/IC86_I_Point_Source_Analysis/Analysis_Variables#Paraboloid_Pull_Tests_and_Performance
#    #I'll keep this part separate from loading in data, because it's computationally inexpensive and might need to change.
#
#def RescaledSigma_IC86_SplineMPE( Sigma, Energy):
# x = np.log10(Energy);
# return Sigma * pullCorrect * (79.0 - 86.7 * pow(x,1.) + 38.45 * pow(x,2.) - 8.673 * pow(x,3.) + 1.056 * pow(x,4.) - 0.0658 * pow(x,5.) + 0.00165 * pow(x,6.))

#IC86II-IV are taken from epinat's pull corrected files - I see no reason to redo this at the moment, though may choose to later.
#MESE files are taken from stefan's files - they match the number of events used, and I'll use these for now.

# skylab (not importing right now)
from mstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH
from mstacking.ps_model import  ClassicLLH, EnergyLLH, PowerLawLLH

#The following exp and mc events have been split into dictionaries in accordance with their respective  prepdata folder

##### IC86I ########

##### IC79 #######


##### IC59 #######

##### IC40 #######

##Mike Sez: try 40 bins equal in sindec, energy pdf 40 bins sindec, 30 in log10 energy, range is from 1 - 10 (linearly spaced, logarithm inside))
def init40(energy=False, stacking=True, mode='all', nopull=False, **kwargs):
    if nopull:
      arr_exp = cache.load(filename_pickle+"nopull/sirin_IC40/exp.pickle")
      arr_mc = cache.load(filename_pickle+"nopull/sirin_IC40/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"IC40/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC40/mc.pickle")

    dec_bins = np.unique(np.concatenate([
                           np.linspace(-1., -0.25, 10 + 1),
                           np.linspace(-0.25, 0.0, 5 + 1),
                           np.linspace(0.0, 1., 10 + 1),
                           ]))    

    dec_bins_logE = np.unique(np.concatenate([
                           np.linspace(-1., -0.25, 10 + 1),
                           np.linspace(-0.25, 0.0, 10 + 1),
                           np.linspace(0.0, 1., 10 + 1),
                           ]))     
    #These binnings are done, year specifically, in load.py from stefan.
    energy_bins = [np.linspace(2., 9., 75 + 1), dec_bins_logE]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
        
    else:

        llh_model = ClassicLLH(sinDec_bins=dec_bins,
                               sinDec_range=[-1., 1.])
    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC40, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC40, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init59(energy=False, mode='all', stacking=True, nopull=False, **kwargs):
    if nopull:
      arr_exp = cache.load(filename_pickle+"nopull/sirin_IC59/exp.pickle")
      arr_mc = cache.load(filename_pickle+"nopull/sirin_IC59/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"IC59/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC59/mc.pickle")

    dec_bins = np.unique(np.concatenate([
                           np.linspace(-1., -0.95, 2 + 1),
                           np.linspace(-0.95, -0.25, 25 + 1),
                           np.linspace(-0.25, 0.05, 15 + 1),
                           np.linspace(0.05, 1., 10 + 1),
                           ]))    

    dec_bins_logE = np.unique(np.concatenate([
                           np.linspace(-1., -0.05, 20 + 1),
                           np.linspace(0.05, 1., 10 + 1),
                           ]))     

    energy_bins = [np.linspace(2., 9.5, 67 + 1), dec_bins_logE]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
        
        #llh_model = EnergyLLH(sinDec_bins=dec_bins(50),#min(50, Nexp // 50),
        #                      sinDec_range=[-1., 1.],
        #                      logE_bins=min(50, Nexp // 50), #energybins(50)
        #                      logE_range=[[0.9 * min(arr_exp["logE"].min(),
        #                                               arr_mc["logE"].min()),
        #                                     1.1 * max(arr_exp["logE"].max(),
        #                                               arr_mc["logE"].max())],
        #                                     [-1., 1.]]) 
    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC59, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC59, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init79_sirin(energy=False, decorr=False, stacking=True, nopull=False, mode='all', **kwargs):
    if decorr:
      arr_exp = cache.load(filename_pickle+"sirin_IC79/noMESE/exp.pickle")
      arr_mc = cache.load(filename_pickle+"sirin_IC79/noMESE/mc.pickle")
    elif nopull:
      arr_exp = cache.load(filename_pickle+"nopull/sirin_IC79/exp.pickle")
      arr_mc = cache.load(filename_pickle+"nopull/sirin_IC79/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"sirin_IC79/exp.pickle")
      arr_mc = cache.load(filename_pickle+"sirin_IC79/mc.pickle")

    dec_bins = np.unique(np.concatenate([
                           np.linspace(-1., -0.75, 10 + 1),
                           np.linspace(-0.75, 0.0, 15 + 1),
                           np.linspace(0.0, 1., 20 + 1),
                           ]))    

    energy_bins = [np.linspace(2., 9., 67 + 1), dec_bins]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
        
        #llh_model = EnergyLLH(sinDec_bins=dec_bins(50),#min(50, Nexp // 50),
        #                      sinDec_range=[-1., 1.],
        #                      logE_bins=min(50, Nexp // 50), #energybins(50)
        #                      logE_range=[[0.9 * min(arr_exp["logE"].min(),
        #                                               arr_mc["logE"].min()),
        #                                     1.1 * max(arr_exp["logE"].max(),
        #                                               arr_mc["logE"].max())],
        #                                     [-1., 1.]]) 
    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC79_sirin, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC79_sirin, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init79(energy=False, stacking=True, decorr=False, mode='all', **kwargs):
    if decorr:
      arr_exp = cache.load(filename_pickle+"IC79/noMESE/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC79/noMESE/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"IC79/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC79/mc.pickle")

    dec_bins = np.unique(np.concatenate([
                           np.linspace(-1., -0.75, 10 + 1),
                           np.linspace(-0.75, 0.0, 15 + 1),
                           np.linspace(0.0, 1., 20 + 1),
                           ]))    

    energy_bins = [np.linspace(2., 9., 67 + 1), dec_bins]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
        
        #llh_model = EnergyLLH(sinDec_bins=dec_bins(50),#min(50, Nexp // 50),
        #                      sinDec_range=[-1., 1.],
        #                      logE_bins=min(50, Nexp // 50), #energybins(50)
        #                      logE_range=[[0.9 * min(arr_exp["logE"].min(),
        #                                               arr_mc["logE"].min()),
        #                                     1.1 * max(arr_exp["logE"].max(),
        #                                               arr_mc["logE"].max())],
        #                                     [-1., 1.]]) 
    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC79, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC79, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    return llh

def init86I(energy=False, decorr=False, stacking=True, nopull=False, mode='all', **kwargs):
    if decorr:
      arr_exp = cache.load(filename_pickle+"IC86I/noMESE/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC86I/noMESE/mc.pickle")
    elif nopull:
      arr_exp = cache.load(filename_pickle+"nopull/sirin_IC86I/exp.pickle")
      arr_mc = cache.load(filename_pickle+"nopull/sirin_IC86I/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"IC86I/exp.pickle")
      arr_mc = cache.load(filename_pickle+"IC86I/mc.pickle")
    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.2, 10 + 1),
                         np.linspace(-0.2, hem, 4 + 1),
                         np.linspace(hem, 0.2, 5 + 1),
                         np.linspace(0.2, 1., 10),
                         ]))
    #check this later, BEN!
    energy_bins = [np.linspace(1., 10., 67 + 1), dec_bins]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
        
        #llh_model = EnergyLLH(sinDec_bins=dec_bins(50),#min(50, Nexp // 50),
        #                      sinDec_range=[-1., 1.],
        #                      logE_bins=min(50, Nexp // 50), #energybins(50)
        #                      logE_range=[[0.9 * min(arr_exp["logE"].min(),
        #                                               arr_mc["logE"].min()),
        #                                     1.1 * max(arr_exp["logE"].max(),
        #                                               arr_mc["logE"].max())],
        #                                     [-1., 1.]]) 
    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC86I, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC86I, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init3yr(energy=False, stacking=True, decorr = False, mode='all', **kwargs):
    if decorr:
      arr_exp = cache.load(filename_pickle+"epinat_3yr/noMESE/exp.pickle")
      arr_mc = cache.load(filename_pickle+"epinat_3yr/noMESE/mc.pickle")
    else:
      arr_exp = cache.load(filename_pickle+"epinat_3yr/exp.pickle")
      arr_mc = cache.load(filename_pickle+"epinat_3yr/mc.pickle")
    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.93, 4 + 1),
                         np.linspace(-0.93, -0.3, 10 + 1),
                         np.linspace(-0.3, 0.05, 9 + 1),
                         np.linspace(0.05, 1., 18 + 1),
                         ]))

    energy_bins = [np.linspace(1., 9.5, 50 + 1), dec_bins]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)

    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime3yr, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetime3yr, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def initMESE3yr(energy=False, stacking=True, mode='all', **kwargs):
    arr_exp = cache.load(filename_pickle+"MESE/exp3yr.pickle")
    arr_mc = cache.load(filename_pickle+"MESE/mc.pickle")

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.92, 3 + 1),
                         np.linspace(-0.92, hem, 10 + 1),
                         ]))

    energy_bins = [np.linspace(2., 8.5, 40 + 1),
                   np.linspace(-1., hem, 4 + 1),
                   ]

    arr_mc = arr_mc[arr_mc["logE"] > 1.]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)

    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetimeMESE3yr, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetimeMESE3yr, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def initMESEfollowup(energy=False, stacking=True, mode='all', **kwargs):
    arr_exp = cache.load(filename_pickle+"MESE/expfollowup.pickle")
    arr_mc = cache.load(filename_pickle+"MESE/mc.pickle")

    dec_bins = np.unique(np.concatenate([np.linspace(-1., -0.93, 4 + 1),
                                        np.linspace(-0.93, hem, 12 + 1),
                                        ]))

    energy_bins = [np.linspace(2., 8.5, 67 + 1), np.linspace(-1., hem, 4 + 1)]

    arr_mc = arr_mc[arr_mc["logE"] > 1.]


    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)

    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    #Why a different sindec binning? Curious...
    if stacking:
      llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetimeMESEfollowup, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)
    else:
      llh = PointSourceLLH(arr_exp, arr_mc, livetime=livetimeMESEfollowup, llh_model=llh_model,
                         mode=mode, 
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh


def multi_init(n, stacking=True, **kwargs):
    #energy = kwargs.pop("energy", False)

    #Now it looks like we don't need llhmodel ahead of time? Maybe because we can use init to get it from each sample...
    if stacking:
      llh = StackingMultiPointSourceLLH(
                              seed=np.random.randint(2**32),
                              **kwargs)
    else:
      llh = MultiPointSourceLLH(
                              seed=np.random.randint(2**32),
                              **kwargs)
    
    for i in xrange(len(n)):
        llh.add_sample(str(i), n[i])
        print ("Adding sample: {}".format(str(n[i])))
    return llh

def initMC_40():
  mc = cache.load(filename_pickle+"IC40/mc.pickle")
  return mc

def initMC_59():
  mc = cache.load(filename_pickle+"IC59/mc.pickle")
  return mc

def initMC_79(decorr=False):
  if decorr:
    mc = cache.load(filename_pickle+"IC79/noMESE/mc.pickle")
  else:
    mc = cache.load(filename_pickle+"IC79/mc.pickle")
  return mc

def initMC_86I(decorr=False):
  if decorr:
    mc = cache.load(filename_pickle+"IC86I/noMESE/mc.pickle")
  else: 
    mc = cache.load(filename_pickle+"IC86I/mc.pickle")
  return mc

def initMC_3yr(decorr=False):
  if decorr:
    mc = cache.load(filename_pickle+"epinat_3yr/noMESE/mc.pickle")
  else: 
    mc = cache.load(filename_pickle+"epinat_3yr/mc.pickle")
  return mc

def loadMC_3yr_no40():
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"IC79/mc.pickle")
  mc86I = cache.load(filename_pickle+"IC86I/mc.pickle")
  MC = {0:mc59, 1:mc79, 2: mc86I}
  return MC

def loadMC_3yr():
  mc40 = cache.load(filename_pickle+"IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"IC79/mc.pickle")
  MC = {0:mc40, 1:mc59, 2:mc79}
  return MC

def loadMC_3yr_nopull():
  mc40 = cache.load(filename_pickle+"nopull/sirin_IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"nopull/sirin_IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"nopull/sirin_IC79/mc.pickle")
  MC = {0:mc40, 1:mc59, 2:mc79}
  return MC

def loadMC_4yr_nopull():
  mc40 = cache.load(filename_pickle+"nopull/sirin_IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"nopull/sirin_IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"nopull/sirin_IC79/mc.pickle")
  mc86I = cache.load(filename_pickle+"nopull/sirin_IC86I/mc.pickle")
  MC = {0:mc40, 1:mc59, 2:mc79, 3: mc86I}
  return MC

def loadMC_4yr():
  mc40 = cache.load(filename_pickle+"IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"IC79/mc.pickle")
  mc86I = cache.load(filename_pickle+"IC86I/mc.pickle")

  MC = {0:mc40, 1:mc59, 2:mc79, 3: mc86I}
  return MC

def loadMC7yr(sirin=False):
  mc40 = cache.load(filename_pickle+"IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  if sirin:
    mc79 = cache.load(filename_pickle+"sirin_IC79/mc.pickle")
  else:
    mc79 = cache.load(filename_pickle+"IC79/mc.pickle")
  mc86I = cache.load(filename_pickle+"IC86I/mc.pickle")
  mc3yr = cache.load(filename_pickle+"epinat_3yr/mc.pickle")

  MC = {0:mc40, 1:mc59, 2:mc79, 3: mc86I, 4:mc3yr}
  return MC

def loadMC7yr_mese(sirin=False):
  print ("Loading decorrelated 7yr data plus MESE...")
  mc40 = cache.load(filename_pickle+"IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  if sirin:
    mc79 = cache.load(filename_pickle+"IC79/noMESE/mc.pickle")
  else:
    mc79 = cache.load(filename_pickle+"sirin_IC79/noMESE/mc.pickle")
  mc86I = cache.load(filename_pickle+"IC86I/noMESE/mc.pickle")
  mc3yr = cache.load(filename_pickle+"epinat_3yr/noMESE/mc.pickle")
  mcMESE = cache.load(filename_pickle+"MESE/mc.pickle")
  #The llh object has the livetimes for mcMESE 3yr and followup - I think I just use mcMESE for both llh samples.
  MC = {0:mc40, 1:mc59, 2:mc79, 3: mc86I, 4:mc3yr, 5:mcMESE, 6:mcMESE}
  return MC

def load3yr(energy=False, stacking=True, sirin=False, nopull=False, mode = 'all', **kwargs):
    print ("Energy: " + str(energy))
    print ("Mode: " + str(mode))
    print ("No pull?" + str(nopull))
    if stacking:
      llh = StackingMultiPointSourceLLH()
    else:
      llh=MultiPointSourceLLH()
    print("loading IC40...")
    IC40 = init40(energy=energy,stacking=stacking,mode=mode,nopull=nopull)
    print("loading IC59...")
    IC59 = init59(energy=energy,stacking=stacking,mode=mode, nopull=nopull)
    if sirin:
      print("loading IC79 - sirin's version...")
      IC79 = init79_sirin(energy=energy,stacking=stacking,mode=mode, nopull=nopull)
    else:
      print("loading IC79...")
      #Nopull is not used for kai's sample - the correction happened after his dataset started being used.
      IC79 = init79(energy=energy,stacking=stacking,mode=mode)
    samples = [IC40,IC59,IC79]

    for i in xrange(len(samples)):
        llh.add_sample(str(i), samples[i])
        print ("Adding sample" + str(samples[i]))
    return llh

def load3yr_no40(energy=False, sirin=False, stacking=True, nopull=False, mode = 'all', **kwargs):
    print ("Energy: " + str(energy))
    print ("Mode: " + str(mode))
    print ("No pull?" + str(nopull))
    print ("Stacking?" + str(stacking))
    if stacking:
      llh = StackingMultiPointSourceLLH()
    else:
      llh=MultiPointSourceLLH()
    print("loading IC59...")
    IC59 = init59(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    if sirin:
      print("loading IC79 - sirin's version...")
      IC79 = init79_sirin(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    else:
      print("loading IC79...")
      #Nopull is not used for kai's sample - the correction happened after his dataset started being used.
      IC79 = init79(energy=energy,stacking=stacking,mode=mode)
    print("loading IC86I...")
    IC86I = init86I(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    samples = [IC59,IC79,IC86I]

    for i in xrange(len(samples)):
        llh.add_sample(str(i), samples[i])
        print ("Adding sample" + str(samples[i]))
    return llh

def load4yr(energy=False, sirin=False, stacking=True, nopull=False, mode = 'all', **kwargs):
    print ("Energy: " + str(energy))
    print ("Mode: " + str(mode))
    print ("No pull?" + str(nopull))
    print ("Stacking?" + str(stacking))
    if stacking:
      llh = StackingMultiPointSourceLLH()
    else:
      llh=MultiPointSourceLLH()
    print("loading IC40...")
    IC40 = init40(energy=energy,mode=mode,stacking=stacking,nopull=nopull)
    print("loading IC59...")
    IC59 = init59(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    if sirin:
      print("loading IC79 - sirin's version...")
      IC79 = init79_sirin(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    else:
      print("loading IC79...")
      #Nopull is not used for kai's sample - the correction happened after his dataset started being used.
      IC79 = init79(energy=energy,stacking=stacking,mode=mode)
    print("loading IC86I...")
    IC86I = init86I(energy=energy,mode=mode,stacking=stacking, nopull=nopull)
    samples = [IC40,IC59,IC79,IC86I]

    for i in xrange(len(samples)):
        llh.add_sample(str(i), samples[i])
        print ("Adding sample" + str(samples[i]))
    return llh

def load7yr(energy=False, sirin=False, stacking=True,mode = 'all', **kwargs):
    print ("Energy: " + str(energy))
    print ("Mode: " + str(mode))
    print ("Stacking?" + str(stacking))
    if stacking:
      llh = StackingMultiPointSourceLLH()
    else:
      llh=MultiPointSourceLLH()
    print("loading IC40...")
    IC40 = init40(energy=energy,mode=mode, stacking=stacking)
    print("loading IC59...")
    IC59 = init59(energy=energy,mode=mode, stacking=stacking)
    if sirin:
      print("loading IC79 - sirin's version...")
      IC79 = init79_sirin(energy=energy,mode=mode, stacking=stacking)
    else:
      print("loading IC79...")
      IC79 = init79(energy=energy,mode=mode, stacking=stacking)
    print("loading IC86I...")
    IC86I = init86I(energy=energy,mode=mode, stacking=stacking)
    print("loading IC86II - IC86IV...")
    IC86II_III_IV = init3yr(energy=energy,mode=mode, stacking=stacking)
    samples = [IC40,IC59,IC79,IC86I,IC86II_III_IV]

    for i in xrange(len(samples)):
        llh.add_sample(str(i), samples[i])
        print ("Adding sample" + str(samples[i]))
    return llh

def load7yr_mese(energy=False,sirin=False,stacking=True,mode='all',**kwargs):
    print ("Energy: " + str(energy))
    print ("Mode: " + str(mode))
    print ("Stacking?" + str(stacking))
    if stacking:
      llh = StackingMultiPointSourceLLH()
    else:
      llh=MultiPointSourceLLH()
    print("loading IC40...")
    IC40 = init40(energy=energy,mode=mode,stacking=stacking)
    print("loading IC59...")
    IC59 = init59(energy=energy,mode=mode,stacking=stacking)
    if sirin:
      print("loading IC79 - sirin's version...")
      IC79 = init79_sirin(energy=energy,mode=mode,stacking=stacking,decorr=True)
    else:
      print("loading IC79...")
      IC79 = init79(energy=energy,mode=mode,stacking=stacking,decorr=True)
    print("loading IC86I...")
    IC86I = init86I(energy=energy,mode=mode,stacking=stacking,decorr=True)
    print("loading IC86II - IC86IV...")
    IC86II_III_IV = init3yr(energy=energy,mode=mode,stacking=stacking,decorr=True)
    print("loading MESE 3yr...")
    MESE3yr = initMESE3yr(energy=energy,mode=mode,stacking=stacking)
    print("loading MESE followup...")
    MESEfollowup = initMESEfollowup(energy=energy,mode=mode,stacking=stacking)
    samples = [IC40,IC59,IC79,IC86I,IC86II_III_IV,MESE3yr,MESEfollowup]

    for i in xrange(len(samples)):
        llh.add_sample(str(i), samples[i])
        print ("Adding sample" + str(samples[i]))
    return llh

def multi_init(n, stacking=True, **kwargs):
    #energy = kwargs.pop("energy", False)

    #Now it looks like we don't need llhmodel ahead of time? Maybe because we can use init to get it from each sample...
    if stacking:
      llh = StackingMultiPointSourceLLH(#hemispheres=dict(Full=[-np.inf, np.inf]),
                              seed=np.random.randint(2**32),
                              **kwargs)
    else:
      llh = MultiPointSourceLLH(#hemispheres=dict(Full=[-np.inf, np.inf]),
                              seed=np.random.randint(2**32),
                              **kwargs)
    #Not sure why it's written this way... I'll try to make it work to preserve as much functionality that I don't understand as possible.
    #for i in xrange(n):
    #    llh.add_sample(str(i), init(energy=energy))
    #On second thought I'll make n a list of the samples I want to put in.

    for i in xrange(len(n)):
        llh.add_sample(str(i), n[i])
    return llh
