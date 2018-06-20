# -*-coding:utf8-*-
from __future__ import print_function
import numpy as np
import scipy as sp
import tables
import sys
import os
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite, astro, neutrinoflux, NewNuFlux
r"""Load different years of data in and load it correctly to the LLH classes."""

#Livetimes (in days)
livetime_IC86I = 332.61
livetime_IC79  = 315.506
livetime_IC59  = 348.138
livetime_IC40  = 375.539

#Loading Zone#
projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/AGN_Core_Sample/'
filename_pickle = datafolder+'pickle/'

mhuber_datafolder = filename_pickle + 'mhuber_7yr/'
sys.path.append(mhuber_datafolder)

import load_bstacking

#End Loading Zone#

#Weighting
honda = NewNuFlux.makeFlux('honda2006')
honda.knee_reweighting_model = 'gaisserH3a_elbert' 
flux = honda.getFlux

# skylab (not importing right now)
from bstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH
from bstacking.ps_model import  ClassicLLH, EnergyLLH, PowerLawLLH

#Here, let's reinstitute the randomness, so that I can redo it in specific things.
#np.random.seed(1)

#Nexp40=len(ra40_data)
#Nmc40=len(ra40_sim)

#Nexp59=len(ra59_data)
#Nmc59=len(ra59_sim)

#Nexp79=len(ra79_data)
#Nmc79=len(ra79_sim)

#Nexp86I=len(ra86I_data)
#Nmc86I=len(ra86I_sim)

#The following exp and mc events have been split into dictionaries in accordance with their respective  prepdata folder

##### IC86I ########

##### IC79 #######


##### IC59 #######

##### IC40 #######

##Mike Sez: try 40 bins equal in sindec, energy pdf 40 bins sindec, 30 in log10 energy, range is from 1 - 10 (linearly spaced, logarithm inside))

def init86I(energy=False,  mode='all', **kwargs):
    arr_exp = cache.load(filename_pickle+"IC86I/exp.pickle")
    arr_mc = cache.load(filename_pickle+"IC86I/mc.pickle")
    #This hem... not sure about what it means, but stefan uses it in load.py.
    #Obviously it's for the hemisphere line - but why doesn't he use it for every year?
    hem = np.sin(np.radians(-5.))
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

    #Why a different sindec binning? Curious...
    llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC86I, llh_model=llh_model,
                         mode=mode, #hemisphere=dict(Full=[-np.inf, np.inf]),
                         #nsource=Nexp / 100.,
                         #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                         #               if not energy else (0., Nexp / 2.),
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh


def init79(energy=False,  mode='all', **kwargs):
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

    #Why a different sindec binning? Curious...
    llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC79, llh_model=llh_model,
                         mode=mode, #hemispheres=dict(Full=[-np.inf, np.inf]),
                         #nsource=Nexp / 100.,
                         #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                         #               if not energy else (0., Nexp / 2.),
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init59(energy=False, mode='all', **kwargs):
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

    #Why a different sindec binning? Curious...
    llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC59, llh_model=llh_model,
                         mode=mode, #hemispheres=dict(Full=[-np.inf, np.inf]),
                         #nsource=Nexp / 100.,
                         #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                         #               if not energy else (0., Nexp / 2.),
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh

def init40(energy=False, mode='all', **kwargs):
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

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    #Why a different sindec binning? Curious...
    llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC40, llh_model=llh_model,
                         mode=mode, #hemispheres=dict(Full=[-np.inf, np.inf]),
                         #nsource=Nexp / 100.,
                         #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                         #               if not energy else (0., Nexp / 2.),
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh
##Okay, so now let's figure out how to use... this.
##It looks like it makes an LLH without using weighting or stipulating an event sample... how? is it related to what n is?
##What is 'n'? Looks like it gets used in llh.add_sample as a string object - is this how the llh is found? they add samples separately after defining the llh object first?

def multi_init(n, energy=False, **kwargs):
    #energy = kwargs.pop("energy", False)

    #Now it looks like we don't need llhmodel ahead of time? Maybe because we can use init to get it from each sample...

    llh = StackingMultiPointSourceLLH(#hemispheres=dict(Full=[-np.inf, np.inf]),
                              #nsource=Nexp / 100.,
                              #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                              #               if not energy else (0., Nexp / 2.),
                              seed=np.random.randint(2**32),
                              **kwargs)
    #Not sure why it's written this way... I'll try to make it work to preserve as much functionality that I don't understand as possible.
    #for i in xrange(n):
    #    llh.add_sample(str(i), init(energy=energy))
    #On second thought I'll make n a list of the samples I want to put in.

    for i in xrange(len(n)):
        llh.add_sample(str(i), n[i])

    return llh

def initMC_40():
  mc = cache.load(filename_pickle+"IC40/mc.pickle")
  return mc

def initMC_59():
  mc = cache.load(filename_pickle+"IC59/mc.pickle")
  return mc

def initMC_79():
  mc = cache.load(filename_pickle+"IC79/mc.pickle")
  return mc

def initMC_86I():
  mc = cache.load(filename_pickle+"IC86I/mc.pickle")

  return mc


def initMC_4yr():
  mc40 = cache.load(filename_pickle+"IC40/mc.pickle")
  mc59 = cache.load(filename_pickle+"IC59/mc.pickle")
  mc79 = cache.load(filename_pickle+"IC79/mc.pickle")
  mc86I = cache.load(filename_pickle+"IC86I/mc.pickle")

  MC = {0:mc40, 1:mc59, 2:mc79, 3: mc86I}
  return MC

