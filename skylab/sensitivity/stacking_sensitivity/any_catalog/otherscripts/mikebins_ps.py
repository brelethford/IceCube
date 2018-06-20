#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from bstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH #What I'm using - uses energy info
from bstacking.ps_injector import PointSourceInjector, StackingPointSourceInjector
from scipy.stats import chi2
import itertools
from bstacking.utils import poisson_weight, delta_chi2
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")


Gamma =2.0 
catalog = 'SwiftBAT70m' 
years = 4

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))


src_ra, src_dec =  params['ra'], params['dec']
src_ra = 0.0
src_dec = np.arcsin(0.5)

TSval = 0 

#Now I've defined all the possible modelweights for each catalog (I can add new ones later). Now I'll read in the years of data.
import mikebins_single as datascript

llh40= datascript.init40(energy=True,mode='all')
llh59= datascript.init59(energy=True,mode='all')
llh79 = datascript.init79(energy=True,mode='all')
llh86I= datascript.init86I(energy=True,mode='all')
## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##

#pop out segments of MC at the same time we select our llh samples.
MC = datascript.initMC_4yr()
samples = [llh40,llh59,llh79,llh86I]
limit = np.round(np.pi/2,7)

llhmodel = datascript.multi_init(samples,energy=True)
src_dec = np.clip(np.arcsin(np.linspace(-1,1,21)),-limit,limit)

inj = PointSourceInjector(2.0, sinDec_bandwidth = 0.05, src_dec = np.atleast_1d(src_dec))
flux_array = []
mu_array = []
#inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)
result =  MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=np.atleast_1d(0.),src_dec=np.atleast_1d(src_dec),alpha=.5,beta=.9,inj=inj,mc=MC,TSval=0.,eps=0.04,n_iter=20)
flux_array.append (inj.mu2flux(1)*1000)
mu_array.append(result['n_inj'])


#DONT change the rounding on limit - I chose it specifically because it rounds down.

#fill the injarrays:

for i in range(len(src_dec)):
  print ('for declination = ' + str(src_dec[i]) + ', mu = ' +str(mu_array[i]) + ', flux = ' + str(flux_array[i]))

def get_results():
  results = {'dec':src_dec, 'mu':mu_array, 'flux':flux_array}
  return results
print 'Saving...'

results = cache.get('/data/user/brelethford/Output/stacking_sensitivity/testing/onesource_mikebins.pickle',get_results)


