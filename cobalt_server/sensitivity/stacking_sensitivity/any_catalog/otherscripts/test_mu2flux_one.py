#!/usr/bin/env python
#NOTE: This whole script is simply a way to test the mu2flux function in mhuber's old version of the code.
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from skylab.psLLH_stack import PointSourceLLH #What I'm using - uses energy info
from skylab.ps_model_stack import ClassicLLH  #calculates llh from spatial information only
from scipy.stats import chi2
import healpy as hp
from skylab.ps_injector_stack import PointSourceInjector
from skylab.psLLH_stack import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'


###Now I'll read in the background TS stuff.###
#This is just the most recent 4yr_data I have. the obsolete should go away once I can access condor again and do the background TS in tandem.
datafolder='/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/bstacking_4yr/flux/background_trials/obsolete/'
#I'm using these simply to fill the injectors through the weighted sens function.
files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
n_inj=[]
nsources=[]
TS=[]
beta=(0.5) #For background ts
TS_beta=[] #Calculated from the total TS median after we get all the TS.
beta_err=[]
gamma=[]

for file in files:
  for item in range(len(file['n_inj'])):
    n_inj.append(file['n_inj'][item])
    nsources.append(file['nsources'][item])
    TS.append(file['TS'][item])
    gamma.append(file['gamma'][item])

TSs=TS
TS_beta = np.percentile(TSs, 100.*(1. - beta))
m=np.count_nonzero(np.asarray(TSs) > (TS_beta))
i = len(TSs)
fraction = float(m)/float(i)
beta_err = (np.sqrt(fraction * (1. - fraction) / float(i)) if 0 < beta < 1 else 1.)

##Now we have all the pieces of the original dictionary. Time to glue bckg_trials back in place, in their proper file type.##
bckg_trials = {'n_inj':n_inj,'nsources':np.asarray(nsources), 'TS':np.asarray(TS), 'beta':beta, 'beta_err':beta_err, 'TS_beta':TS_beta, 'gamma':np.asarray(gamma)}



## Now to import my llh model framework. Both the llhmodel and injector have to be the old version, but the data is all the stuff used by the new version - totally up to date.##
import data_multi

## Like in the background trials, we have to define which llhmodel to use.
llh79 = data_multi.init79(energy=True)
llh86I= data_multi.init86I(energy=True)
llh59= data_multi.init59(energy=True)
llh40= data_multi.init40(energy=True)

#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)
samples = [llh40,llh59,llh79,llh86I]

llhmodel = data_multi.multi_init(samples,energy=True)

##Remember, weighted sensitivity requires src dec in radians.#

##Now, I'll input my injection weighting scheme from the commandline. APPARENTLY, I need these as radians, not sin(dec).##
#Also, I need to clip these just before they hit the pole - otherwise I'm gonna get out of bounds errors.
limit = np.round(np.pi/2,7)
#DONT change the rounding on limit - I chose it specifically because it rounds down.
src_dec = np.clip(np.arcsin(np.linspace(-1,1,21)),-limit,limit)

inj_array = [PointSourceInjector(2.0, sinDec_bandwidth = 0.05, src_dec = np.atleast_1d(src_dec[i])) for i in range(len(src_dec))]
flux_array = []
#fill the injarrays:
for i in range(len(src_dec)):
  inj = inj_array[i]
  MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=np.atleast_1d(0.),src_dec=np.atleast_1d(src_dec[i]),alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=10.99,n_iter=1)
  flux_array.append (inj.mu2flux(1)*1000)

for i in range(len(src_dec)):
  print ('for declination = ' + str(src_dec[i]) + ', flux = ' + str(flux_array[i]))

def get_results():
  results = {'dec':src_dec, 'flux':flux_array}
  return results

results = cache.get('/data/user/brelethford/Output/stacking_sensitivity/testing/onesource.pickle',get_results)

