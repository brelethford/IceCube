#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
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
import itertools
from scipy.signal import convolve2d
from skylab.ps_injector_stack_100 import PointSourceInjector
from skylab.psLLH_stack import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse
## I'll have to add an argument later to do specific weighting schemes beyond this script, but for now I'll just add the ability manually. ##
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")
if os.uname()[1]=='icecreamcone':
    projfolder='/home/relethford/Documents/IceCube_Research/'
    datafolder=projfolder+'Data/AGN_Core_Sample/'
elif 'wisc.edu' in os.uname()[1]:
    projfolder='/home/brelethford/Documents/IceCube_Research/'
    datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'
filename_pickle = projfolder+'Scripts/AGN_Core/sensitivity/pickle/'

if not os.path.exists(filename_pickle):
	os.mkdir(filename_pickle)

def get_SwiftBAT_params():
    filename=open(datafolder + 'SwiftBAT70m/SeyfertPure.csv','r')
    rawdata=filename.read()
    filename.close()
    table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
    data=table[1:]
    for i in range(len(data)):
      for j in range(len(data[i])):
        data[i][j]=data[i][j].split(',')
    catalog=[list(itertools.chain.from_iterable(data[i])) for i in range(len(data))]
    points=[(float(catalog[i][2]),float(catalog[i][3])) for i in range(len(catalog))]
    ra,dec=zip(*points)
    ra=list(ra)
    dec=list(dec)
    redshift=[float(catalog[i][15]) for i in range(len(catalog))]
    lum =[float(catalog[i][16]) for i in range(len(catalog))]
    gamma=[float(catalog[i][11]) for i in range(len(catalog))]
    flux = [float(catalog[i][7]) for i in range(len(catalog))]
    params={'ra':ra,'dec':dec,'redshift':redshift,'gamma':gamma,'flux':flux,'lum':lum}
    return params

params = cache.get (datafolder + 'SwiftBAT70m/pickle/params.pickle', get_SwiftBAT_params)

###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/background_trials/'

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

TSs=file['TS']
TS_beta = np.percentile(TSs, 100.*(1. - beta))
m=np.count_nonzero(np.asarray(TSs) > (TS_beta))
i = len(TSs)
fraction = float(m)/float(i)
beta_err = (np.sqrt(fraction * (1. - fraction) / float(i)) if 0 < beta < 1 else 1.)

##Now we have all the pieces of the original dictionary. Time to glue bckg_trials back in place, in their proper file type.##
bckg_trials = {'n_inj':n_inj,'nsources':np.asarray(nsources), 'TS':np.asarray(TS), 'beta':beta, 'beta_err':beta_err, 'TS_beta':TS_beta, 'gamma':np.asarray(gamma)}

## These params contain everything we should need to weight our sources. Remember to convert src_dec to sindec ##

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], np.sin(params['dec']), params['redshift'], params['gamma'], params['flux'], params['lum']

import datatest_stack_bins

#llh_uniform = datatest_stack.init(energy=True)

llh_flux = datatest_stack_bins.init(energy=True, weighting = flux)

#llh_redshift = datatest_stack.init(energy=True, weighting = np.power(redshift,-2))

#For this part it's unnecessary to populate signal, but remember how to do so down the line...

# init likelihood class
   # llh = datatest_stack.init(energy=True)
   # mc = datatest_stack.MC()
   # extra=datatest_stack.extra()
   # dpsi=extra["dpsi"]
   # print llh

### Sensitivity ###

##This one's gonna work a little differently than the single source sensitivity. First off, we need to calculate the background scrambles ahead of time, with the definition provided in psLLH_stack.py. I'll try to do this all within the same function:##

## Background Trials have the following keys:##
##['beta', 'TS_beta', 'beta_err', 'n_inj', 'nsources', 'TS', 'gamma']##
## Let's use a uniform weight (none) first to yield our bckg trials. ##

##We also need to input an argument, 'trials', which takes a number of mean_mu as an argument (optional).

##I'll seed some common values for testing for now, but take these out in the final thing.##
parser = OptionParser ( usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

opts,args = parser.parse_args ()
Gamma=opts.gamma

##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)
## NOTE: ADD WEIGHTS HERE FOR THE INJECTED EVENTS!! (Naturally only for flux, redshift. (tried with theo_weight, haven't tested it yet)
inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec,theo_weight = flux, seed=0) 

results = PointSourceLLH.weighted_sensitivity(llh_flux,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
print results

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/E0_100/disc/')

# save the output
outfile = out_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(results, outfile)





