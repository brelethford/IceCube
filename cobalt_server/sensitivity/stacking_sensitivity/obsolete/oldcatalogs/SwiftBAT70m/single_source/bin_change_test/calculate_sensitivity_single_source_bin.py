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
from skylab.ps_injector_stack import PointSourceInjector
from skylab.psLLH_stack import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse
## I'll have to add an argument later to do specific weighting schemes beyond this script, but for now I'll just add the ability manually. ##
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")
projfolder='/home/brelethford/Documents/IceCube_Research/'

## Read in background trials previously logged ##

import datatest_stack_bins

llh_single = datatest_stack_bins.init(energy=True)

mc = datatest_stack_bins.MC()

##I have to load in whichever declination I'm working on at the moment...##

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--dec', dest = 'dec', type = float,
                default = 0., metavar = 'DEC',
                help = 'sin of the source declination.')

opts, args = parser.parse_args ()
dec_deg = np.arcsin(opts.dec) * 180./np.pi

datafolder='/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/dec{0:+010.5}/'.format(dec_deg)

##I'll load the pieces here and then combine them.##

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

## Now, we need to do the injection trials ahead of time, because doing an estimation is taking multiple hours... way too long!

### Sensitivity ###

##This one's gonna work a little differently than the single source sensitivity. First off, we need to calculate the background scrambles ahead of time, with the definition provided in psLLH_stack.py. I'll try to do this all within the same function:##

## Background Trials have the following keys:##
##['beta', 'TS_beta', 'beta_err', 'n_inj', 'nsources', 'TS', 'gamma']##
## Let's use a uniform weight (none) first to yield our bckg trials. ##

##We also need to input an argument, 'trials', which takes a number of mean_mu as an argument (optional).

##I'll seed some common values for testing for now, but take these out in the final thing.##

Gamma=2

##I need to have this 
src_dec=[np.radians(dec_deg)]
src_ra=[0.0]
## NOTE: ADD WEIGHTS HERE FOR THE INJECTED EVENTS!! (Naturally only for flux, redshift.

##For now the weighted sensitivity function only works if there are at least two sources. to push it through for a single source, I'll copy the source location.##

### Mike - This is the step that allows the function weighted_sensitivity to process! ###

#src_dec=[src_dec[0],src_dec[0]]
#src_ra=[src_ra[0],src_ra[0]]

inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec,seed=0) 

results = PointSourceLLH.weighted_sensitivity(llh_single,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250, miniter=2500)# maxtrial=1000)
#Currently have n_iter down from 1000 to reduce estimation time. Also lowering maxtrial from 1000 to 500.
print results

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir (datafolder+'sensitivity/')

# save the output
outfile = out_dir + 'dec{0:+010.5}.array'.format(dec_deg)
print 'Saving', outfile, '...'
cache.save(results, outfile)





