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

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'

## Input my commandline parameters ##

parser = OptionParser ( usage = '%prog [options]')

parser.add_option ('--dec', dest = 'dec', type = float,
                default = 0., metavar = 'DEC',
                help = 'sin of the source declination.')

parser.add_option ('--year', dest = 'year', type = str,
                default = '86', metavar = 'YEAR',
                help = 'Single year of data')

opts,args = parser.parse_args ()
dec_deg = np.arcsin(opts.dec) * 180./np.pi
year = opts.year

##Not sure if this background folder needs to be redone... but it worked for single source so I don't think so.
datafolder='/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/single_year/IC{0}/dec{1:+010.5}/'.format(year,dec_deg)

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

##I'll seed some common values for testing for now, but take these out in the final thing.##

Gamma=2

##I need to have this 
src_dec=[np.radians(dec_deg)]
src_ra=[0.0]


## Now to import my llh model framework. ##
import data_multi

## Like in the background trials, we have to define which llhmodel to use.
##Now, I'll input my injections scheme (same for cut and uncut) and likelihood model (mc varies with injection range). ##
inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec,seed=0)
if year == '40':
  llhmodel = data_multi.init40(energy=True,mode='box')
elif year == '79':
  llhmodel = data_multi.init79(energy=True,mode='box')
elif year == '86':
  llhmodel = data_multi.init86I(energy=True,mode='box')
elif year == '59':
  llhmodel = data_multi.init59(energy=True,mode='box')  

#If I change the injection range I'll redefine the _e_range variable in ps_injector_stack.py.

sensitivity = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.02,n_iter=250)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

sens_dir = datafolder
  

#disc_dir =  misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}/{1}_inj/disc/'.format(llhweight, injweight))

# save the output
outfile_sens = sens_dir + 'sensitivity.array'
#outfile_disc = disc_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

