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
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'uniform', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--longrun', dest = 'longrun',
                default = False, action = 'store_true',
                help = 'Allows longer runtime for sensitivity calculation (>12 hours)')

opts,args = parser.parse_args ()
Gamma=opts.gamma

llhweight = opts.llhweight

params = cache.load (datafolder + '2LAC/pickle/params.pickle')



###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/2LAC/flux_3yr/background_trials/'

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

## These params contain everything we should need to weight our sources. Remember to convert src_dec to sindec ##

src_ra, src_dec, flux =  params['src_ra'], params['src_dec'], params['flux']


## Now to import my llh model framework. ##
import data_triple

## Like in the background trials, we have to define which llhmodel to use. But for 2LAC let's stick to flux..
llhmodel = data_triple.init(energy=True, weighting = flux)

##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)

##Now, I'll input my injection weighting scheme from the commandline. For 2LAC blazars we weight injections by flux. ##

inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, theo_weight = flux, seed=0) 

sensitivity = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},longrun=opts.longrun,bckg_trials=bckg_trials,eps=0.02,n_iter=250)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/2LAC/flux_3yr/sensitivity/')
#disc_dir =  misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}/{1}_inj/disc/'.format(llhweight, injweight))

# save the output
outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)
#outfile_disc = disc_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

