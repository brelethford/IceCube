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
from bstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH #What I'm using - uses energy info
from bstacking.utils import poisson_weight, delta_chi2
from bstacking.ps_injector import PointSourceInjector, StackingPointSourceInjector
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

opts,args = parser.parse_args ()
Gamma=opts.gamma

params = np.load (datafolder + '30youngSNR/source_params.npy')


###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/30youngSNR/background_trials/'

files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
#We only need to get the appropriate chi_fit such that we can get a TSval for weighted sensitivity fcn.

fitfun = delta_chi2

TS=[]

for file in files:
  for item in range(len(file['n_inj'])):
    if file['n_inj'][item] ==0:
      TS.append(file['TS'][item])

TSs=TS
fit = fitfun(TSs, df=2., floc=0., fscale=1.)
#This gives us the median bckg TS distribution, which is useful for calculating sensitivities (will I need alpha = 0.9 for discovery potentials? unclear as of now...)
TSval = np.asscalar(fit.isf(0.5))

src_ra, src_dec, weights =  params[0], params[1], params[2]

## Now to import my llh model framework. ##
import new_data_multi

llh79 = new_data_multi.init79(energy=True,mode='all')
llh86I= new_data_multi.init86I(energy=True,mode='all')
llh59= new_data_multi.init59(energy=True,mode='all')
llh40= new_data_multi.init40(energy=True,mode='all')

samples = [llh40,llh59,llh79,llh86I]

llhmodel = new_data_multi.multi_init(samples,energy=True)
MC = new_data_multi.initMC_4yr()

inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.1, src_dec= src_dec, seed=0)
#CURRENT PROBLEM  -never goes to the 'Estimation sens in region' phase. problem of too many events? wrong n_iter variable?
sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,mc=MC,TSval = TSval,w_theoMC=weights,w_theo=weights,w_theo_fit = weights,eps=0.02,n_iter=250)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/30youngSNR/sensitivity/')

# save the output
outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

