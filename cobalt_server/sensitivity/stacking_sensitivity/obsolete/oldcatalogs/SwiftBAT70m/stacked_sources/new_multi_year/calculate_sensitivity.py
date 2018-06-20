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

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'uniform', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = 'uniform', metavar = 'INJWEIGHT',
                help = 'Sets the weighting used for signal injection for point source searches.')

parser.add_option ('--n', dest = 'n', type = int, 
                default = 4, metavar = 'N',
                help = 'Number of years of data')


opts,args = parser.parse_args ()
Gamma=opts.gamma

llhweight = opts.llhweight
injweight = opts.injweight
n = opts.n

params = cache.load (datafolder + 'SwiftBAT70m/pickle/params.pickle')


###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_{0}yr/{1}/background_trials/'.format(str(n),llhweight)

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

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

## Time to define my modelweights in a dictionary. ##

modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(np.atleast_1d(src_dec))))}

## Now to import my llh model framework. ##
import new_data_multi

#I'm pretty sure I don't have to define the weights in the llh model this time - I think I can do it during sensitivity calculation.
llh79 = new_data_multi.init79(energy=True,mode='all')
llh86I= new_data_multi.init86I(energy=True,mode='all')
llh59= new_data_multi.init59(energy=True,mode='all')
llh40= new_data_multi.init40(energy=True,mode='all')

#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)
if n == 2:
  samples = [llh79,llh86I]
elif n == 3:
  samples = [llh59,llh79,llh86I]
elif n == 4:
  samples = [llh40,llh59,llh79,llh86I]

llhmodel = new_data_multi.multi_init(samples,energy=True)
MC = new_data_multi.initMC_4yr()
##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)


inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.1, src_dec= src_dec, seed=0)
#CURRENT PROBLEM  -never goes to the 'Estimation sens in region' phase. problem of too many events? wrong n_iter variable?
sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,mc=MC,TSval = TSval,w_theoMC=modelweights['{}'.format(llhweight)],w_theo=modelweights['{}'.format(injweight)],w_theo_fit = modelweights['{}'.format(injweight)],eps=0.05,n_iter=100)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_{0}yr/{1}/{2}_inj/sensitivity/'.format(str(n),llhweight, injweight))

# save the output
outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

