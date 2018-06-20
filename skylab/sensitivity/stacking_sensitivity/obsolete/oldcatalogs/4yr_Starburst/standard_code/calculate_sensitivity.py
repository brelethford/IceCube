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
                default = 'flux', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = 'flux', metavar = 'INJWEIGHT',
                help = 'Sets the weighting used for signal injection for point source searches.')

parser.add_option ('--longrun', dest = 'longrun', 
                default = False, action = 'store_true',
                help = 'Allows longer runtime for sensitivity calculation (>40 hours)')

parser.add_option ('--n', dest = 'n', type = int, 
                default = 4, metavar = 'N',
                help = 'Number of years of data')


opts,args = parser.parse_args ()
Gamma=opts.gamma

llhweight = opts.llhweight
injweight = opts.injweight
n = opts.n

params = cache.load (datafolder + '4yr_Starburst/pickle/params.pickle')


###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/bstacking_{0}yr/{1}/background_trials/obsolete/'.format(str(n),llhweight)

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

src_ra, src_dec, redshift, flux =  params['ra'], params['dec'], params['z'], params['S60m']

## Time to define my modelweights in a dictionary. ##

modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2))}

## Now to import my llh model framework. ##
import data_multi

## Like in the background trials, we have to define which llhmodel to use.
if llhweight == 'uniform':
  llh79 = data_multi.init79(energy=True)
  llh86I= data_multi.init86I(energy=True)
  llh59= data_multi.init59(energy=True)
  llh40= data_multi.init40(energy=True)

else:
  llh40= data_multi.init40(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh79 = data_multi.init79(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh86I= data_multi.init86I(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh59= data_multi.init59(energy=True, weighting = modelweights['{}'.format(llhweight)])

#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)
if n == 2:
  samples = [llh79,llh86I]
elif n == 3:
  samples = [llh59,llh79,llh86I]
elif n == 4:
  samples = [llh40,llh59,llh79,llh86I]

llhmodel = data_multi.multi_init(samples,energy=True)

##Remember, weighted sensitivity requires src dec in radians.#

##Now, I'll input my injection weighting scheme from the commandline. ##
if injweight == 'uniform':
  inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)
else:
  inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, theo_weight = modelweights['{}'.format(injweight)], seed=0) 

sensitivity = MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},longrun=opts.longrun,bckg_trials=bckg_trials,eps=0.2,n_iter=250)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/{0}yr/{1}/{2}_inj/sensitivity/'.format(str(n),llhweight, injweight))

# save the output
outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

