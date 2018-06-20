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
from scipy.stats import chi2
import healpy as hp
import itertools
from scipy.signal import convolve2d
#Check this out - see if there's a new version specifically for multi / stacking
from bstacking.ps_injector import PointSourceInjector, StackingPointSourceInjector
from bstacking.utils import poisson_weight
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

params = cache.load (datafolder + 'SwiftBAT70m/pickle/params.pickle')


###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}yr/{1}/background_trials/'.format(str(n),llhweight)

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

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

## Time to define my modelweights in a dictionary. ##

modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(np.atleast_1d(src_dec))))}

## Now to import my llh model framework. ##
import new_data_multi

#I'm pretty sure I don't have to define the weights in the llh model this time - I think I can do it during sensitivity calculation. But wait - if that's the case, how do I include for it in the background TS?
llh79 = new_data_multi.init79(energy=True,mode='box')
llh86I= new_data_multi.init86I(energy=True,mode='box')
llh59= new_data_multi.init59(energy=True,mode='box')
llh40= new_data_multi.init40(energy=True,mode='box')

MC = new_data_multi.initMC_4yr()
#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)
if n == 2:
  samples = [llh79,llh86I]
elif n == 3:
  samples = [llh59,llh79,llh86I]
elif n == 4:
  samples = [llh40,llh59,llh79,llh86I]

llhmodel = new_data_multi.multi_init(samples,energy=True)

##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)

##Now, I'll input my injection weighting scheme from the commandline. I'll also need my mc from StackingMultiPointSourceLLH? It doesn't seem to require it, though the script says otherwise... also, why does line 637 indicate that there's a kwarg of 'w_theo'_ but I can't seem to input one? Maybe it's only when I call it with inj.fill? So it'd be a kwarg for weighted_sensitivity? ##
#Actually, doesn't appear to need an mc ahead of time? by using weighted sensitivity with a llhmodel, we get mc information into the inj object. could also be the case for w_theo.
inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)#w_theo = modelweights['{}'.format(injweight)], seed=0) 


#I think I need an mc sample in this? But how do I get one for a multipointsourcellh?
sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj, mc=MC, w_theoMC=modelweights['{}'.format(llhweight)],eps=0.02,n_iter=250)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
#this_dir = os.path.dirname(os.path.abspath (__file__))
#sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}yr/{1}/{2}_inj/sensitivity/'.format(str(n),llhweight, injweight))

# save the output
#outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)

#print 'Saving', outfile_sens, '...'
#cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

