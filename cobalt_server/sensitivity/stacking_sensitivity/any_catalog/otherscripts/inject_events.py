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

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--n_inj', dest = 'n_inj', type = int,
                default = 10, metavar = 'N_INJ',
                help = 'number of events to inject')

opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
years = opts.n
n_inj = opts.n_inj

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == '4yr_Starburst':
        flux = np.array(params['S60m'])
        modelweights = {'flux':flux}
elif catalog == '30youngSNR':
        weights = np.array(params['weights'])
        modelweights = {'weights':weights}
elif catalog == '2LAC':
        flux = np.array(params['flux'])
        modelweights = {'flux':flux}
#Rather than use background trials, we'll forgo the estimation - we just want to inject a set number of events and observe.

#Now I've defined all the possible modelweights for each catalog (I can add new ones later). Now I'll read in the years of data.
if catalog == '4yr_Starburst':
    import starburst_data as datascript
else:
    import new_data_multi as datascript

llh40= datascript.init40(energy=True,mode='all')
llh59= datascript.init59(energy=True,mode='all')
llh79 = datascript.init79(energy=True,mode='all')
llh86I= datascript.init86I(energy=True,mode='all')
## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##

#pop out segments of MC at the same time we select our llh samples.
MC = datascript.initMC_4yr()

if years == 1:
  samples = [llh40]
  MC.pop(1,2,3)
elif years == 2:
  samples = [llh40,llh59]
  MC.pop(2,3)
elif years == 3:
  samples = [llh40,llh59,llh79]
  MC.pop(3)
elif years == 4:
  samples = [llh40,llh59,llh79,llh86I]

#For the 2LAC comparison, they don't use year 40 of data.
if catalog == '2LAC':
  MC = datascript.initMC_4yr()
  MC.pop(0)
  samples = [llh59,llh79,llh86I]

llhmodel = datascript.multi_init(samples,energy=True)
inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)
inj.fill(src_dec, MC, llhmodel.livetime, w_theo=np.array(flux))


result = StackingMultiPointSourceLLH.do_trials(llhmodel,src_ra,src_dec,w_theo=np.array(modelweights['{}'.format(llhweight)]),n_iter=1000,mu=inj.sample(src_ra,mea
n_mu=n_inj)) 

#choose an output dir, and make sure it exists


this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/bstacking_{1}yr/{2}/{3}_inj/sensitivity/'.format(catalog,str(years),llhweight,injweight))

# save the output
outfile = sens_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





