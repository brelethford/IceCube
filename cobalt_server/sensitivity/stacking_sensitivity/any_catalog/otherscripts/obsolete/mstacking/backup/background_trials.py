#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from mstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH #What I'm using - uses energy info
from mstacking.ps_injector import PointSourceInjector, StackingPointSourceInjector
from scipy.stats import chi2
import itertools
from mstacking.utils import poisson_weight
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 30000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--years', dest = 'years', type = int,
                default = 4, metavar = 'YEARS',
                help = 'Number of years of data')


opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
years = opts.years

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':flux, 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == '4yr_Starburst':
        flux = params['S60m']
        modelweights = {'flux':flux}
elif catalog == '30youngSNR':
        weights = params['weights']
        modelweights = {'weights':weights}
elif catalog == '2LAC':
        flux = params['flux']
        modelweights = {'flux':flux}

#Now I've defined all the possible modelweights for each catalog (I can add new ones later). Now I'll read in the years of data.
if catalog == '4yr_Starburst':
        import mstarburst_data as datascript
        print ("Importing starburst data")
else:
        import mstacking_data_multi as datascript


##Here I add an argument for the sole purpose of having separate background files differently named. ##

## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##

llh79 = datascript.init79(energy=True,mode='all')
llh86I= datascript.init86I(energy=True,mode='all')
llh59= datascript.init59(energy=True,mode='all')
llh40= datascript.init40(energy=True,mode='all')

# We use however many years of data are necessary.
if years == 1:
  samples = [llh40]
if years == 2:
  samples = [llh40,llh59]
elif years == 3:
  samples = [llh40,llh59,llh79]
elif years == 4:
  samples = [llh40,llh59,llh79,llh86I]

llhmodel = datascript.multi_init(samples,energy=True)


bckg_trials = StackingMultiPointSourceLLH.do_trials(llhmodel,src_ra,src_dec,w_theo=np.array(modelweights['{}'.format(llhweight)]),n_iter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/background_trials/'.format(catalog,str(years),llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





