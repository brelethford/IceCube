#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from jstacking.ps_llh import PointSourceLLH, MultiPointSourceLLH
from jstacking.ps_injector import PointSourceInjector
from jstacking.llh_models import ClassicLLH, EnergyLLH
from jstacking.utils import poisson_weight
from scipy.stats import chi2
import itertools
from optparse import OptionParser
import argparse


sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")

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


opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog

datafolder='/data/user/coenders/data/MultiYearPointSource/npz'
catfolder='/data/user/brelethford/stacktest/{}'.format(catalog)

params = cache.load (catfolder + '/pickle/params.pickle')

years = 2 #for now...

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##
src_ra, src_dec =  params['ra'], params['dec']

exp79 = np.load(datafolder + '/IC79_exp.npy')
mc79 = np.load(datafolder + '/IC79_MC.npy')

exp86 = np.load(datafolder + '/IC86_exp.npy')
mc86 = np.load(datafolder + '/IC86_MC.npy')

#need llhmodel
dec_bins = np.unique(np.linspace(-1., 1, 100 + 1))

energy_bins = [np.linspace( 1.,10.,36),dec_bins]

llh_model = EnergyLLH(energy_bins, sinDec_bins=dec_bins)
llh79 = PointSourceLLH(exp79, mc79, livetime=315.506, llh_model=llh_model,
                         seed=0)#np.random.randint(2**32),

llh86 = PointSourceLLH(exp86, mc86, livetime=332.61, llh_model=llh_model,
                         seed=0)#np.random.randint(2**32),
samples = [llh79,llh86]
years = len(samples)

llh=MultiPointSourceLLH()
for i in xrange(len(samples)):
   llh.add_sample(str(i), samples[i])
   print ("Adding sample" + str(samples[i]))

bckg_trials = MultiPointSourceLLH.do_trials(llh,src_ra,src_dec,n_iter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





