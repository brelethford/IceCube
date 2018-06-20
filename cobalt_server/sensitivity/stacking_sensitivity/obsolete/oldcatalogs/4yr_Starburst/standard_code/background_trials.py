#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
#import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from skylab.psLLH_stack import PointSourceLLH #What I'm using - uses energy info
from skylab.ps_model_stack import ClassicLLH  #calculates llh from spatial information only
#from scipy.stats import chi2
#import healpy as hp
#import itertools
#from scipy.signal import convolve2d
from skylab.ps_injector_stack import PointSourceInjector
from skylab.psLLH_stack import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse
## I'll have to add an argument later to do specific weighting schemes beyond this script, but for now I'll just add the ability manually. ##
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'


params = cache.load (datafolder + 'starburst/pickle/params.pickle')

## These params contain everything we should need to weight our sources. note, these declinations and RAs are in degrees, not radians or sins of either. ##
#I call S60m 'flux' here in order to keep all the things that reference 'flux' the same.

src_ra, src_dec, redshift, flux =  params['ra'], params['dec'], params['z'], params['S60m']

#We must switch these to radians in order to correctly calculate background scrambles.
src_dec = np.radians(src_dec)
src_ra = np.radians(src_ra)

## There are three modelweights I can use, so lets put them in a dictionary for easy access. ##
###Note: for this check there's no reason to check any weight except flux. I'll keep the options open though just in case I want to look at it later.
modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2))}

import data_multi

##Here I add an argument for the sole purpose of having separate background files differently named. ##
parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 1000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'flux', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--years', dest = 'years', type = int,
                default = 4, metavar = 'YEARS',
                help = 'Number of years of data')


opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
llhweight = opts.llhweight
years = opts.years

## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##

if llhweight == 'uniform':
  llh79 = data_multi.init79(energy=True)
  llh86I= data_multi.init86I(energy=True)
  llh59= data_multi.init59(energy=True)
  llh40= data_multi.init40(energy=True)
else:
  llh79 = data_multi.init79(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh86I= data_multi.init86I(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh59= data_multi.init59(energy=True, weighting = modelweights['{}'.format(llhweight)])
  llh40= data_multi.init40(energy=True, weighting = modelweights['{}'.format(llhweight)])
#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)

if years == 2:
  samples = [llh79,llh86I]
elif years == 3:
  samples = [llh59,llh79,llh86I]
elif years == 4:
  samples = [llh40,llh59,llh79,llh86I]

llhmodel = data_multi.multi_init(samples,energy=True)

bckg_trials = MultiPointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/{0}yr/{1}/background_trials/'.format(len(samples),llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





