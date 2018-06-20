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

params = np.load (datafolder + '30youngSNR/source_params.npy')

## These params contain everything we should need to weight our sources. ##

src_ra, src_dec, weights =  params[0], params[1], params[2]


import data_multi

##Here I add an argument for the sole purpose of having separate background files differently named. ##
parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 1000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')


opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize

## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##

llh79 = data_multi.init79(energy=True, weighting = weights)
llh86I= data_multi.init86I(energy=True, weighting = weights)
llh59= data_multi.init59(energy=True, weighting = weights)
llh40= data_multi.init40(energy=True, weighting = weights)
#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)

samples = [llh40,llh59,llh79,llh86I]

llhmodel = data_multi.multi_init(samples,energy=True)

bckg_trials = MultiPointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/30youngSNR/old_background_trials/')

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





