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


params = np.load(datafolder + '30youngSNR/source_params.npy')

## These params contain everything we should need to weight our sources. ##

src_ra, src_dec, weights =  params[0], params[1], params[2]
#these are already in terms of radians... do I need to make them sinDec?


import new_data_multi

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

llh79 = new_data_multi.init79(energy=True,mode='box')
llh86I= new_data_multi.init86I(energy=True,mode='box')
llh59= new_data_multi.init59(energy=True,mode='box')
llh40= new_data_multi.init40(energy=True,mode='box')

#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)

samples = [llh40,llh59,llh79,llh86I]

llhmodel = new_data_multi.multi_init(samples,energy=True)

bckg_trials = StackingMultiPointSourceLLH.do_trials(llhmodel,src_ra,src_dec,w_theo=np.array(weights),n_iter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/30youngSNR/background_trials/')

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





