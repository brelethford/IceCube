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


params = cache.load (datafolder + '2LAC/pickle/params.pickle')

## These params contain everything we should need to weight our sources. ##

deg_dec, deg_ra, flux =  params['src_dec'], params['src_ra'], params['flux']

#These are the background locations in degrees. We need to convert them into radians in order to use the software.

src_dec= np.radians(deg_dec)
src_ra = np.radians(deg_ra)

## There are three modelweights I can use, so lets put them in a dictionary for easy access. ##

import data_triple

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

##For this check we'll use the gamma weighting scheme for llh method and injection:
llhmodel = data_triple.init(energy=True, weighting = flux)

bckg_trials = PointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/2LAC/flux_3yr/background_trials/')

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





