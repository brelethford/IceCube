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


datafolder='/data/user/brelethford/Data/'

parser = OptionParser ( usage = '%prog [options]')
parser.add_option ('--sky', dest = 'sky', type = str,
                default = None, metavar = 'SKY',
                help = 'tells which sources to use (and which folders to reference)')

##Here I add an argument for the sole purpose of having separate background files differently named. ##
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 1000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

opts,args = parser.parse_args ()

batch = opts.batch
batchsize = opts.batchsize
sky = opts.sky

params = cache.load (datafolder + 'SwiftBAT70m/pickle/params.pickle')

## These params contain everything we should need to weight our sources. I've got to select the sources I need to use in this example. ##

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

##Choose the one we're going to use...

if sky == 'onenorth':
  src_ra, src_dec, redshift, flux =  [params['ra'][2]], [params['dec'][2]], [params['redshift'][2]], [params['flux'][2]]
elif sky == 'twonorth':
  src_ra, src_dec, redshift, flux =  [params['ra'][2],params['ra'][6]], [params['dec'][2],params['dec'][6]], [params['redshift'][2],params['redshift'][6]], [params['flux'][2],params['flux'][6]]
elif sky == 'onesouth':
  src_ra, src_dec, redshift, flux =  [params['ra'][0]], [params['dec'][0]], [params['redshift'][0]], [params['flux'][0]]
elif sky == 'twosouth':
  src_ra, src_dec, redshift, flux =  [params['ra'][0],params['ra'][7]], [params['dec'][0],params['dec'][7]], [params['redshift'][0],params['redshift'][7]], [params['flux'][0],params['flux'][7]]
elif sky == 'both':
  src_ra, src_dec, redshift, flux =  [params['ra'][0],params['ra'][6]], [params['dec'][0],params['dec'][6]], [params['redshift'][0],params['redshift'][6]], [params['flux'][0],params['flux'][6]]

print ('my sources are at declination(s):')

## There are three modelweights I can use, so lets put them in a dictionary for easy access. ##

modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2))}

import data


## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##
llhmodel = data.init(energy=True, weighting = modelweights['flux'])

bckg_trials = PointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

print (bckg_trials)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth_one_bckg/{}/background_trials/'.format(sky))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





