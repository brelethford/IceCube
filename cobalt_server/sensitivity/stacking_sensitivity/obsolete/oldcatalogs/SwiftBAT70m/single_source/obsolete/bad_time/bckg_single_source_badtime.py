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
##This entire practice script is just to try to do single_source sensitivities with the stacking code, in a vain hope that it'll come out correct.

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'
filename_plots=projfolder+'Plots/AGNCore/Stacking/'
filename_pickle = projfolder+'Scripts/AGN_Core/sensitivity/pickle/'

## I only need one source, the src_dec of which will be determined by the submitter scripts. ##
parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--dec', dest = 'dec', type = float,
                default = 0., metavar = 'DEC',             
                help = 'sin of the source declination.')

parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 1000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

opts, args = parser.parse_args ()
dec_deg = np.arcsin(opts.dec) * 180./np.pi

src_ra=[0.0]
src_dec=[np.radians(opts.dec)]

batch = opts.batch
batchsize = opts.batchsize
import datatest_stack_badtime

llh_uniform = datatest_stack_badtime.init(energy=True)

##Okay, so the following is the part where we need to split this up into parallel processing. I think the pertinant variable to use here is n_iter... let's test a background scramble with n_iter=5 to see how fast it goes. Though, maybe it's max_iter? check with previously pickled results to see which number of bckg trials we got.

bckg_trials_single = PointSourceLLH.background_scrambles(llh_uniform,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

## Background Trials have the following keys:##
##['beta', 'TS_beta', 'beta_err', 'n_inj', 'nsources', 'TS', 'gamma']##
## Let's use a uniform weight (none) first to yield our bckg trials. ##


#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bad_time/dec{0:+010.5}/'.format(dec_deg))

# save the output
outfile = out_dir + 'background_dec_{0:+010.5}_batch_{1:03}.array'.format(dec_deg,batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials_single, outfile)





