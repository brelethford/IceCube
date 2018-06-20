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

parser.add_option ('--years', dest = 'years', type = int,
                default = 2, metavar = 'YEARS',
                help = 'Number of years of data')

opts, args = parser.parse_args ()
dec_deg = np.arcsin(opts.dec) * 180./np.pi
years = opts.years
src_ra=[0.0]
src_dec=[np.radians(dec_deg)]

batch = opts.batch
batchsize = opts.batchsize
import data_multi

llh40 = data_multi.init40(energy=True, mode='box')
samples=[llh40]

if years>1:
  llh59 = data_multi.init59(energy=True, mode='box')
  samples.append(llh59)
if years>2:
  llh79 = data_multi.init79(energy=True, mode='box')
  samples.append(llh79)
if years>3:
  llh86 = data_multi.init86I(energy=True, mode='box')
  samples.append(llh86)
print(len(samples))

llhmodel = data_multi.multi_init(samples,energy=True)

bckg_trials_single = PointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

## Background Trials have the following keys:##
##['beta', 'TS_beta', 'beta_err', 'n_inj', 'nsources', 'TS', 'gamma']##
## Let's use a uniform weight (none) first to yield our bckg trials. ##

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/multi_year/{0}/dec{1:+010.5}/'.format(str(years),dec_deg))

# save the output
outfile = out_dir + 'batch_{0:03}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials_single, outfile)



