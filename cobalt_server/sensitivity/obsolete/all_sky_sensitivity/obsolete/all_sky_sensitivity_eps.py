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
from skylab.psLLH import PointSourceLLH #What I'm using - uses energy info
from skylab.ps_model import ClassicLLH  #calculates llh from spatial information only
from scipy.stats import chi2
import healpy as hp
from scipy.signal import convolve2d

from skylab.ps_injector import PointSourceInjector
from skylab.psLLH import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse
##This script is used to test the effect of varying eps on the all_sky sensitivity.##
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

if os.uname()[1]=='icecreamcone':
    projfolder='/home/relethford/Documents/IceCube_Research/'
    datafolder=projfolder+'Data/AGN_Core_Sample/'
elif 'wisc.edu' in os.uname()[1]:
    projfolder='/home/brelethford/Documents/IceCube_Research/'
    datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'
filename_pickle = projfolder+'Scripts/AGN_Core/sensitivity/pickle/'

if not os.path.exists(filename_pickle):
	os.mkdir(filename_pickle)

import datatest

# init likelihood class
if __name__=="__main__":
    llh = datatest.init(ncpu=4, energy=True)
    if isinstance(llh, MultiPointSourceLLH):
        mc = dict([(key, datatest.MC()) for key in llh._enum.iterkeys()])
    else:
        mc = datatest.MC()
    extra=datatest.extra()
    dpsi=extra["dpsi"]
    print llh
### Sensitivity ###

#For sensitivity: if we want to increase the number of trials for a given mu, increase n_iter.
def sensitivity_skylab(Gamma=2,src_dec=0, n_iter=1000, n_bckg=100000,eps=0.005):
  inj = PointSourceInjector(Gamma, sinDec_bandwidth=.1, src_dec= src_dec,seed=0)
  return PointSourceLLH.weighted_sensitivity(llh,src_ra=0.0,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,mc=mc,n_iter=n_iter,eps=eps)

#We give it dec as a commandline argument.
parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--dec', dest = 'dec', type = float,
		default = 0., metavar = 'DEC',
		help = 'sin of the source declination.')
#parser = argparse.ArgumentParser(description='Produce a sensitivity')
#parser.add_argument('dec', type=float,
#                    help='sin of the source declination.')
#args = parser.parse_args()
#dec_deg = np.arcsin(args.dec)*180./np.pi
opts, args = parser.parse_args ()
dec_deg = np.arcsin(opts.dec) * 180./np.pi
#Now we do the work...
sensitivity_array = sensitivity_skylab(src_dec=opts.dec, eps=0.003)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
## NOTE: THE OUT DIRECTORY CHANGES DEPENDING ON WHICH SUBMITTER I'M TESTING!:
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/results/eps/')

# save the output
outfile = out_dir + 'sensitivity__dec_{0:+010.5}.array'.format (dec_deg)
print 'Saving', outfile, '...'
cache.save(sensitivity_array, outfile)





