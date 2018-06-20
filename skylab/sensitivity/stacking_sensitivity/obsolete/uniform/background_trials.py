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
## I'll have to add an argument later to do specific weighting schemes beyond this script, but for now I'll just add the ability manually. ##
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

def get_SwiftBAT_params():
    filename=open(datafolder + 'SwiftBAT70m/SeyfertPure.csv','r')
    rawdata=filename.read()
    filename.close()
    table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
    data=table[1:]
    for i in range(len(data)):
      for j in range(len(data[i])):
        data[i][j]=data[i][j].split(',')
    catalog=[list(itertools.chain.from_iterable(data[i])) for i in range(len(data))]
    points=[(float(catalog[i][2]),float(catalog[i][3])) for i in range(len(catalog))]
    ra,dec=zip(*points)
    ra=list(ra)
    dec=list(dec)
    redshift=[float(catalog[i][15]) for i in range(len(catalog))]
    lum =[float(catalog[i][16]) for i in range(len(catalog))]
    gamma=[float(catalog[i][11]) for i in range(len(catalog))]
    flux = [float(catalog[i][7]) for i in range(len(catalog))]
    params={'ra':ra,'dec':dec,'redshift':redshift,'gamma':gamma,'flux':flux,'lum':lum}
    return params

params = cache.get (datafolder + 'SwiftBAT70m/pickle/params.pickle', get_SwiftBAT_params)

## These params contain everything we should need to weight our sources. ##

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

import datatest_stack_bins

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

llh_uniform = datatest_stack_bins.init(energy=True)

#llh_flux = datatest_stack.init(energy=True, weighting = flux)

#llh_redshift = datatest_stack.init(energy=True, weighting = np.power(redshift,-2))

#For this part it's unnecessary to populate signal, but remember how to do so down the line...

# init likelihood class
   # llh = datatest_stack.init(energy=True)
   # mc = datatest_stack.MC()
   # extra=datatest_stack.extra()
   # dpsi=extra["dpsi"]
   # print llh

##Okay, so the following is the part where we need to split this up into parallel processing. I think the pertinant variable to use here is n_iter... let's test a background scramble with n_iter=5 to see how fast it goes. Though, maybe it's max_iter? check with previously pickled results to see which number of bckg trials we got.


bckg_trials_uniform = PointSourceLLH.background_scrambles(llh_uniform,src_ra,src_dec,alpha=0.5,maxiter=batchsize) 


## And let's also cache bckg_trials for the other weighting schemes. ##

#
#def get_flux_bckg_trials():
#    return PointSourceLLH.background_scrambles(llh_flux,src_ra,src_dec,alpha=0.5)
#
#bckg_trials_flux = cache.get (datafolder + 'SwiftBAT70m/pickle/bckg_trials_flux.pickle', get_flux_bckg_trials)
#
#
#def get_redshift_bckg_trials():
#    return PointSourceLLH.background_scrambles(llh_redshift,src_ra,src_dec,alpha=0.5)
#
#bckg_trials_redshift = cache.get (datafolder + 'SwiftBAT70m/pickle/bckg_trials_redshift.pickle', get_redshift_bckg_trials)
#

##This one's gonna work a little differently than the single source sensitivity. First off, we need to calculate the background scrambles ahead of time, with the definition provided in psLLH_stack.py. I'll try to do this all within the same function:##

## Background Trials have the following keys:##
##['beta', 'TS_beta', 'beta_err', 'n_inj', 'nsources', 'TS', 'gamma']##
## Let's use a uniform weight (none) first to yield our bckg trials. ##


#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/background_trials/')

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials_uniform, outfile)





