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

## There are three modelweights I can use, so lets put them in a dictionary for easy access. ##

modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2))}

import data

##Here I add an argument for the sole purpose of having separate background files differently named. ##
parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 1000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'uniform', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
llhweight = opts.llhweight

## We'll assign the proper weighting scheme for the search, then use it to calculate and cache the associated bckg trials: ##
if llhweight == 'uniform':
  llhmodel = data.init(energy=True)
else:
  llhmodel = data.init(energy=True, weighting = modelweights['{}'.format(llhweight)])

bckg_trials = PointSourceLLH.background_scrambles(llhmodel,src_ra,src_dec,alpha=0.5,maxiter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{}/background_trials/'.format(llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





