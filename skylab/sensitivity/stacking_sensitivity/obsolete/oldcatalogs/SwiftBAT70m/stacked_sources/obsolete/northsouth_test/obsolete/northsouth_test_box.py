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

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'

## Input my commandline parameters ##

parser = OptionParser ( usage = '%prog [options]')
parser.add_option ('--sky', dest = 'sky', type = str,
                default = None, metavar = 'SKY',
                help = 'tells which sources to use (and which folders to reference)')

opts,args = parser.parse_args ()
Gamma=2

sky = opts.sky

if sky == None:
  raise NameError('Error - did not establish which part of the sky to test.')

#The sources are all going to be loaded in - what changes is which ones I'll use.
params = cache.load(datafolder + 'SwiftBAT70m/pickle/params.pickle')

###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/{}/background_trials/'.format(sky)

bckg_trials = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]

## These params contain everything we should need to weight our sources. Remember to convert src_dec to sindec ##
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
print (str(src_dec))
## Time to define my modelweights in a dictionary. ##

modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2))}

## Now to import my llh model framework. ##
import data_box

## Like in the background trials, we have to define which llhmodel to use. For this check I'll use flux because it's easiest.
llhmodel = data_box.init(energy=True, weighting = modelweights['flux'])

##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)

##Now, I'll input my injection weighting scheme from the commandline. ##

inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, theo_weight = modelweights['flux'], seed=0) 
#inj.e_range = [1.e4,1.e8]

sensitivity = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.02,n_iter=250)
print (sensitivity)


#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/')
#disc_dir =  misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}/{1}_inj/disc/'.format(llhweight, injweight))

# save the output
outfile_sens = sens_dir + '{}.array_box'.format(sky)
#outfile_disc = disc_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

