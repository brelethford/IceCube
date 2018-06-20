#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from mstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH, StackingMultiPointSourceLLH #What I'm using - uses energy info
from mstacking.ps_injector import PointSourceInjector, StackingPointSourceInjector
from scipy.stats import chi2
import itertools
from bstacking.utils import poisson_weight
from optparse import OptionParser
import argparse
#extended is just a toggle that triggers using an extended source llh.

extended=False

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")

catalog = 'SNR_noPWN'
n=4
llhweight = 'weight'
injweight = 'weight'
Gamma = 2.

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

weights = np.array(params['weight'])
modelweights = {'weight':np.array(weights)}
ext = np.array(params['extension'])

import load_PS7yr as datascript
llhmodel = datascript.load3yr(energy=True, mode='all', sirin=True, nopull=True)
MC = datascript.loadMC_3yr_nopull()
alpha = 0.5
beta = 0.9
TSval = 0.
calctype = "sens"

inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)

print(type(src_ra))
print(type(src_dec))
print(type(ext))

sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,ext=ext,alpha=alpha,beta=beta,inj=inj,mc=MC,TSval = TSval,w_theoMC=modelweights['{}'.format(llhweight)],w_theo=modelweights['{}'.format(injweight)],eps=0.04,n_iter=10)

print ("Background Trials:")
print (bckg_trials)

