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
from scipy.stats import chi2, norm
import itertools
from mstacking.utils import poisson_weight, delta_chi2
from optparse import OptionParser
import argparse
#This toggle will cause me to use EnergyLLH unless an extension is specified (SNR) in which case ExtendedEnergyLLH is used.
ext=None

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")


parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--fom', dest = 'fom', type = float,
                default = 0.0, metavar = 'FOM',
                help = 'Figure of Merit cut for WHSP blazars.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--sirin', dest = 'sirin', type = int,
                default = 0, metavar = 'SIRIN',
                help = 'Tells whether the unsplined IC79 is used')

parser.add_option ('--mhubertest', dest = 'mhubertest', type = int, 
                default = 0, metavar = 'MHUBERTEST',
                help = 'Toggles mhuber code to take precedence')

parser.add_option ('--disc', dest = 'disc', type = int,
                default = 0, metavar = 'DISC',
                help = 'Toggles calculation of discovery potential instead of sensitivity.')

parser.add_option ('--mese', dest = 'mese', type = int,
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
fom = opts.fom
years = opts.n
sirin = opts.sirin
mhubertest = opts.mhubertest
disc = opts.disc
mese = opts.mese

##I'm gonna put these here as defaults just in case I need to do a quick check in ipython.
catalog = 'SNR_noPWN'
gammas = np.linspace(1.0,4.0,31)
llhweight = 'weight'
injweight = 'weight'
years = 7
mese = 1
sirin = 1

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == '4yr_Starburst':
        flux = np.array(params['S60m'])
        modelweights = {'flux':flux}
elif catalog == 'SNR_cloud' or catalog == 'SNR_noPWN' or catalog == 'SNR_PWN':
        weights = np.array(params['weight'])
        modelweights = {'weight':weights}
        ext = np.array(params['extension'])
        #For some reason using an extension requires src_dec and src_ra to explicitly be in arrays, not lists.
        src_dec = np.array(src_dec)
        src_ra = np.array(src_ra)
elif catalog == 'WHSP_Blazars':
        #Here we'll also apply the fom cut to the weights and source locations.
        FOMmask = (params['weight']>=fom)
        weights = np.array(params['weight'][FOMmask])
        src_ra=src_ra[FOMmask]
        src_dec=src_dec[FOMmask]
        modelweights = {'weights':weights, 'equal':np.ones_like(weights)}
elif catalog == '2LAC':
        flux = np.array(params['flux'])
        modelweights = {'flux':flux}



##load llhmodel and assoc. MC
if catalog == 'WHSP_Blazars' or mhubertest:
  import load_mstacking
  if mese:
    if sirin:
      llhmodel = load_mstacking.load_7yr_mese_nospline_79()
    else:
      llhmodel = load_mstacking.load_7yr_mese()
  else:
    if sirin:
      llhmodel = load_mstacking.load_7yr_nospline_79()
    else:
      llhmodel = load_mstacking.load_7yr()
  MC = load_mstacking.monte_carlo(llhmodel) 
else:
  import load_PS7yr as datascript
  if years == 7:
    if mese:
      llhmodel = datascript.load7yr_mese(energy=True, mode='all', sirin=sirin)
      MC = datascript.loadMC7yr_mese(sirin=sirin)
  #The following catalogs used sirins IC79 and no pull correct of 1.1774 for the first 4yrs of data.
  elif catalog == '4yr_Starburst' or catalog == 'SNR_noPWN' or catalog == 'SNR_cloud' or catalog =='SNR_PWN':
    if years == 3:
      llhmodel = datascript.load3yr(energy=True, mode='all', sirin=True, nopull=True)
      MC = datascript.loadMC_3yr_nopull()
    elif years == 4:
      llhmodel = datascript.load4yr(energy=True, mode='all', sirin=True, nopull=True)
      MC = datascript.loadMC_4yr_nopull()
  else:
    if years == 3:
      if catalog == '2LAC':
        llhmodel = datascript.load3yr_no40(energy=True, mode='all')
        MC = datascript.loadMC_3yr_no40()
      else:
        llhmodel = datascript.load3yr(energy=True, mode='all')
        MC = datascript.loadMC_3yr()
    elif years == 4:
      llhmodel = datascript.load4yr(energy=True, mode='all')
      MC = datascript.loadMC_4yr() 

#inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)


weights = [np.round(StackingMultiPointSourceLLH.source_weights(llhmodel,src_dec, gamma=gamma)[0],4) for gamma in gammas]


for (gamma,weight) in zip (gammas,weights):
  print("Gamma = " + str(gamma) + "     -   weights: " + str(weight))







