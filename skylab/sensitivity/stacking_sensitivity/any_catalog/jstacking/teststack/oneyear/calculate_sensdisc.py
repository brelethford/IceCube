#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from jstacking.ps_llh import PointSourceLLH, MultiPointSourceLLH
from jstacking.ps_injector import PointSourceInjector
from jstacking.llh_models import ClassicLLH, EnergyLLH
from jstacking.utils import poisson_weight, delta_chi2
from scipy.stats import chi2, norm
import itertools
from optparse import OptionParser
import argparse

corrected=False

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")

parser = OptionParser (usage = '%prog [options]')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--disc', dest = 'disc', type = int,
                default = 0, metavar = 'DISC',
                help = 'Toggles calculation of discovery potential instead of sensitivity.')


opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
disc = opts.disc

datafolder='/data/user/coenders/data/MultiYearPointSource/npz'
catfolder='/data/user/brelethford/Data/{}'.format(catalog)

params = cache.load (catfolder + '/pickle/params.pickle')

years = 1 #for now...

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##
src_ra, src_dec =  params['ra'], params['dec']

#get TS median, chi2 from background:
def getTSval():
  files = [cache.load(backfolder+file) for file in os.listdir(backfolder) if file.endswith('.array')]
  fitfun = delta_chi2
  TS=[]
  for file in files:
    for item in range(len(file['n_inj'])):
      if file['n_inj'][item] ==0:
        TS.append(file['TS'][item])
  TSs=TS
  TSval = np.median(TSs)
  return TSval

def getTSval5sig():
  files = [cache.load(backfolder+file) for file in os.listdir(backfolder) if file.endswith('.array')]
  fitfun = delta_chi2
  TS=[]
  for file in files:
    for item in range(len(file['n_inj'])):
      if file['n_inj'][item] ==0:
        TS.append(file['TS'][item])
  TSs=TS
  eta = float(np.count_nonzero(TSs))/len(TSs)
  fit = fitfun(eta = eta, df=2.)
  TSval5sig = fit.isf(norm.sf(5))
  return TSval5sig

if corrected:
  backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr_corrected/{2}/background_trials/'.format(catalog, str(years), llhweight)
else:
  backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight)


TSval_median = cache.get(backfolder + 'TSval.pickle', getTSval)
TSval5sig = cache.get(backfolder + 'TSval5sig.pickle', getTSval5sig)

arr_exp = np.load(datafolder + '/IC86_exp.npy')
if corrected:
  arr_mc = np.load(datafolder + '/IC86_corrected_MC.npy')
else:
  arr_mc = np.load(datafolder + '/IC86_MC.npy')

#need llhmodel
dec_bins = np.unique(np.linspace(-1., 1, 100 + 1))

energy_bins = [np.linspace( 2.5,8.5,24+1),dec_bins]

llh_model = EnergyLLH(energy_bins, sinDec_bins=dec_bins)

llh = PointSourceLLH(arr_exp, arr_mc, livetime=332.61, llh_model=llh_model,
                         seed=0)#np.random.randint(2**32),

#need injector

inj = PointSourceInjector(Gamma, seed=0) 
inj.fill(src_dec,arr_exp,arr_mc,livetime = 332.61)

if disc:
  print("Calculating Discovery Potential...")
  alpha = norm.sf(5)
  beta = 0.5
  TSval = TSval5sig
  calctype = "disc"
else:
  print("Calculating Sensitivity...")
  alpha = 0.5
  beta = 0.9
  TSval = TSval_median
  calctype = "sens"

sensitivity = PointSourceLLH.weighted_sensitivity(llh,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,eps=0.02,n_iter=100)


#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
if corrected:
  out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr_corrected/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))
else:
  out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))
# save the output
outfile = out_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





