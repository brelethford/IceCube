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

corrected = False

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

years = 4 #for now...

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

exp40 = np.load(datafolder + '/IC40_exp.npy')
exp59 = np.load(datafolder + '/IC59_exp.npy')
exp79 = np.load(datafolder + '/IC79b_exp.npy')
exp86 = np.load(datafolder + '/IC86_exp.npy')
if corrected:
  mc40 = np.load(datafolder + '/IC40_corrected_MC.npy')
  mc59 = np.load(datafolder + '/IC59_corrected_MC.npy')
  mc79 = np.load(datafolder + '/IC79b_corrected_MC.npy')
  mc86 = np.load(datafolder + '/IC86_corrected_MC.npy')
else:
  mc40 = np.load(datafolder + '/IC40_MC.npy')
  mc59 = np.load(datafolder + '/IC59_MC.npy')
  mc79 = np.load(datafolder + '/IC79b_MC.npy')
  mc86 = np.load(datafolder + '/IC86_MC.npy')

exp = {0:exp40,1:exp59,2:exp79,3:exp86}
MC = {0:mc40,1:mc59,2:mc79,3:mc86}

#need llhmodel
dec_bins = np.unique(np.linspace(-1., 1, 40 + 1))

energy_bins = [np.linspace( 2.5,8.5,24+1),dec_bins]

llh_model40 = EnergyLLH(energy_bins, sinDec_bins=dec_bins)
llh_model59 = EnergyLLH(energy_bins, sinDec_bins=dec_bins)
llh_model79 = EnergyLLH(energy_bins, sinDec_bins=dec_bins)
llh_model86 = EnergyLLH(energy_bins, sinDec_bins=dec_bins)

llh40 = PointSourceLLH(exp40, mc40, livetime=375.539, llh_model=llh_model40, mode = 'all',
                         seed=0)#np.random.randint(2**32),

llh59 = PointSourceLLH(exp59, mc59, livetime=348.138, llh_model=llh_model59, mode = 'all',
                         seed=0)#np.random.randint(2**32),

llh79 = PointSourceLLH(exp79, mc79, livetime=315.506, llh_model=llh_model79, mode = 'all',
                         seed=0)#np.random.randint(2**32),

llh86 = PointSourceLLH(exp86, mc86, livetime=332.61, llh_model=llh_model86, mode = 'all',
                         seed=0)#np.random.randint(2**32),
samples = [llh40,llh59,llh79,llh86]
years = len(samples)

llh=MultiPointSourceLLH()
for i in xrange(len(samples)):
   llh.add_sample(str(i), samples[i])
   print ("Adding sample" + str(samples[i]))

#need injector


inj = PointSourceInjector(Gamma, seed=0) 
inj.fill(src_dec, exp, MC, livetime = llh.livetime)

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
print ("Starting sensitivity calculation.")
sensitivity = MultiPointSourceLLH.weighted_sensitivity(llh,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,eps=0.04,n_iter=100)


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





