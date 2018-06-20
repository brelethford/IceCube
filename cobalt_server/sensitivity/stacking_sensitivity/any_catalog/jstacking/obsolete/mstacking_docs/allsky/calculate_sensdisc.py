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
##I'm gonna put these here as defaults just in case I need to do a quick check in ipython.

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--sindec', dest = 'sindec', type = float,
                default = 0.0, metavar = 'SINDEC',
                help = 'Sindec of the hypothesized source')

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
years = opts.n
sirin = opts.sirin
mhubertest = opts.mhubertest
disc = opts.disc
mese = opts.mese
datafolder='/data/user/brelethford/Data/'
sindec= opts.sindec

##  I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  np.array(0.0), np.array(np.arcsin(sindec))

if mhubertest:
  if mese:
    backfolder = '/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/mhubertest/background_trials/'.format(str(years))
  else:
    backfolder = '/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/mhubertest/background_trials/'.format(str(years))
else:
  if mese:
    backfolder = '/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/background_trials/'.format(str(years))
  else:
    backfolder = '/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/background_trials/'.format(str(years))

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
  fit = fitfun(TSs, df=2., floc=0., fscale=1.)
  TSval5sig = np.asscalar(fit.isf(norm.sf(5)))
  return TSval5sig

TSval_median = cache.get(backfolder + 'TSval.pickle', getTSval)
TSval5sig = cache.get(backfolder + 'TSval5sig.pickle', getTSval5sig)

##load llhmodel and assoc. MC
import load_PS7yr as datascript
if years == 7:
  if mese:
    llhmodel = datascript.load7yr_mese(energy=True, stacking = False, mode='all')
    MC = datascript.loadMC7yr_mese()
elif years == 3:
  llhmodel = datascript.load3yr(energy=True, stacking = False, mode='all')
  MC = datascript.loadMC_3yr()
elif years == 4:
  llhmodel = datascript.load4yr(energy=True, stacking = False, mode='all')
  MC = datascript.loadMC_4yr() 

inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec=src_dec, seed=0)

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

sensitivity = MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,ext=ext,alpha=alpha,beta=beta,inj=inj,mc=MC,TSval = TSval,eps=0.04,n_iter=100)
print sensitivity

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))



if mhubertest:
  if mese:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/mhubertest/{1}/'.format(str(years), calctype))
  else:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/mhubertest/{1}/'.format(str(years), calctype))
else:
  if mese:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/{1}/'.format(str(years),calctype))
  else:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/{1}/'.format(str(years),calctype))

# save the output
outfile = sens_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





