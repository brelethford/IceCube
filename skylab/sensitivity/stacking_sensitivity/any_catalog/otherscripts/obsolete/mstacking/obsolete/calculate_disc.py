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

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")
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


opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
years = opts.n
fom = opts.fom

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':flux, 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == '4yr_Starburst':
        flux = np.array(params['S60m'])
        modelweights = {'flux':flux}
elif catalog == '30youngSNR':
        weights = np.array(params['weights'])
        modelweights = {'weights':weights}
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

#We'll need the background trials.
if catalog == 'WHSP_Blazars':
  backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/FOM_{2}/{3}/background_trials/'.format(catalog, str(years), llhweight, str(fom))
else:
  backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight)

backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight)

def getTSval():
  files = [cache.load(backfolder+file) for file in os.listdir(backfolder) if file.endswith('.array')]
  fitfun = delta_chi2
  TS=[]
  for file in files:
    for item in range(len(file['n_inj'])):
      if file['n_inj'][item] ==0:
        TS.append(file['TS'][item])
  TSs=TS
  fit = fitfun(TSs, df=2., floc=0., fscale=1.)
  TSval = np.asscalar(fit.isf(norm.sf(5)))
  return TSval

TSval = cache.get(backfolder + 'TSval5sig.pickle', getTSval)

#Now I've defined all the possible modelweights for each catalog (I can add new ones later). Now I'll read in the years of data.
if catalog == 'WHSP_Blazars':
  import load_mstacking
  llhmodel = load_mstacking.load_7yr()
  MC = load_mstacking.monte_carlo(llhmodel)
else:
  if catalog == '4yr_Starburst':
    import load_starburst as datascript
  else:
    import load_PS7yr as datascript
  llh40= datascript.init40(energy=True,mode='all')
  llh59= datascript.init59(energy=True,mode='all')
  llh79 = datascript.init79(energy=True,mode='all')
  llh86I= datascript.init86I(energy=True,mode='all')
  MC = datascript.initMC_4yr()
  if years == 1:
    samples = [llh40]
    MC.pop(1,2,3)
  elif years == 2:
    samples = [llh40,llh59]
    MC.pop(2,3)
  elif years == 3:
    #2LAC blazars doesn't use IC40.
    if catalog == '2LAC':
      samples = [llh59,llh79,llh86I]
      MC.pop(0)
    else:
      samples = [llh40,llh59,llh79]
      MC.pop(3)
  elif iyears == 4:
    samples = [llh40,llh59,llh79,llh86I]
  llhmodel = datascript.multi_init(samples,energy=True)

# We use however many years of data are necessary, and pop out the extra MC.

inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)

sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=norm.sf(5),beta=.5,inj=inj,mc=MC,TSval = TSval,w_theoMC=modelweights['{}'.format(llhweight)],w_theo=modelweights['{}'.format(injweight)],eps=0.03,n_iter=100)
print sensitivity

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

if catalog == 'WHSP_Blazars':
  sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/FOM{3}/{4}_inj/sensitivity/'.format(catalog,str(years),llhweight,str(fom),injweight))
else:
  sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/{3}_inj/discovery/'.format(catalog,str(years),llhweight,injweight))

# save the output
outfile = sens_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





