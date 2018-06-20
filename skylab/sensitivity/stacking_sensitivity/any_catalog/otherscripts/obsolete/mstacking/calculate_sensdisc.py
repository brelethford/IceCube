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
catalog = '4yr_Starburst'
Gamma = 2.0
llhweight = 'flux'
injweight = 'flux'
years = 7


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


#Load the background trials and calculate the background TS thresholds.
if mhubertest:
  if mese:
    backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}_mhubertest/background_trials/'.format(catalog,str(years),llhweight)
  else:
    backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}_mhubertest/background_trials/'.format(catalog,str(years),llhweight)
elif catalog == 'WHSP_Blazars':
  backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/FOM_{2}/{3}/background_trials/'.format(catalog, str(years), str(fom), llhweight)
else:
  if mese:
    backfolder = '/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}/background_trials/'.format(catalog, str(years), llhweight)
  else:
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

inj = StackingPointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)

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

sensitivity = StackingMultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,ext=ext,alpha=alpha,beta=beta,inj=inj,mc=MC,TSval = TSval,w_theoMC=modelweights['{}'.format(llhweight)],w_theo=modelweights['{}'.format(injweight)],eps=0.04,n_iter=100)
print sensitivity

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))



if mhubertest:
  if mese:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}_mhubertest/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight, calctype))
  else:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}_mhubertest/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight, calctype))
else:
  if catalog == 'WHSP_Blazars':
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/FOM_{2}/{3}/{4}_inj/{5}/'.format(catalog,str(years),str(fom),llhweight,injweight, calctype))
  elif mese:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))
  else:
    sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))

# save the output
outfile = sens_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





