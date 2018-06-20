#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from skylab.ps_llh import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector import PointSourceInjector
from scipy.stats import chi2, norm
import itertools
from skylab.utils import poisson_weight, delta_chi2
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")
##I'm gonna put these here as defaults just in case I need to do a quick check in ipython.

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--sirin', dest = 'sirin', type = int,
                default = 0, metavar = 'SIRIN',
                help = 'Tells whether the unsplined IC79 is used')

parser.add_option ('--disc', dest = 'disc', type = int,
                default = 0, metavar = 'DISC',
                help = 'Toggles calculation of discovery potential instead of sensitivity.')

parser.add_option ('--mese', dest = 'mese', type = int,
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

#opts for allsky
parser.add_option ('--sindec', dest = 'sindec', type = float,
                default = 0, metavar = 'SINDEC',
                help = 'The sin of the declination.')

#opts for stacking
parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--datatype', dest = 'datatype', type = str,
                default = 'my_data', metavar = 'DATA',
                help = 'Chooses my datatype or npz datatype')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
years = opts.n
datatype = opts.datatype
#temp:
sirin = opts.sirin
disc = opts.disc
mese = opts.mese

datafolder='/data/user/brelethford/Data/'

if catalog:
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
  elif catalog == 'teststack' or catalog == 'teststack50' or catalog == 'teststack300':
          weight = np.ones_like(src_dec)
          modelweights = {'equal':weight}
  elif catalog == 'blackhole':
          flux = np.array(params['flux2m'])
          modelweights = {'flux2m':flux}
  elif catalog == 'Milagro17':
          weight = np.array(params['weight'])
          modelweights = {'weight':weight}
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

else:
  #setting up a single source version of the calculation below.
  src_dec = np.arcsin(opts.sindec)
  src_ra = 0.0

#Load the background trials and calculate the background TS thresholds.
outputfolder = '/data/user/brelethford/Output/{}/'.format(datatype)

if catalog:
  if mese:
    backfolder = outputfolder+'stacking_sensitivity/{0}/{1}yr_mese/{2}/background_trials/'.format(catalog, str(years), llhweight)
  else:
    backfolder = outputfolder+'stacking_sensitivity/{0}/{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight)
else:
  if mese:
    backfolder = outputfolder+'allsky_sensitivity/{0}yr_mese/background_trials/dec_{1:4f}/'.format(str(years),np.degrees(src_dec))
  else:
    backfolder = outputfolder+'allsky_sensitivity/{0}yr/background_trials/dec_{1:4f}/'.format(str(years),np.degrees(src_dec))

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

TSval_median = cache.get(backfolder + 'TSval.pickle', getTSval)
TSval5sig = cache.get(backfolder + 'TSval5sig.pickle', getTSval5sig)

##load llh and assoc. MC
if datatype == 'my_data':
  import load_j7yr as datascript
elif datatype == 'npz_data':
  import load_j7yr_npz as datascript

catalogs_4yr = ['4yr_Starburst', 'blackhole', 'SNR_noPWN', 'SNR_PWN', 'SNR_PWN', 'SwiftBAT70m', 'teststack', 'teststack50', 'teststack300']

if catalog in catalogs_4yr:
    llh = datascript.load4yr(mode='box', sirin=True, nopull=True)
elif catalog == 'Milagro17':
    llh = datascript.init86I(mode='box')
else:
    if years == 7:
      if mese:
        llh = datascript.load7yr_mese(mode='box', sirin=sirin)
      else:
        llh = datascript.load7yr(mode='box', sirin=sirin)
    elif years == 4:
      llh = datascript.load4yr(mode='box', sirin=True, nopull=True)
    elif years == 1:
      print ("Using IC86")
      llh = datascript.init86I(mode='box')



llhmodel, exp, MC = llh[0], llh[1], llh[2]
inj = PointSourceInjector(Gamma, seed=0)
inj.fill(src_dec, exp, MC, livetime = llhmodel.livetime)

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

if catalog:
  if years == 1:
    sensitivity = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,src_w=modelweights['{}'.format(injweight)],eps=0.04,n_iter=100)
  else:
    sensitivity = MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,src_w=modelweights['{}'.format(injweight)],eps=0.04,n_iter=100)
else:
  if years == 1:
    sensitivity = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,eps=0.04,n_iter=100)
  else:
    sensitivity = MultiPointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=alpha,beta=beta,inj=inj,TSval = TSval,eps=0.04,n_iter=100)

print sensitivity

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

if catalog:
  if mese:
    sens_dir = misc.ensure_dir (outputfolder+'stacking_sensitivity/{0}/{1}yr_mese/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))
  else:
    sens_dir = misc.ensure_dir (outputfolder+'stacking_sensitivity/{0}/{1}yr/{2}/{3}_inj/{4}/'.format(catalog,str(years),llhweight,injweight,calctype))
else:
  if mese:
    sens_dir = misc.ensure_dir (outputfolder+'allsky_sensitivity/{0}yr_mese/{1}/dec_{2:4f}/'.format(str(years),calctype,np.degrees(src_dec)))
  else:
    sens_dir = misc.ensure_dir (outputfolder+'allsky_sensitivity/{0}yr/{1}/dec_{2:4f}/'.format(str(years),calctype,np.degrees(src_dec)))

# save the output
outfile = sens_dir + 'gamma{}.array'.format(Gamma)
print 'Saving', outfile, '...'
cache.save(sensitivity, outfile)





