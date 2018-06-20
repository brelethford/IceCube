#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from skylab.ps_llh import PointSourceLLH, MultiPointSourceLLH, FastMultiPointSourceLLH
from skylab.ps_injector import PointSourceInjector
from scipy.stats import chi2, norm
import itertools
from skylab.utils import poisson_weight, delta_chi2
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")
##I'm gonna put these here as defaults just in case I need to do a quick check in ipython.

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 30000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

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

parser.add_option ('--mese', dest = 'mese', type = int,
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
catalog = opts.catalog
batch = opts.batch
batchsize = opts.batchsize
fom = opts.fom
years = opts.n
sirin = opts.sirin
mese = opts.mese
datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == 'Milagro17':
        weight = np.array(params['weight'])
        modelweights = {'weight':weight}
elif catalog == '4yr_Starburst':
        flux = np.array(params['S60m'])
        modelweights = {'flux':flux}
elif catalog == 'blackhole':
        flux = np.array(params['flux2m'])
        modelweights = {'flux2m':flux}
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

##load llh and assoc. MC
if catalog == 'WHSP_Blazars':
  import load_mstacking
  if mese:
    if sirin:
      llh = load_mstacking.load_7yr_mese_nospline_79()
    else:
      llh = load_mstacking.load_7yr_mese()
  else:
    if sirin:
      llh = load_mstacking.load_7yr_nospline_79()
    else:
      llh = load_mstacking.load_7yr()
  MC = load_mstacking.monte_carlo(llh) 
else:
  import load_j7yr as datascript
  if years == 7:
    if mese:
      llh = datascript.load7yr_mese(mode='all', sirin=sirin)
    else:
      llh = datascript.load7yr(mode='all', sirin=sirin)
  #The following catalogs used sirins IC79 and no pull correct of 1.1774 for the first 4yrs of data.
  elif catalog == '4yr_Starburst' or catalog == 'blackhole' or catalog == 'SNR_noPWN' or catalog == 'SNR_cloud' or catalog =='SNR_PWN':
    if years == 3:
      llh = datascript.load3yr(mode='box', sirin=True, nopull=True)
    elif years == 4:
      llh = datascript.load4yr(mode='box', sirin=True, nopull=True)
  elif catalog == 'Milagro17':
    llh = datascript.init86I(mode='box')
  else:
    if years == 3:
      if catalog == '2LAC':
        llh = datascript.load3yr_no40(mode='all')
      else:
        llh = datascript.load3yr(mode='all')
    elif years == 4:
      llh = datascript.load4yr(mode='all')

if years == 1:
  bckg_trials = PointSourceLLH.do_trials(llh[0],n_iter = batchsize,src_ra=src_ra,src_dec=src_dec,src_w=modelweights['{}'.format(llhweight)])
else:
  bckg_trials = FastMultiPointSourceLLH.do_trials(llh[0],n_iter = batchsize,src_ra=src_ra,src_dec=src_dec,src_w=modelweights['{}'.format(llhweight)])

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

if catalog == 'WHSP_Blazars':
    out_dir = misc.ensure_dir('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/FOM_{2}/{3}/background_trials/'.format(catalog, str(years), str(fom), llhweight))
else:
  if mese:
    out_dir = misc.ensure_dir('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr_mese/{2}/background_trials/'.format(catalog, str(years), llhweight))
  else:
    out_dir = misc.ensure_dir('/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





