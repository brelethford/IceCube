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
#ext is just a toggle that triggers using an extended source llh.
ext=None

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 30000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--fom', dest = 'fom', type = float,
                default = 0.0, metavar = 'FOM',
                help = 'Figure of Merit cut for WHSP blazars.')

parser.add_option ('--sirin', dest = 'sirin', type=int,
                metavar = 'SIRIN', default=0,
                help = 'Determines if the unsplined IC79 is used.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--years', dest = 'years', type = int,
                default = 4, metavar = 'YEARS',
                help = 'Number of years of data')

parser.add_option ('--mhubertest', dest = 'mhubertest', type=int,
                default = 0, metavar = 'MHUBERTEST',
                help = 'Toggles mhuber code to take precedence')

parser.add_option ('--mese', dest = 'mese', type = int,
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
llhweight = opts.llhweight
injweight = opts.injweight
catalog = opts.catalog
years = opts.years
fom = opts.fom
sirin = opts.sirin
mhubertest = opts.mhubertest
mese = opts.mese

datafolder='/data/user/brelethford/Data/'

params = cache.load (datafolder + '{}/pickle/params.pickle'.format (catalog))

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##

src_ra, src_dec =  params['ra'], params['dec']

## Here's were we put the weights in depending on the catalog. If I do this right, I should only need one version of the script afterwards.

if catalog == 'SwiftBAT70m':
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':flux, 'redshift': np.array(list(np.power(redshift,-2))), 'uniform':np.array(list(np.ones_like(src_dec)))}
elif catalog == '4yr_Starburst':
        flux = params['S60m']
        modelweights = {'flux':flux}
elif catalog == 'SNR_cloud' or catalog == 'SNR_noPWN' or catalog == 'SNR_PWN':
        weights = np.array(params['weight'])
        modelweights = {'weight':weights}
        ext = np.array(params['extension'])
elif catalog == 'WHSP_Blazars':
        #Here we'll also apply the fom cut to the weights and source locations.
        FOMmask = (params['weight']>=fom)
        weights = np.array(params['weight'][FOMmask])
        src_ra=src_ra[FOMmask]
        src_dec=src_dec[FOMmask]
        modelweights = {'weights':weights, 'equal':np.ones_like(weights)}
elif catalog == '2LAC':
        flux = params['flux']
        modelweights = {'flux':flux, 'uniform':np.array(list(np.ones_like(src_dec)))}

#Now I've defined all the possible modelweights for each catalog (I can add new ones later). Now I'll read in the years of data.
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
else:
  import load_PS7yr as datascript
  if years == 7:
    if mese:
      llhmodel = datascript.load7yr_mese(energy=True, mode='all', sirin=sirin)
    else:
      llhmodel = datascript.load7yr(energy=True, mode='all', sirin=sirin)
  #The following catalogs used sirins IC79 and no pull correct of 1.1774 for the first 4yrs of data.
  elif catalog == '4yr_Starburst' or catalog == 'SNR_noPWN' or catalog == 'SNR_cloud' or catalog == 'SNR_PWN':
    if years == 3:
      llhmodel = datascript.load3yr(energy=True, mode='all', sirin=True, nopull=True)
    elif years == 4:
      llhmodel = datascript.load4yr(energy=True, mode='all', sirin=True, nopull=True)
  else:
    if years == 3:
      if catalog == '2LAC':
        llhmodel = datascript.load3yr_no40(energy=True, mode='all')
      else:
        llhmodel = datascript.load3yr(energy=True, mode='all')
    elif years == 4:
      llhmodel = datascript.load4yr(energy=True, mode='all')

bckg_trials = StackingMultiPointSourceLLH.do_trials(llhmodel,src_ra,src_dec,ext=ext,w_theo=np.array(modelweights['{}'.format(llhweight)]),n_iter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
if mhubertest:
  if mese:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}_mhubertest/background_trials/'.format(catalog,str(years),llhweight))
  else:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}_mhubertest/background_trials/'.format(catalog,str(years),llhweight))
elif catalog == 'WHSP_Blazars':
  out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/FOM_{2}/{3}/background_trials/'.format(catalog,str(years),str(fom),llhweight))
else:
  if mese:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/{2}/background_trials/'.format(catalog, str(years), llhweight))
  else:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight))

# save the output
outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)





