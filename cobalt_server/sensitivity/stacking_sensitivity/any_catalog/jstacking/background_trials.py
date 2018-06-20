#!/usr/bin/env python
import numpy as np
import time
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

start_time = time.time()

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

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--sirin', dest = 'sirin', type = int,
                default = 0, metavar = 'SIRIN',
                help = 'Tells whether the unsplined IC79 is used')

parser.add_option ('--mese', dest = 'mese', type = int,
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

#Opts for allsky
parser.add_option ('--sindec', dest = 'sindec', type = float,
                default = 0, metavar = 'SINDEC',
                help = 'The sin of the declination.')

#Opts for stacking
parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG', default = None,
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--datatype', dest = 'datatype', type = str,
                default = 'standard', metavar = 'DATA',
                help = 'determines whether my_data or npz_data is used')

opts, args = parser.parse_args ()
Gamma = opts.gamma
llhweight = opts.llhweight
catalog = opts.catalog
batch = opts.batch
datatype = opts.datatype
#just for now...:
batchsize = opts.batchsize
years = opts.n
sirin = opts.sirin
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
  elif catalog == 'Milagro17':
          weight = np.array(params['weight'])
          modelweights = {'weight':weight}
  elif catalog == 'teststack' or catalog == 'teststack50' or catalog == 'teststack300':
          weight = np.ones_like(src_dec)
          modelweights = {'equal':weight}
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

else:
  #setting up a single source version of the calculation below.
  src_dec = np.arcsin(opts.sindec)
  src_ra = 0.0

##load llh and assoc. MC
if datatype == 'standard':
  import load_j7yr_standard as datascript
elif datatype == 'npz':
  import load_j7yr_npz as datascript


if years == 7:
  if mese:
    llh = datascript.load7yr_mese()
  else:
    llh = datascript.load7yr()
elif years == 4:
  llh = datascript.load4yr()
elif years == 3:
  llh = datascript.load3yr()
elif years == 1:
  if catalog == 'Milagro17':
    print ("Using IC40")
    llh = datascript.loadIC40()
  else:
    print ("Using IC86")
    llh = datascript.load1yr()

if catalog:
  #if years == 1:
  #  bckg_trials = PointSourceLLH.do_trials(llh[0],n_iter = batchsize,src_ra=src_ra,src_dec=src_dec,src_w=modelweights['{}'.format(llhweight)])
  #else:
    bckg_trials = llh[0].do_trials(n_iter = batchsize,src_ra=src_ra,src_dec=src_dec,src_w=modelweights['{}'.format(llhweight)])
else:
    bckg_trials = llh[0].do_trials(n_iter = batchsize,src_ra=src_ra,src_dec=src_dec)
  

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))

outputfolder = '/data/user/brelethford/Output/{}/'.format(datatype)

if catalog:
  if mese:
    out_dir = misc.ensure_dir(outputfolder + 'stacking_sensitivity/{0}/{1}yr_mese/{2}/background_trials/'.format(catalog, str(years), llhweight))
  else:
    out_dir = misc.ensure_dir(outputfolder + 'stacking_sensitivity/{0}/{1}yr/{2}/background_trials/'.format(catalog, str(years), llhweight))
else:
  if mese:
    out_dir = misc.ensure_dir(outputfolder + 'allsky_sensitivity/{0}yr_mese/background_trials/dec_{1:4f}/'.format(str(years),np.degrees(src_dec)))
  else:
    out_dir = misc.ensure_dir(outputfolder + 'allsky_sensitivity/{0}yr/background_trials/dec_{1:4f}/'.format(str(years),np.degrees(src_dec)))

# save the output

outfile = out_dir + 'background_batch_{}.array'.format(batch)
print 'Saving', outfile, '...'
cache.save(bckg_trials, outfile)

print("Elapsed time: {}".format(time.time() - start_time))



