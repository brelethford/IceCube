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

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload")

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--batch', dest = 'batch', type = int,
                default = 0, metavar = 'BATCH',
                help = 'Assigns a number to each batch of background trials.')

parser.add_option ('--batchsize', dest = 'batchsize', type = int,
                default = 30000, metavar = 'BATCHSIZE',
                help = 'Assigns how many background trials are used in each batch.')

parser.add_option ('--sirin', dest = 'sirin', type=int,
                metavar = 'SIRIN', default=0,
                help = 'Determines if the unsplined IC79 is used.')

parser.add_option ('--sindec', dest = 'sindec', type = float,
                default = 0, metavar = 'SINDEC',
                help = 'The sin of the declination.')

parser.add_option ('--years', dest = 'years', type = int,
                default = 4, metavar = 'YEARS',
                help = 'Number of years of data')

parser.add_option ('--mhubertest', dest = 'mhubertest', type=int,
                default = 0, metavar = 'MHUBERTEST',
                help = 'Toggles mhuber code to take precedence')

parser.add_option ('--mese', dest = 'mese', type=int, 
                default = 0, metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
batch = opts.batch
batchsize = opts.batchsize
years = opts.years
sirin = opts.sirin
mhubertest = opts.mhubertest
mese = opts.mese
sindec = opts.sindec

datafolder='/data/user/brelethford/Data/'

src_ra, src_dec = 0.0, np.arcsin(sindec)

import load_PS7yr as datascript
if years == 7:
  if mese:
    llhmodel = datascript.load7yr_mese(energy=True, stacking=False, mode='all', sirin=sirin)
  else:
    llhmodel = datascript.load7yr(energy=True, stacking=False, mode='all', sirin=sirin)
  #The following catalogs used sirins IC79 and no pull correct of 1.1774 for the first 4yrs of data.
elif years == 3:
  llhmodel = datascript.load3yr(energy=True, stacking=False,mode='all')
elif years == 4:
  llhmodel = datascript.load4yr(energy=True, stacking=False, mode='all')

bckg_trials = MultiPointSourceLLH.do_trials(llhmodel,src_ra,src_dec,n_iter=batchsize)

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
if mhubertest:
  if mese:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/mhubertest/background_trials/'.format(str(years)))
  else:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/mhubertest/background_trials/'.format(str(years)))
else:
  if mese:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr_mese/background_trials/'.format(str(years)))
  else:
    out_dir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/background_trials/'.format(str(years)))


# save the output
outdir = misc.ensure_dir(out_dir + 'dec_{0:4f}/'.format(np.degrees(src_dec)))
outfile = outdir + 'background_batch_{0}.array'.format(batch))
print ('Saving', outfile, '...')
cache.save(bckg_trials, outfile)





