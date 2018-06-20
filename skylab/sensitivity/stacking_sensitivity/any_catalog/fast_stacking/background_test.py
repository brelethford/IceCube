#!/usr/bin/env python
from __future__ import print_function

import argparse, os, time
import numpy as np
from icecube.umdtools import cache,misc
from config                   import config, fit_background, TeV
from skylab.sensitivity_utils import estimate_sensitivity
import subprocess

###############################################################################

#############
# ARGUMENTS #
#############

p = argparse.ArgumentParser(description="Calculates Sensitivity and Discovery"
                            " Potential Fluxes for Point Source Analysis.",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--sample", default="PointSourceTracks4yr", type=str,
                help="Sample Name (default=\"PointSourceTracks4yr\")")
p.add_argument("--nscramble", default=10, type=int,
                help="Number of background only scrambles used to measure TS distribution (default=1000)")
p.add_argument("--batch", default=0, type=int,
                help="batch of background trials")
p.add_argument("--seed", default=1, type=int,
                help="Seed for RNG (default=1)")
p.add_argument("--ra", default=180., type=float,
                help="Right ascension of source in degrees (default=180)")
p.add_argument("--sindec", default=0., type=float,
                help="Declination of source in sindec of radians (default=0)")
p.add_argument("--cat", default=None, type=str,
                help="Catalog used in stacking (default=None)")
p.add_argument("--weight", default='equal', type=str,
                help="Optional argument for determining the likelihood model weighting scheme.")

args = p.parse_args()
cat = args.cat
batch = args.batch

## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##
if cat:
  #Grab the right catalog.
  params = cache.load ('/data/user/brelethford/Data/{}/pickle/params.pickle'.format (cat))
  src_ra, src_dec =  params['ra'], params['dec']
  try:
    src_ext = params['ext']
  except:
    src_ext = np.zeros_like(src_dec)

else:
  src_ra, src_dec = np.radians(args.ra), np.arcsin(args.sindec)
  

if cat == 'SwiftBAT70m':
        weight = args.weight
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'equal':np.array(list(np.ones_like(src_dec)))}
        src_w  = modelweights['{}'.format(weight)]
elif cat == 'Milagro17':
        weight = 'weight'
        src_w = np.array(params['weight'])
elif cat == '4yr_Starburst':
        weight = 'S60m'
        src_w = np.array(params['S60m'])
elif cat == 'teststack' or cat == 'teststack50' or cat == 'teststack300':
        weight  = 'equal'
        src_w = np.ones_like(src_dec)
elif cat == 'SNR_noPWN' or cat == 'SNR_PWN' or cat == 'SNR_cloud':
        weight  = 'weight'
        src_w = np.array(params['weight'])
elif cat == 'blackhole':
        weight = 'flux2m'
        src_w = np.array(params['flux2m'])
else:
        raise ValueError('Must provide a catalog')

if cat:
  tag = "%s_cat%s" % (args.sample, str(cat))
else:
  tag = "%s_dec%+06.2f" % (args.sample, np.degrees(src_dec))

msg  = "# output from \\`python background_catalog.py"
msg += " --sample %s --nscramble %d " % \
       (args.sample, args.nscramble)
msg += " --seed %d --dec %.3f\\`" % \
       (args.seed, np.degrees(args.sindec))

os.system("echo \"" + msg + " on $HOSTNAME " + time.strftime("(%b %d, %Y)\""))

#############
# ARGUMENTS #
#############

###############################################################################

##########
# SKYLAB #
##########
if cat == 'Milagro17':
  seasons = [ ["PointSourceTracks", "IC40"]]
  n = 1
elif args.sample == "PointSourceTracks4yr":
  seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"]]
  n = 4
elif args.sample == "PointSourceTracks7yr":
  seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"],
            ["PointSourceTracks", "IC86, 2012-2014"]]
  n = 7

llh, inj = config(seasons, gamma=2.0, seed=args.seed, dec=src_dec, src_w = src_w)
'''
##########
# SKYLAB #
##########

outputfolder = '/data/user/brelethford/Output/fast/'
sens_dir = misc.ensure_dir (outputfolder+'stacking_sensitivity/{0}yr/{1}/{2}/bg/'.format(n,cat,weight))

###############################################################################

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

print ("\nRunning background only trials ...")
def get_background():
  t0     = time.time()
  trials = llh.do_trials(args.nscramble, src_ra=src_ra, src_dec=src_dec, src_extension = src_ext, src_w = src_w)
  dt     = time.time() - t0
  print (" Completed %d trials in %.2f sec (%.2f trials/sec)" % (trials.size, dt, trials.size/dt))
  return trials

#save trials and fit params
outfile = sens_dir + 'trials_batch_{}.array'.format(batch)
trials = cache.get(outfile, get_background)

print("\nBackground batch {} finished.".format(batch))

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

'''
