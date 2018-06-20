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
cat = '4yr_Starburst'

params = cache.load ('/data/user/brelethford/Data/{}/pickle/params.pickle'.format (cat))
src_ra, src_dec =  params['ra'], params['dec']
  

weight = 'S60m'
src_w = np.array(params['S60m'])

seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"]]#,
            #["PointSourceTracks", "IC86, 2012-2014"] ]


llh, inj = config(seasons, gamma=2.0, seed=args.seed, dec=src_dec, src_w = src_w)


##########
# SKYLAB #
##########


###############################################################################

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

print ("\nRunning background only trials ...")
trials = llh.do_trials(1, src_ra=src_ra, src_dec=src_dec, src_w = src_w)

print("Background over - starting injections:")

trials = llh.do_trials(1, injector=inj, mean_signal=5, src_ra=src_ra, src_dec=src_dec, src_w = src_w)

