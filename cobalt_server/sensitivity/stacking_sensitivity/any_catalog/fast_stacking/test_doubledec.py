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
p.add_argument("--nscramble", default=10000, type=int,
                help="Number of background only scrambles used to measure TS distribution (default=1000)")
p.add_argument("--nsample", default=100, type=int,
                help="Number of signal samples used to compute discovery potential (default=100)")
p.add_argument("--ni-bounds", default=[0, 25], nargs=2, type=float, dest="ni_bounds",
                help="Range in signal events to search for discovery potential flux (default= 0 25)")
p.add_argument("--nstep", default=11, type=int,
                help="Number of signal injection steps (default=11)")
p.add_argument("--seed", default=1, type=int,
                help="Seed for RNG (default=1)")
p.add_argument("--ra", default=180., type=float,
                help="Right ascension of source in degrees (default=180)")
p.add_argument("--sindec", default=0., type=float,
                help="Declination of source in sindec of radians (default=0)")
p.add_argument("--discovery_thresh", default=5, type=float,
                help="Discovery threshold in sigma (default=5)")
p.add_argument("--cat", default=None, type=str,
                help="Catalog used in stacking (default=None)")
p.add_argument("--weight", default='equal', type=str,
                help="Optional argument for determining the likelihood model weighting scheme.")

args = p.parse_args()
cat = 'doubledec'


## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##
src_ra, src_dec = np.radians([0.,0.]), np.arcsin([0.5,-0.5])
  

src_w = None

if cat:
  tag = "%s_cat%s" % (args.sample, str(cat))
else:
  tag = "%s_dec%+06.2f" % (args.sample, np.degrees(src_dec))

msg  = "# output from \\`python sensitivity_multiyear.py"
msg += " --sample %s --nscramble %d --nsample %d --discovery_thresh %.1f" % \
       (args.sample, args.nscramble, args.nsample, args.discovery_thresh)
msg += " --ni-bounds %d %d --nstep %d --seed %d --dec %.3f\\`" % \
       (args.ni_bounds[0], args.ni_bounds[1], args.nstep, args.seed, np.degrees(args.sindec))

os.system("echo \"" + msg + " on $HOSTNAME " + time.strftime("(%b %d, %Y)\""))

#############
# ARGUMENTS #
#############

###############################################################################

##########
# SKYLAB #
##########

seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"]]#,
            #["PointSourceTracks", "IC86, 2012-2014"] ]


if args.sample == "PointSourceTracks7yr+GFU":
  seasons.append(["GFU", "IC86, 2015-2017"])

llh, inj = config(seasons, gamma=2.0, seed=args.seed, dec=src_dec)

##########
# SKYLAB #
##########

outputfolder = '/data/user/brelethford/Output/fast/'
if cat:
  sens_dir = misc.ensure_dir (outputfolder+'stacking_sensitivity/{0}yr/{1}/'.format(len(seasons),cat))
else:
  sens_dir = misc.ensure_dir (outputfolder+'allsky_sensitivity/{0}yr/{1:.3f}/'.format(len(seasons),np.degrees(src_dec)))
###############################################################################

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

print ("\nRunning background only trials ...")
def get_background():
  if cat:
    print ("Grabbing background trials from {}...".format(sens_dir))
    batches = [cache.load(sens_dir+file) for file in os.listdir(sens_dir) if file.endswith('.array')] 
    trials = np.concatenate(batches,axis=1)
  else:
    t0     = time.time()
    trials = llh.do_trials(args.nscramble, src_ra=src_ra, src_dec=src_dec, src_w = src_w)
    dt     = time.time() - t0
    print (" Completed %d trials in %.2f sec (%.2f trials/sec)" % (trials.size, dt, trials.size/dt))
  return trials

if cat:
    dir_name = "/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/".format(len(seasons),cat)
    #dir_name = "/home/brelethford/public_html/skylab/fast/plots/stacking/{}yr/{}/{}/".format(len(seasons),cat,weight)
    subprocess.call(['mkdir', '-p', '{0}'.format(dir_name)])
else:
    dir_name = "/data/user/brelethford/AGN_Core/Plots/fast/plots/allsky/{0}yr/dec_{1:.3f}/".format(len(seasons),np.degrees(src_dec))
    #dir_name = "/home/brelethford/public_html/skylab/fast/plots/allsky/{0}yr/dec_{1:.3f}/".format(len(seasons),np.degrees(src_dec))
    subprocess.call(['mkdir', '-p', '{0}'.format(dir_name)])

def get_bg_fit():
    ts_parameters = fit_background(trials, dir_name + "background_trials_%s.png" % tag)
    print ("Finished background fit.")
    return ts_parameters

#save trials and fit params
outfile = sens_dir + 'trials.array'
trials = cache.get(outfile, get_background)

outfile = sens_dir + 'bg_fit.array'
ts_parameters = cache.get(outfile, get_bg_fit)

median_ts, eta, ndf, scale = ts_parameters



print("\nBackground only TS > 0 described by:")
print(" Median TS: %6.4f" % median_ts)
print(" PDF(TS>0): %6.4f * chi2(ts / %.4f, ndf = %.4f)" % (eta, scale, ndf))

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

###############################################################################

########################
# ESTIMATE SENSITIVITY #
########################

def get_results():
  results = estimate_sensitivity(llh, inj,
                               nstep = args.nstep, 
                           ni_bounds = args.ni_bounds, 
                       ts_parameters = ts_parameters, # see BACKGROUND ONLY TS DISTRIBUTION
                             nsample = args.nsample, 
                              src_ra = src_ra, 
                             src_dec = src_dec,
                             src_w   = src_w,
                         disc_thresh = args.discovery_thresh,
                                path = dir_name + "sensitivity_%s_" % tag)
  return results

outfile = sens_dir + 'results.array'
results = cache.get(outfile, get_results)


sensitivity_str = ("| Sensitivity flux @ %.1f TeV is %.2e +/- %.2e TeV^-1cm^-2s^-1 |" % 
                      (inj.E0/TeV, results['sensitivity_flux']*TeV,  results['sensitivity_flux_error']*TeV))
discovery_str =   ("| Discovery flux   @ %.1f TeV is %.2e +/- %.2e TeV^-1cm^-2s^-1 |" % 
                      (inj.E0/TeV, results['discovery_flux']*TeV, results['discovery_flux_error']*TeV))

print ("\n*" + "-" * (len(sensitivity_str)-2) + "*")
print (sensitivity_str)
print (discovery_str)
print ("*" + "-" * (len(sensitivity_str)-2) + "*")

########################
# ESTIMATE SENSITIVITY #
########################

###############################################################################
