
from __future__ import print_function

import argparse, os, time
import numpy as np

from config                   import config, fit_background, TeV
from skylab.sensitivity_utils import estimate_sensitivity

###############################################################################

#############
# ARGUMENTS #
#############

p = argparse.ArgumentParser(description="Calculates Sensitivity and Discovery"
                            " Potential Fluxes for Point Source Analysis.",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--sample", default="PointSourceTracks7yr_epinat", type=str,
                help="Sample Name (default=\"PointSourceTracks7yr_epinat\")")
p.add_argument("--season", default="IC86, 2014", type=str,
                help="Season Name (default=\"IC86, 2014\")")
p.add_argument("--nscramble", default=10000, type=int,
                help="Number of background only scrambles used to measure TS distribution (default=10000)")
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
p.add_argument("--dec", default=0., type=float,
                help="Declination of source in degrees (default=0)")
p.add_argument("--erange", default=[0, np.inf], type=float, nargs=2,
                help="Energy range (default=0 np.inf)")

args = p.parse_args()

src_ra  = np.radians(args.ra)
src_dec = np.radians(args.dec)

tag = "%s_%s_dec%+06.2f" % (args.sample, args.season.replace(", ","_"), args.dec)
if args.erange != np.inf:
  tag += "_erange%.2f_%.2f" % (args.erange[0], args.erange[1])

msg  = "# output from \\`python sensitivity_1year.py"
msg += " --sample %s --season '%s'" % (args.sample, args.season)
msg += " --nscramble %d --nsample %d" % (args.nscramble, args.nsample)
msg += " --ni-bounds %d %d --nstep %d --seed %d --erange %.2f %.2f --dec %.3f\\`" % \
       (args.ni_bounds[0], args.ni_bounds[1], args.nstep, args.seed, args.erange[0], args.erange[1], args.dec)

os.system("echo \"" + msg + " on $HOSTNAME " + time.strftime("(%b %d, %Y)\""))

#############
# ARGUMENTS #
#############

###############################################################################

##########
# SKYLAB #
##########

llh, inj = config([[args.sample, args.season]], gamma=2.0, seed=args.seed,
                   dec=src_dec, e_range=(10**args.erange[0], 10**args.erange[1]))
##########
# SKYLAB #
##########

###############################################################################

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

print ("\nRunning background only trials ...")

t0     = time.time()
trials = llh.do_trials(args.nscramble, src_ra=src_ra, src_dec=src_dec)
dt     = time.time() - t0
print (" Completed %d trials in %.2f sec (%.2f trials/sec)" % (trials.size, dt, trials.size/dt))

ts_parameters = fit_background(trials, "figures/background_trials_%s.png" % tag)
median_ts, eta, ndf, scale = ts_parameters

'''
# GFU 2015 #
median_ts = 0.00000
eta       = 0.49929 
scale     = 0.99336
ndf       = 1.27374
ts_parameters = [median_ts, eta, ndf, scale]

# EPINAT 2014 #
median_ts = 0.00000
eta       = 0.48919
scale     = 0.95814
ndf       = 1.38249
ts_parameters = [median_ts, eta, ndf, scale]
'''

print("\nBackground only TS > 0 described by:")
print(" Median TS: %7.5f" % median_ts)
print(" PDF(TS>0): %7.5f * chi2(ts / %.5f, ndf = %.5f)" % (eta, scale, ndf))

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

###############################################################################

########################
# ESTIMATE SENSITIVITY #
########################

results = estimate_sensitivity(llh, inj,
                               nstep = args.nstep, 
                           ni_bounds = args.ni_bounds, 
                       ts_parameters = ts_parameters, # see BACKGROUND ONLY TS DISTRIBUTION
                             nsample = args.nsample, 
                              src_ra = src_ra, 
                             src_dec = src_dec,
                                path = "figures/sensitivity_%s_" % tag)

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
