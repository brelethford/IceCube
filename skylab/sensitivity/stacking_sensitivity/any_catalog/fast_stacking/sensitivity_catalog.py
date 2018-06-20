#!/usr/bin/env python
from __future__ import print_function

import argparse, os, time
import numpy as np
import scipy
from icecube.umdtools import cache,misc
from config                   import config, fit_background, TeV
from skylab.sensitivity_utils import sensitivity_flux, fit
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
p.add_argument("--gamma", default=2.0, type=float,
                help="Injection spectrum.")
p.add_argument("--n", default=4, type=int,
                help="number of years.")
p.add_argument("--discovery_thresh", default=5, type=float,
                help="Discovery threshold in sigma (default=5)")
p.add_argument("--cat", default=None, type=str,
                help="Catalog used in stacking (default=None)")
p.add_argument("--weight", default='equal', type=str,
                help="Optional argument for determining the likelihood model weighting scheme.")

args = p.parse_args()
cat = args.cat
gamma = args.gamma
disc_thresh = args.discovery_thresh
n = int(args.n)
print(n)
## These params contain everything we should need to weight our sources. I've made sure the dec and ra are in radians. ##
params = cache.load ('/data/user/brelethford/Data/{}/pickle/params.pickle'.format (cat))
src_ra, src_dec =  params['ra'], params['dec']
  

if cat == 'SwiftBAT70m':
        weight = args.weight
        flux, redshift = params['flux'], params['redshift']
        modelweights = {'flux':np.array(flux), 'redshift': np.array(list(np.power(redshift,-2))), 'equal':np.array(list(np.ones_like(src_dec)))}
        src_w  = modelweights['{}'.format(weight)]
elif cat == '4yr_Starburst':
        weight = 'S60m'
        src_w = np.array(params['S60m'])
elif cat == 'SNR_noPWN' or cat == 'SNR_PWN' or cat == 'SNR_cloud':
        weight  = 'weight'
        src_w = np.array(params['weight'])
elif cat == 'teststack' or cat == 'teststack50' or cat == 'teststack300':
        weight  = 'equal'
        src_w = np.ones_like(src_dec)
elif cat == 'blackhole':
        weight = 'flux2m'
        src_w = np.array(params['flux2m'])
elif cat == 'Milagro17':
        weight = 'weight'
        src_w = np.array(params['weight'])
else:
        raise ValueError("This script is for stacking only - choose an appropriate catalog.")

if cat:
  tag = "%s_cat%s" % (args.sample, str(cat))
else:
  tag = "%s_dec%+06.2f" % (args.sample, np.degrees(src_dec))

msg  = "# sensitivity calculation for catalog {}".format(cat)

os.system("echo \"" + msg + " on $HOSTNAME " + time.strftime("(%b %d, %Y)\""))

#############
# ARGUMENTS #
#############

###############################################################################

##########
# SKYLAB #
##########

#determine sample by number of years:

if cat == 'Milagro17':
  n=1
  seasons = [ ["PointSourceTracks", "IC40"]]
elif n == 3:
  seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"]]
elif n == 4:
  seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"]]
elif n == 7:
  seasons = [ ["PointSourceTracks", "IC40"],
            ["PointSourceTracks", "IC59"],
            ["PointSourceTracks", "IC79"],
            ["PointSourceTracks", "IC86, 2011"],
            ["PointSourceTracks", "IC86, 2012-2014"] ]

llh, inj = config(seasons, gamma=gamma, seed=args.seed, dec=src_dec, src_w = src_w)
print ("mu2flux (1) calculation: {}".format(inj.mu2flux(1)))
##########
# SKYLAB #
##########

outputfolder = '/data/user/brelethford/Output/fast/'
sens_dir = misc.ensure_dir (outputfolder+'stacking_sensitivity/{0}yr/{1}/{2}/'.format(n,cat,weight))
bg_dir = sens_dir+'bg/'
sig_dir = misc.ensure_dir(sens_dir+'sig/gamma_{}/'.format(gamma))

#For Plots:
dir_name = misc.ensure_dir("/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/{}/{}/".format(n,cat,weight,gamma))
subprocess.call(['mkdir', '-p', '{0}'.format(dir_name)])
###############################################################################

###################################
# BACKGROUND ONLY TS DISTRIBUTION #
###################################

print ("\nRunning background only trials ...")
print ("Grabbing background trials from {}...".format(bg_dir))
batches = [cache.load(bg_dir+file) for file in os.listdir(bg_dir) if file.endswith('.array')] 
trials = np.concatenate(batches,axis=1)


#save trials and fit params
outfile = sens_dir + 'trials.array'
cache.save(trials,outfile)

ts_parameters = fit_background(trials, dir_name + "background_trials_%s.png" % tag)
print ("Finished background fit.")

outfile = sens_dir + 'bg_fit.array'
cache.save(ts_parameters,outfile)

#These get used below for fitting the disc. pot.
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

print ("Concatenating injection trials from {}...".format(sig_dir))
n_injs = []
trials = []
for n_inj in [ns for ns in os.listdir(sig_dir) if ns not in ['inj_trials.array', 'sens.array','disc.array']]:
  n_injs.append(np.float(n_inj))
  batches = [cache.load(sig_dir+n_inj+'/'+file) for file in os.listdir(sig_dir+n_inj+'/') if file.endswith('.array')] 
  trials.append(np.concatenate(batches,axis=1))
inj_trials = dict(zip(n_injs,trials))
print ("Finished concatenation.")

#save trials
outfile = sig_dir + 'inj_trials.array'
cache.save(inj_trials, outfile)

#grab the function to fit TS for disc. pot.
def sigma2ts(s, ndf = 1.0, eta = 0.5, scale = 1.0):
  p = sigma2p(s)
  return p2ts(p, ndf, eta, scale)

def p2ts(p, ndf = 1.0, eta = 0.5, scale = 1.0):
  r""" Convert 1-sided p-value (probability for observing >= x)

  Parameters
  ----------
  p : float
    1-sided p-value
  ndf : float
    effective degrees of freedom, should match number of free parameters
    in the likelihood in the limit of large statistics
  eta : float
    fraction of TS > 0
  """
  return scale*scipy.stats.chi2.isf(p/eta, ndf)


def sigma2p(s):
  r""" Convert 1-sided significance defined by integral of normal distribution
       from s to infinity to 1-sided p-value.
  """
  return 0.5*scipy.special.erfc(s/np.sqrt(2))

ts_sens = median_ts 
ts_disc = sigma2ts(disc_thresh,ndf,eta,scale)

print ("The TS median is: {}".format(ts_sens))
print ("The TS 5 sigma threshold is: {}".format(ts_disc))

#Now just copy what's in the sensitivity flux calculation in sensitivity_utils.py.
results_sens = []
results_disc = []

#Apply thresholds for each one by hand (consult the relevant plots to do this:)

if cat == 'Milagro17':
  sens_thresh = 20
  disc_thresh = 80
elif cat == '4yr_Starburst':
  sens_thresh = 30
  disc_thresh = 100
elif cat == 'blackhole':
  sens_thresh = 40
  disc_thresh = 100
elif cat in ['SNR_cloud','SNR_PWN','SNR_noPWN']:
  sens_thresh = 16
  disc_thresh = 60
elif cat == 'SwiftBAT70m':
  sens_thresh = 60
  disc_thresh = 180

for mu in inj_trials.viewkeys():
  if mu < sens_thresh:
    trials = inj_trials[mu]
    n    = trials['TS'][ trials['TS'] > ts_sens ].size
    ntot = trials['TS'].size
    results_sens.append([mu, n, ntot])

for mu in inj_trials.viewkeys():
  if mu < disc_thresh:
    trials = inj_trials[mu]
    n    = trials['TS'][ trials['TS'] > ts_disc ].size
    ntot = trials['TS'].size
    results_disc.append([mu, n, ntot])

results_sens = np.transpose(results_sens)
results_disc = np.transpose(results_disc)

#I don't want to have to deal with the rest... let's just input these results into the est. sensitivity fcn halfway through.

results_sens = sensitivity_flux(ts_sens, 0.9, llh, inj, fguess=0,
                               factor = 2, npar_fit = 2, par_fit = None,
                               results = results_sens,
                               nstep = args.nstep, 
                           ni_bounds = args.ni_bounds, 
                       ts_parameters = ts_parameters, # see BACKGROUND ONLY TS DISTRIBUTION
                             nsample = args.nsample, 
                                name = "split_sensitivity",
                                path = dir_name + "sensitivity_%s_" % tag)

outfile = sig_dir + 'sens.array'
cache.save(results_sens,outfile)

results_disc = sensitivity_flux(ts_disc, 0.5, llh, inj, fguess=0,
                               factor = 2, npar_fit = 2, par_fit = None,
                               results = results_disc,
                               nstep = args.nstep, 
                           ni_bounds = args.ni_bounds, 
                       ts_parameters = ts_parameters, # see BACKGROUND ONLY TS DISTRIBUTION
                             nsample = args.nsample, 
                                name = "split_sensitivity",
                                path = dir_name + "discovery_potential_%s_" % tag)

outfile = sig_dir + 'disc.array'
cache.save(results_disc,outfile)

mu2flux_str = ("| mu2flux for 1 event @ %.1f TeV is %.2e |" % 
                      (inj.E0/TeV, inj.mu2flux(1)*TeV))
sensitivity_str = ("| Sensitivity flux @ %.1f TeV is %.2e +/- %.2e TeV^-1cm^-2s^-1 |" % 
                      (inj.E0/TeV, results_sens[1]*TeV,  results_sens[2]*TeV))
discovery_str =   ("| Discovery flux   @ %.1f TeV is %.2e +/- %.2e TeV^-1cm^-2s^-1 |" % 
                      (inj.E0/TeV, results_disc[1]*TeV, results_disc[2]*TeV))

print ("\n*" + "-" * (len(sensitivity_str)-2) + "*")
print (mu2flux_str)
print (sensitivity_str)
print (discovery_str)
print ("*" + "-" * (len(sensitivity_str)-2) + "*")

########################
# ESTIMATE SENSITIVITY #
########################

#copy folder to public

#skylab_folder = "/home/public/brelethford/skylab/fast/plots/stacking/{}yr/{}/{}/{}/".format(n,cat,weight,gamma)
#subprocess.call(['mkdir', '-p', '{0}'.format(skylab_folder)])
#subprocess.call(['cp', '*', '{}'.format(dir_name), '{}'.format(skylab_folder)])
#print ('Saving images to {}'.format(skylab_folder))

###############################################################################
