import time, argparse
import numpy as np
from config import config

###############################################################################

#############
# ARGUMENTS #
#############

p = argparse.ArgumentParser(description="Calculates TS for a given source"
                            "location using 7 year point source sample",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--name", default="J1908+06", type=str,
               help="Source Name (default='J1908+06'")
p.add_argument("--ra",   default=286.98, type=float,
               help="Right ascension [degrees] (default=286.98)")
p.add_argument("--dec",  default=6.27, type=float,
               help="Declination [degrees] (default=6.27)")
p.add_argument("--nscramble",  default=5000, type=int,
               help="Number of background scrambles to use when computing p-value (default=5000)")

args = p.parse_args()

name      = args.name
ra        = np.radians(args.ra)
dec       = np.radians(args.dec)
nscramble = args.nscramble

#############
# ARGUMENTS #
#############

###############################################################################

##########
# RESULT #
##########

t0 = time.time()

seasons  = ["IC40", "IC59", "IC79b", "IC86, 2011", "IC86, 2012-2014"]
llh, inj = config(seasons, scramble = False, dec = dec)

TS, Xmin = llh.fit_source(src_ra = ra, src_dec = dec)

ns    = Xmin['nsources']
gamma = Xmin['gamma']

print("\n*" + 78*"-" + "*")
print("| Source: %s" % name)
print("| Sample: All-Sky Point Source Tracks (7yr)")
print("| Location: RA %.2f deg, Dec %.2f deg" % (args.ra, args.dec))
print("|")
print("| best fit TS = %.2f {'nsources': %.2f, 'gamma': %.2f}" % (TS, ns, gamma))

##########
# RESULT #
##########

###############################################################################

###########
# P-VALUE #
###########

# estimate p-value using nscramble background scrambles this location

t1 = time.time()
trials = llh.do_trials(src_ra = ra, src_dec = dec, n_iter = nscramble)
dt = time.time()-t1

p = trials[ trials['TS'] >= TS ].size / float(trials.size)

if p > 0: msg = "| p-pretrial  = %.3f" % p
else:     msg = "| p-pretrial  < %.1e" % (1.0/trials.size)
msg += " (est. from %.1e scrambles in %.2f sec)" % (trials.size, dt)
print(msg)
print("*" + 78*"-" + "*")

############
# P-VALUES #
############

###############################################################################

print("\n Completed in %.2f sec" % (time.time()-t0))

