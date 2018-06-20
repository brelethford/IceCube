import time
import numpy as np
from config import config

###############################################################################

###########
# SOURCES #
###########

# list of notable sources from: https://arxiv.org/pdf/1609.04981.pdf

srcs = np.empty(11, dtype=[('name', np.unicode_, 16), ('ra_deg', float),
                ('dec_deg', float), ('ns', float), ('gamma', float)])

srcs[ 0] = ('Northern Hotspot',  32.20,  62.10, 32.6, 2.8)
srcs[ 1] = ('Southern Hotspot', 174.60, -39.30, 15.4, 2.9)
srcs[ 2] = ('J1908+06',         286.98,   6.27,  4.5, 2.0)
srcs[ 3] = ('Cyg A',            299.87,  40.73,  2.1, 1.4)
srcs[ 4] = ('Cyg X-3',          308.11,  40.96, 12.8, 4.0)
srcs[ 5] = ('1ES 1959+650',     300.00,  65.15, 15.4, 3.1)
srcs[ 6] = ('PKS 1406-076',     212.24,  -7.87,  7.3, 2.6)
srcs[ 7] = ('HESS J1303-631',   195.74, -63.20,  4.5, 2.3)
srcs[ 8] = ('PKS 2005-489',     302.37, -48.82,  0.9, 1.0)
srcs[ 9] = ('HESS J1616-508',   243.78, -51.40,  2.4, 4.0)
srcs[10] = ('HESS J1614-518',   243.58, -51.82,  2.2, 4.0)

###########
# SOURCES #
###########

###############################################################################

###########
# RESULTS #
###########

results = srcs[['name','ra_deg','dec_deg','ns']].copy()
results.dtype.names = ('name','ra_deg','dec_deg','TS')

t0 = time.time()

seasons = ["IC40", "IC59", "IC79b", "IC86, 2011", "IC86, 2012-2014"]
llh, inj = config(seasons, scramble=False)

print("\n results:")
for i, src in enumerate(srcs):

  TS, Xmin = llh.fit_source(src_ra  = np.radians(src['ra_deg']),
                            src_dec = np.radians(src['dec_deg']))

  ns = Xmin['nsources']
  gamma = Xmin['gamma']

  msg  = "   - % 16s ra %6.2f dec %6.2f ns %4.1f gamma %.1f" % tuple(src)
  msg += " (7yr paper) | ns %4.1f gamma %.1f (this code)" % (ns, gamma)
  print(msg)

  results['TS'][i] = TS

###########
# RESULTS #
###########

###############################################################################

############
# P-VALUES #
############

# estimate p-values using 5k background scrambles at each location

print("\n p-values:")
for result in results:

  t1 = time.time()
  trials = llh.do_trials(src_ra  = np.radians(result['ra_deg']),
                         src_dec = np.radians(result['dec_deg']),
                         n_iter  = 5000)
  dt = time.time()-t1

  p = trials[ trials['TS'] >= result['TS'] ].size / float(trials.size)

  if p > 0:
    vals = (result['name'], result['TS'], p, dt)
    print("   - % 16s ts %5.2f p %.3f (%6.2f sec)" % vals)
  else:
    vals = (result['name'], result['TS'], 1/float(trials.size), dt)
    print("   - % 16s ts %5.2f p < %.1e (%6.2f sec)" % vals)

  np.save(result['name'].replace(" ","_") + ".npy", trials)

print("\n Completed in %.2f sec" % (time.time()-t0))

############
# P-VALUES #
############

###############################################################################
