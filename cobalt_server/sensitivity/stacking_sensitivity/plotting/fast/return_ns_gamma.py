import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from skylab import utils
from scipy.stats import chi2, norm, exponweib
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite
from optparse import OptionParser
import argparse
fitfun = utils.FitDeltaChi2()

def getfig (fignum=None, aspect=None, width=None, figsize=None):
    aspect = aspect or 4/3.
    width = width or 7
    if figsize is None:
        figsize = (width, width / aspect)
    out = plt.figure (num=fignum, figsize=figsize)
    plt.clf ()
    return out

def pfig (*a, **kw):
    fig = getfig (*a, **kw)
    ax = fig.add_subplot (111)
    return fig, ax

##This script should be able to plot the fit ns and gamma for bg only.

##arguments:

parser = OptionParser (usage = '%prog [options]')


parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--cat', dest = 'cat', type = str,
                help = 'Source catalog')

parser.add_option ('--weight', dest = 'weight', type = str,
                help = 'Source weight')

opts, args = parser.parse_args ()
years = opts.n
cat = opts.cat

if cat == 'SwiftBAT70m':
        weight = opts.weight
elif cat == '4yr_Starburst':
        weight = 'S60m'
elif cat == 'teststack' or cat == 'teststack50' or cat == 'teststack300':
        weight  = 'equal'
elif cat == 'blackhole':
        weight = 'flux2m'
elif cat == 'Milagro17':
        weight = 'weight'

##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

sigdata = '/data/user/brelethford/Output/fast/stacking_sensitivity/{0}yr/{1}/{2}/sig/'.format(str(years),cat,weight)

sig_dict = cache.load(sigdata+'inj_trials.array')

mu = []
ns_recover = []
gamma_recover = []

for mu_i in sig_dict.viewkeys():
  mu.append(mu_i)
  trials = sig_dict[mu_i]
  n_return = trials['n_inj']
  gamma_return = trials['gamma']
  ns_recover.append(n_return)
  gamma_recover.append(gamma_return)

  fig,ax = pfig()
  plt.scatter(n_return, gamma_return, alpha=0.5, label = "{} injected events".format(str(mu_i)),color = 'blue')
  plt.ylim(1,4)
  plt.legend(loc='upper right', prop=propsmall)
  plt.xlabel('ns')
  plt.ylabel(r'$\gamma$')
  plt.title('ns-gamma fit for {} injected events'.format(mu_i))

  print ("plotting n_inj = {}".format(mu_i))
  misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/{}/inj_test/'.format(str(years),cat,weight))
  fig.savefig('/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/{}/inj_test/{}_inj_ns_gamma_fit.png'.format(str(years),cat,weight,str(mu_i)))
  misc.ensure_dir('/home/brelethford/public_html/skylab/fast/plots/stacking/{}yr/{}/{}/inj_test/'.format(str(years),cat,weight))
  fig.savefig('/home/brelethford/public_html/skylab/fast/plots/stacking/{}yr/{}/{}/inj_test/{}_inj_ns_gamma_fit.png'.format(str(years),cat,weight,str(mu_i)))
  plt.close()

#Now plot the whole thing.
#First ns#

fig,ax = pfig()
#plt.scatter(mu, [np.median(ns) for ns in ns_recover], alpha=0.5, label = "{} trials per point".format(len(ns_recover[0])),color = 'blue')
plt.errorbar(mu,[np.median(ns) for ns in ns_recover], alpha = 0.5, yerr = [np.std(ns) for ns in ns_recover], fmt = 'o', label = "{} trials per point".format(len(ns_recover[0])),color = 'blue')
plt.legend(loc='upper right', prop=propsmall)
plt.xlabel('ns injected')
plt.ylabel('ns recovered')
if cat == '4yr_Starburst':
  plt.title('ns recovered for Starburst')
else:
  plt.title('ns recovered for {}'.format(str(cat)))
print ("plotting summary:")
fig.savefig('/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/{}/inj_test/ns_recovered.png'.format(str(years),cat,weight))
fig.savefig('/home/brelethford/public_html/skylab/fast/plots/stacking/{}yr/{}/{}/inj_test/ns_recovered.png'.format(str(years),cat,weight))
plt.close()
print("Complete!")

# Now  ns#

fig,ax = pfig()
plt.errorbar(mu,[np.median(gamma) for gamma in gamma_recover], alpha = 0.5, yerr = [np.std(gamma) for gamma in gamma_recover], fmt = 'o', label = "{} trials per point".format(len(ns_recover[0])),color = 'blue')
plt.legend(loc='upper right', prop=propsmall)
plt.xlabel('ns injected')
plt.ylabel('median gamma recovered')
if cat == '4yr_Starburst':
  plt.title('gamma recovered for Starburst')
else:
  plt.title('gamma recovered for {}'.format(str(cat)))
print ("plotting summary:")
fig.savefig('/data/user/brelethford/AGN_Core/Plots/fast/plots/stacking/{}yr/{}/{}/inj_test/gamma_recovered.png'.format(str(years),cat,weight))
fig.savefig('/home/brelethford/public_html/skylab/fast/plots/stacking/{}yr/{}/{}/inj_test/gamma_recovered.png'.format(str(years),cat,weight))
plt.close()
print("Complete!")
