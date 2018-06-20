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

opts, args = parser.parse_args ()
years = opts.n

##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

alldata = '/data/user/brelethford/Output/standard/allsky_sensitivity/{0}yr/'.format(str(years))

bckg = alldata + 'background_trials/'

decbins=35
decs = np.linspace(-85.,85.,decbins)
decfolders = [bckg+'dec_{0:4f}/'.format(dec) for dec in decs]


def ns_gamma_fit(datafolder):
  files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
  TS=[]
  ns=[]
  gamma = []
  for file in files:
    for entry in file:
      ns.append(entry[2])
      gamma.append(entry[3])
      if entry[0]==0:
        if entry[1]<0:
          TS.append(0)
        else:
          TS.append(entry[1])
  TSs=TS
  chi2fit_flux= fitfun.fit(TSs)
  print('background only: median TS = {}'.format(str(chi2fit_flux.isf(0.5))))
  print ('max TS = {}'.format(str(max(TSs))))
  print ('Number of trials is: {}'.format(str(len(TSs))))
  print ('Percentage of TS=0 is: {}'.format(str(1.0-np.float(np.count_nonzero(TSs))/np.float(len(TSs)))))
  eta = chi2fit_flux.eta
  df = chi2fit_flux.params[0]
  
  
  return TSs, ns, gamma, eta, df

TSs, nss, gammas, etas, ndofs = np.array([ns_gamma_fit(dec) for dec in decfolders]).T

## Now to plot. ##

for dec,ns,gamma in zip(decs,nss,gammas):
  fig,ax = pfig()
  plt.scatter(ns, gamma, alpha=0.5, label = "{} bg trials".format(str(len(ns))),color = 'blue')
  plt.ylim(1,4)
  plt.xlim(0.,10.)
  plt.legend(loc='upper right', prop=propsmall)
  plt.xlabel('ns')
  plt.ylabel(r'$\gamma$')
  plt.title('ns-gamma fit')

  misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/allsky/bg_tsds/{}yr/{}/'.format(str(years),str(dec)))
  fig.savefig ('/data/user/brelethford/AGN_Core/Plots/allsky/bg_tsds/{}yr/{}/ns_gamma_fit.png'.format(str(years),str(dec)))
  misc.ensure_dir('/home/brelethford/public_html/skylab/allsky/bg_tsds/{}yr/{}/'.format(str(years),str(dec)))
  fig.savefig('/home/brelethford/public_html/skylab/allsky/bg_tsds/{}yr/{}/ns_gamma_fit.png'.format(str(years),str(dec)))
  plt.close()

