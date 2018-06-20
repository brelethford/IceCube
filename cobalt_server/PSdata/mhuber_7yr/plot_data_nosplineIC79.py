import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sp
import tables
import os
from skylab.psLLH import PointSourceLLH
from skylab.ps_model import ClassicLLH  
from skylab.ps_injector import PointSourceInjector
from skylab.psLLH import MultiPointSourceLLH
from skylab.utils import poisson_weight
from icecube.umdtools import cache, misc
from icecube import icetray, dataclasses, histlite, astro
from scipy.stats import chi2
import healpy as hp
from scipy.signal import convolve2d
#The purpose of this script is to compare the dspi and rescaled sigma.
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/")
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload/")
import load_mstacking 
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size = 'small')
propxsmall = mpl.font_manager.FontProperties (size = 'x-small')

#Livetimes (in days)
livetime_IC86I = 332.61
livetime_IC79  = 315.506
livetime_IC59  = 348.138

#Loading Zone#
projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/AGN_Core_Sample/'
filename_pickle = datafolder+'pickle/'
plotfolder = '/data/user/brelethford/AGN_Core/Plots/data/mhuber_7yr/nospline_79/'


########## mhuber's 7yr sample #########

llh40 = load_mstacking.ic40()
llh59 = load_mstacking.ic59()
llh79 = load_mstacking.ic79()
llh86I = load_mstacking.ic86_I()
llh86II = load_mstacking.ic86_2012()

llh = [llh40,llh59,llh79,llh86I,llh86II]
exp = []
for l in llh:
  exp.append(l.exp)
##Also have to load all at once in order to get the mc in a good format. Probably can do this separately, but this way ensures I use the exact right MC.

llh7yr = load_mstacking.load_7yr_nospline_79()
mc = load_mstacking.monte_carlo(llh7yr)

##Now I have exp and mc for all 7 yrs. plot. The framework for plotting is already set up, so let's just make sure the info is in a manner it's used to.

bins = 80

### Energy ###

fig_energyhist = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()
labels = ['IC40','IC59','IC79','IC86I','IC86II-IV']
colors = ['purple','green','blue','red','orange']

for i,l,c in zip(exp,labels,colors):
  h = histlite.hist(10**(i['logE']), bins=bins, log=True)
  histlite.plot1d(ax,h,histtype='step',label=l, color = c, range = (10,1e7))
#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"Energy Proxy (GeV)")
ax.set_ylabel("Counts per bin")
ax.loglog()
ax.set_xlim(10,1e7)
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propxsmall)
fig_energyhist.savefig(plotfolder + 'energy_hists_7yr.pdf')

###Declination###

fig_dechist = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

for i,l,c in zip(exp,labels,colors):
  h = histlite.hist(i['sinDec'], bins=bins)
  histlite.plot1d(ax,h,histtype='step',label=l, color = c, range = (-1.1))

#ax.set_title("4yr PS data - sindec")
ax.set_xlabel(r"$\sin({\delta})$")
ax.set_ylabel("Counts per bin")
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_xlim(-1.,1.)
ax.legend(loc="best", prop=propsmall)
fig_dechist.savefig(plotfolder + 'dec_hists_7yr.pdf')

### Sigma ###
##need weights for mc for this one.

fig_sigma = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()
for m,i,l,c in zip(mc,exp,labels,colors):
  weight_mc = mc[m]['ow']*mc[m]['trueE']**(-2)
  hexp = histlite.hist(np.degrees(i['sigma']), bins=bins, log=True)
  histlite.plot1d(ax,hexp,histtype='step',label=l, color = c, range = (1e-2,1e3))
  hmc = histlite.hist(np.degrees(mc[m]['sigma']), bins=bins, weights = weight_mc, log=True)
  histlite.plot1d(ax,hmc,histtype='step', color = c, range = (1e-2,1e3),linestyle = ':')
  #print ('number of exp: '+ str(len(i['sigma'])))
  #print ('number of mc: '+ str(len(mc[m]['sigma'])))

#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"Sigma (deg)")
ax.set_ylabel("Counts per bin")
ax.loglog()
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propxsmall)
ax.set_xlim(1e-2,1e3)
fig_sigma.savefig(plotfolder + 'sigma_hists_7yr.pdf')

