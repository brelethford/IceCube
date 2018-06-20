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

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

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
plotfolder = '/data/user/brelethford/AGN_Core/Plots/data/'

ra86I_sim, sindec86I_sim , ra86I_data, sindec86I_data, ra86I_true, sindec86I_true, energy86I_true, muex86I_sim, muex86I_data, sigma86I_sim, sigma86I_data, OneWeight_IC86I, dpsi_IC86I = cache.load (filename_pickle+"IC86I/coords.pickle")

ra79_sim, sindec79_sim , ra79_data, sindec79_data, ra79_true, sindec79_true, energy79_true, muex79_sim, muex79_data, sigma79_sim, sigma79_data, OneWeight_IC79, dpsi_IC79 = cache.load (filename_pickle+"IC79/coords.pickle")

nch79_sim, nch79_data = cache.load (filename_pickle+"IC79/NCh.pickle")

ra59_sim, sindec59_sim , ra59_data, sindec59_data, ra59_true, sindec59_true, energy59_true, mue59_sim, mue59_data, sigma59_sim, sigma59_data, OneWeight_IC59, dpsi_IC59 = cache.load (filename_pickle+"IC59/coords.pickle")

nch59_sim, nch59_data = cache.load (filename_pickle+"IC59/NCh.pickle")

ra40_sim, sindec40_sim , ra40_data, sindec40_data, ra40_true, sindec40_true, energy40_true, mue40_sim, mue40_data, sigma40_sim, sigma40_data, OneWeight_IC40, dpsi_IC40 = cache.load (filename_pickle+"IC40/coords.pickle")

nch40_sim, nch40_data = cache.load (filename_pickle+"IC59/NCh.pickle")

bins = 80

### Energy ###

fig_energyhist = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h_energy_86 = histlite.Hist.normalize(histlite.hist(muex86I_data, bins=bins, log=True))
h_energy_79 = histlite.Hist.normalize(histlite.hist(muex79_data, bins=bins, log=True))
h_energy_59 = histlite.Hist.normalize(histlite.hist(mue59_data, bins=bins, log=True))
h_energy_40 = histlite.Hist.normalize(histlite.hist(mue40_data, bins=bins, log=True))
histlite.plot1d(ax,h_energy_79,histtype='step',label='IC79', color = 'blue')
histlite.plot1d(ax,h_energy_86,histtype='step',label='IC86', color = 'red')
histlite.plot1d(ax,h_energy_59,histtype='step',label='IC59', color = 'green')
histlite.plot1d(ax,h_energy_40,histtype='step',label='IC40', color = 'purple')

ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"logE (GeV)")
ax.set_ylabel("Relative Abundance")
ax.loglog()
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propsmall)
fig_energyhist.savefig(plotfolder + '/energy_hists_normalized.pdf')

###Declination###

fig_dechist = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

h_dec_86 = histlite.Hist.normalize(histlite.hist(sindec86I_data, bins=bins))
h_dec_79 = histlite.Hist.normalize(histlite.hist(sindec79_data, bins=bins))
h_dec_59 = histlite.Hist.normalize(histlite.hist(sindec59_data, bins=bins))
h_dec_40 = histlite.Hist.normalize(histlite.hist(sindec40_data, bins=bins))

histlite.plot1d(ax,h_dec_79,histtype='step',label='IC79', color = 'blue')
histlite.plot1d(ax,h_dec_86,histtype='step',label='IC86', color = 'red')
histlite.plot1d(ax,h_dec_59,histtype='step',label='IC59', color = 'green')
histlite.plot1d(ax,h_dec_40,histtype='step',label='IC40', color = 'purple')

ax.set_title("4yr PS data - sindec")
ax.set_xlabel(r"$\sin({\delta})$")
ax.set_ylabel("Relative Abundance")
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_xlim(-1.,1.)
ax.legend(loc="best", prop=propsmall)
fig_dechist.savefig(plotfolder + '/dec_hists_normalized.pdf')



