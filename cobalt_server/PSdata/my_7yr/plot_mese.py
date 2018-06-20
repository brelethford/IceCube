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
filename_pickle = datafolder+'pickle/MESE/'
plotfolder = '/data/user/brelethford/AGN_Core/Plots/data/'


#Those sigmas up there are w/o pull correction - the following is the fix.
mc = cache.load(filename_pickle+'mc.pickle')
exp = cache.load(filename_pickle+'exp.pickle')

#Adding in a few params from epinat's 3yr
muex_data = 10**(exp['logE'])
true = mc['trueE']
muex_sim = 10**(mc['logE'])
sindec_data = exp['sinDec']
sindec_sim = mc['sinDec']


sigma_data = exp['sigma']
sigma_sim = mc['sigma']

weight = mc['ow']*true**(-2)


bins = 80
bins2d = (80,80)
### Energy ###
energyrange = (1e3,200000)
fig_energyhist = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h_energy = histlite.hist(muex_data, bins=bins, range = energyrange, log=True)
histlite.plot1d(ax,h_energy,histtype='step',label='MESE', color = 'k')

#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"Energy Proxy (GeV)")
ax.set_ylabel("Counts per bin")
ax.loglog()
ax.set_ylim(1,1e2)
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propsmall)
fig_energyhist.savefig(plotfolder + 'energy_hists_MESE.png')

###Declination###
decrange = (-1.,1.)
fig_dechist = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

h_dec = histlite.hist(sindec_data, bins=bins, range=decrange)

histlite.plot1d(ax,h_dec,histtype='step',label='MESE', color = 'k')

#ax.set_title("4yr PS data - sindec")
ax.set_xlabel(r"$\sin({\delta})$")
ax.set_ylabel("Counts per bin")
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_xlim(-1.,0.)
ax.set_ylim(0,100)
ax.legend(loc="best", prop=propsmall)
fig_dechist.savefig(plotfolder + 'dec_hists_MESE.png')

### Sigma ###
#I think this is before pull correction.
sigmarange = (1e-1,1e2)
fig_sigma = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h = histlite.Hist.normalize(histlite.hist(np.degrees(sigma_data), bins=bins, range = sigmarange, log=True))
histlite.plot1d(ax,h,histtype='step',label='MESE', color = 'k')

#h_3yr = histlite.Hist.normalize(histlite.hist(np.degrees(sigma3yr_sim), weights = weight_3yr, bins=bins,range = sigmarange, log=True))
#h_86 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma86I_sim), weights = weight_86, bins=bins,range = sigmarange, log=True))
#h_79 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma79_sim), weights = weight_79, bins=bins,range = sigmarange, log=True))
#h_59 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma59_sim), weights = weight_59, bins=bins,range = sigmarange, log=True))
#h_40 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma40_sim), weights = weight_40, bins=bins,range = sigmarange, log=True))
#histlite.plot1d(ax,h_40,histtype='step', color = 'blue', linestyle = ':')
#histlite.plot1d(ax,h_59,histtype='step', color = 'green', linestyle = ':')
#histlite.plot1d(ax,h_79,histtype='step', color = 'red', linestyle = ':')
#histlite.plot1d(ax,h_86,histtype='step', color = 'cyan', linestyle = ':')
#histlite.plot1d(ax,h_3yr,histtype='step', color = 'magenta', linestyle = ':')
#print ('number of exp: '+str(len(sigma40_data)))
#print ('number of mc: '+str(len(sigma40_sim)))
#print ('number of exp: '+str(len(sigma59_data)))
#print ('number of mc: '+str(len(sigma59_sim)))
#print ('number of exp: '+str(len(sigma79_data)))
#print ('number of mc: '+str(len(sigma79_sim)))
#print ('number of exp: '+str(len(sigma86I_data)))
#print ('number of mc: '+str(len(sigma86I_sim)))

#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r'$\sigma$ [$^\circ$]')
ax.set_ylabel("Probability Density")
ax.loglog()
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylim(1e-3,1)
ax.legend(loc="best", prop=propsmall)
fig_sigma.savefig(plotfolder + 'sigma_hists_MESE.png')


