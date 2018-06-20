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

nch40_sim, nch40_data = cache.load (filename_pickle+"IC40/NCh.pickle")



#Those sigmas up there are w/o pull correction - the following is the fix.
mc40 = cache.load(filename_pickle+'IC40/mc.pickle')
exp40 = cache.load(filename_pickle+'IC40/exp.pickle')
mc59 = cache.load(filename_pickle+'IC59/mc.pickle')
exp59 = cache.load(filename_pickle+'IC59/exp.pickle')
mc79 = cache.load(filename_pickle+'IC79/mc.pickle')
exp79 = cache.load(filename_pickle+'IC79/exp.pickle')
mc86I = cache.load(filename_pickle+'IC86I/mc.pickle')
exp86I = cache.load(filename_pickle+'IC86I/exp.pickle')
mc3yr = cache.load(filename_pickle+'3yr/mc.pickle')
exp3yr = cache.load(filename_pickle+'3yr/exp.pickle')

#Adding in a few params from epinat's 3yr
muex3yr_data = 10**(exp3yr['logE'])
muex3yr_sim = 10**(mc3yr['logE'])
sindec3yr_data = exp3yr['sinDec']
sindec3yr_sim = mc3yr['sinDec']


sigma3yr_data = exp3yr['sigma']
sigma3yr_sim = mc3yr['sigma']
sigma86I_data = exp86I['sigma']
sigma86I_sim = mc86I['sigma']
sigma79_data = exp79['sigma']
sigma79_sim = mc79['sigma']
sigma59_data = exp59['sigma']
sigma59_sim = mc59['sigma']
sigma40_data = exp40['sigma']
sigma40_sim = mc40['sigma']
true40 = mc40['trueE']
true59 = mc59['trueE']
true79 = mc79['trueE']
true86I = mc86I['trueE']
true3yr = mc3yr['trueE']

weight_3yr = mc3yr['ow']*true3yr**(-2)
weight_86 = mc86I['ow']*true86I**(-2)
weight_79 = mc79['ow']*true79**(-2)
weight_59 = mc59['ow']*true59**(-2)
weight_40 = mc40['ow']*true40**(-2)


bins = 80
bins2d = (80,80)
### Energy ###
energyrange = (10,1e7)
fig_energyhist = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h_energy_3yr = histlite.hist(muex3yr_data, bins=bins, range = energyrange, log=True)
h_energy_86 = histlite.hist(muex86I_data, bins=bins, range = energyrange, log=True)
h_energy_79 = histlite.hist(muex79_data, bins=bins, range = energyrange, log=True)
h_energy_59 = histlite.hist(mue59_data, bins=bins, range = energyrange, log=True)
h_energy_40 = histlite.hist(mue40_data, bins=bins, range = energyrange, log=True)
histlite.plot1d(ax,h_energy_40,histtype='step',label='IC40', color = 'blue')
histlite.plot1d(ax,h_energy_59,histtype='step',label='IC59', color = 'green')
histlite.plot1d(ax,h_energy_79,histtype='step',label='IC79', color = 'red')
histlite.plot1d(ax,h_energy_86,histtype='step',label='IC86I', color = 'cyan')
histlite.plot1d(ax,h_energy_3yr,histtype='step',label='IC86II-IV', color = 'magenta')

#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"Energy Proxy (GeV)")
ax.set_ylabel("Counts per bin")
ax.loglog()
ax.set_ylim(1,1e5)
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propsmall)
fig_energyhist.savefig(plotfolder + 'energy_hists_7yr.png')

###Declination###
decrange = (-1.,1.)
fig_dechist = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

h_dec_3yr = histlite.hist(sindec3yr_data, bins=bins, range=decrange)
h_dec_86 = histlite.hist(sindec86I_data, bins=bins, range=decrange)
h_dec_79 = histlite.hist(sindec79_data, bins=bins, range=decrange)
h_dec_59 = histlite.hist(sindec59_data, bins=bins, range=decrange)
h_dec_40 = histlite.hist(sindec40_data, bins=bins, range=decrange)

histlite.plot1d(ax,h_dec_40,histtype='step',label='IC40', color = 'blue')
histlite.plot1d(ax,h_dec_59,histtype='step',label='IC59', color = 'green')
histlite.plot1d(ax,h_dec_79,histtype='step',label='IC79', color = 'red')
histlite.plot1d(ax,h_dec_86,histtype='step',label='IC86I', color = 'cyan')
histlite.plot1d(ax,h_dec_3yr,histtype='step',label='IC86II-IV', color = 'magenta')

#ax.set_title("4yr PS data - sindec")
ax.set_xlabel(r"$\sin({\delta})$")
ax.set_ylabel("Counts per bin")
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_xlim(-1.,1.)
ax.set_ylim(0,7000)
ax.legend(loc="best", prop=propsmall)
fig_dechist.savefig(plotfolder + 'dec_hists_7yr.png')


### Sigma ###
#I think this is before pull correction.
sigmarange = (1e-2,1e3)
fig_sigma = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h_3yr = histlite.Hist.normalize(histlite.hist(np.degrees(sigma3yr_data), bins=bins, range = sigmarange, log=True))
h_86 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma86I_data), bins=bins, range = sigmarange, log=True))
h_79 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma79_data), bins=bins, range = sigmarange, log=True))
h_59 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma59_data), bins=bins, range = sigmarange, log=True))
h_40 = histlite.Hist.normalize(histlite.hist(np.degrees(sigma40_data), bins=bins, range = sigmarange, log=True))
histlite.plot1d(ax,h_40,histtype='step',label='IC40', color = 'blue')
histlite.plot1d(ax,h_59,histtype='step',label='IC59', color = 'green')
histlite.plot1d(ax,h_79,histtype='step',label='IC79', color = 'red')
histlite.plot1d(ax,h_86,histtype='step',label='IC86I', color = 'cyan')
histlite.plot1d(ax,h_3yr,histtype='step',label='IC86II-IV', color = 'magenta')

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
ax.set_ylim(1e-5,1)
ax.legend(loc="best", prop=propsmall)
fig_sigma.savefig(plotfolder + 'sigma_hists_7yr.png')


