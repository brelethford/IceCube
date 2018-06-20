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

dpsi3yr = cache.load(filename_pickle+'epinat_3yr/dpsi.pickle')
#Those sigmas up there are w/o pull correction - the following is the fix.
mc40 = cache.load(filename_pickle+'IC40/mc.pickle')
exp40 = cache.load(filename_pickle+'IC40/exp.pickle')
mc59 = cache.load(filename_pickle+'IC59/mc.pickle')
exp59 = cache.load(filename_pickle+'IC59/exp.pickle')
mc79 = cache.load(filename_pickle+'IC79/mc.pickle')
exp79 = cache.load(filename_pickle+'IC79/exp.pickle')
mc86I = cache.load(filename_pickle+'IC86I/mc.pickle')
exp86I = cache.load(filename_pickle+'IC86I/exp.pickle')
exp3yr = cache.load(filename_pickle+'epinat_3yr/exp.pickle')
mc3yr = cache.load(filename_pickle+'epinat_3yr/mc.pickle')

sigma3yr_data = exp3yr['sigma']
sigma86I_data = exp86I['sigma']
sigma79_data = exp79['sigma']
sigma59_data = exp59['sigma']
sigma40_data = exp40['sigma']
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

weight_86 = OneWeight_IC86I*energy86I_true**(-2)
weight_79 = OneWeight_IC79*energy79_true**(-2)
weight_59 = OneWeight_IC59*energy59_true**(-2)
weight_40 = OneWeight_IC40*energy40_true**(-2)


bins = 80
bins2d = (28,5e3)

### Energy x sigma###
energyrange = (100,1e9)
sigmarange = (0,180)
fig_angular_res = plt.figure (figsize=(w, .75*w)) 
ax=plt.gca()

h_3yr = histlite.hist((true3yr,np.degrees(dpsi3yr['dpsi'])), bins=bins2d, weights = weight_3yr,  range = (energyrange,sigmarange), log=(True,False))
h_86 = histlite.hist((energy86I_true,np.degrees(dpsi_IC86I)), bins=bins2d, weights = weight_86,  range = (energyrange,sigmarange), log=(True,False))
h_79 = histlite.hist((energy79_true,np.degrees(dpsi_IC79)), bins=bins2d,  weights = weight_79, range = (energyrange,sigmarange), log=(True,False))
h_59 = histlite.hist((energy59_true,np.degrees(dpsi_IC59)), bins=bins2d,  weights = weight_59, range = (energyrange,sigmarange), log=(True,False))
h_40 = histlite.hist((energy40_true,np.degrees(dpsi_IC40)), bins=bins2d,  weights = weight_40, range = (energyrange,sigmarange), log=(True,False))


histlite.plot1d(ax,h_40.median(1),drawstyle='lines',label='IC40', color = 'blue')
histlite.plot1d(ax,h_59.median(1),drawstyle='lines',label='IC59', color = 'green')
histlite.plot1d(ax,h_79.median(1),drawstyle='lines',label='IC79', color = 'red')
histlite.plot1d(ax,h_86.median(1),drawstyle='lines',label='IC86', color = 'cyan')
histlite.plot1d(ax,h_3yr.median(1),drawstyle='lines',label='IC86II-IC86IV', color = 'magenta')

#ax.set_title("4yr PS data - energy")
ax.set_xlabel(r"Energy Proxy (GeV)")
ax.set_ylabel("Median angular resolution")
ax.semilogx()
ax.set_ylim(0,2)
ax.set_xlim(100,1e9)
plt.subplots_adjust (left=.2, bottom=.2)
ax.legend(loc="best", prop=propsmall)
fig_angular_res.savefig(plotfolder + 'angular_res_hists.pdf')
fig_angular_res.savefig('/home/brelethford/public_html/wiki_images/angular_res_hists.png')

