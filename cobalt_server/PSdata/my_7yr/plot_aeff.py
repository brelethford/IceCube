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
from icecube.umdtools import cache, misc, arrays
from icecube import icetray, dataclasses, histlite, astro
from scipy.stats import chi2
import healpy as hp
from scipy.signal import convolve2d
#The purpose of this script is to compare the dspi and rescaled sigma.

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/PSdata/mike_plots")
#This one is for colormaps.py

misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size = 'small')
propxsmall = mpl.font_manager.FontProperties (size = 'x-small')

#Loading Zone#
projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/AGN_Core_Sample/'
filename_pickle = datafolder+'pickle/'
plotfolder = '/data/user/brelethford/AGN_Core/Plots/data/'
'''
ra86I_sim, sindec86I_sim , ra86I_data, sindec86I_data, ra86I_true, sindec86I_true, energy86I_true, muex86I_sim, muex86I_data, sigma86I_sim, sigma86I_data, OneWeight_IC86I, dpsi_IC86I = cache.load (filename_pickle+"IC86I/coords.pickle")

ra79_sim, sindec79_sim , ra79_data, sindec79_data, ra79_true, sindec79_true, energy79_true, muex79_sim, muex79_data, sigma79_sim, sigma79_data, OneWeight_IC79, dpsi_IC79 = cache.load (filename_pickle+"IC79/coords.pickle")

nch79_sim, nch79_data = cache.load (filename_pickle+"IC79/NCh.pickle")

ra59_sim, sindec59_sim , ra59_data, sindec59_data, ra59_true, sindec59_true, energy59_true, mue59_sim, mue59_data, sigma59_sim, sigma59_data, OneWeight_IC59, dpsi_IC59 = cache.load (filename_pickle+"IC59/coords.pickle")

nch59_sim, nch59_data = cache.load (filename_pickle+"IC59/NCh.pickle")

ra40_sim, sindec40_sim , ra40_data, sindec40_data, ra40_true, sindec40_true, energy40_true, mue40_sim, mue40_data, sigma40_sim, sigma40_data, OneWeight_IC40, dpsi_IC40 = cache.load (filename_pickle+"IC40/coords.pickle")

nch40_sim, nch40_data = cache.load (filename_pickle+"IC59/NCh.pickle")
'''
ngen_40 = cache.load(filename_pickle+"IC40/n_gen.pickle")
ngen_59 = cache.load(filename_pickle+"IC59/n_gen.pickle")
ngen_79 = cache.load(filename_pickle+"IC79/n_gen.pickle")
ngen_86I = cache.load(filename_pickle+"IC86I/n_gen.pickle")

bins = 80

#To do effective area calculations, I need nugen data. I can get this from my MC data... I think. Following:
###http://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/richman/scripts/trackps_ic861/trunk/tps861.py#L618###


import data_multi

def getfig (fignum=None, aspect=None, width=None, figsize=None):
    aspect = aspect or 4/3.
    width = width or 7
    if figsize is None:
        figsize = (width, width / aspect)
    out = plt.figure (num=fignum, figsize=figsize)
    plt.clf ()
    return out

###Effective Area###
def aeff2d (year):
   """
   Plot and tabulate Aeff.
   """

   import colormaps as cmaps
   plt.register_cmap (name='viridis', cmap=cmaps.viridis)
   plt.set_cmap (cmaps.viridis)

   logEmin, logEmax = 2., 9.
   dlogE = 0.1
   n_bins_E = (logEmax - logEmin) / dlogE
   dcz = 0.01
   dOmega = 2 * np.pi * dcz
   n_bins_cz = 2 / dcz
   #nu = self.nu ..... xn is the mc sample of a certain year. let's start with IC86I and work our way forward.
   if year == 86:
     xn = data_multi.MC86I()
   elif year ==79:
     xn = data_multi.MC79()
   elif year ==59:
     xn = data_multi.MC59()
   elif year ==40:
     xn = data_multi.MC40()

   #Whichever year we have, the following occurs:
   nu = arrays.Arrays (dict (
          (k,xn[k])
          for k in ('ra', 'sinDec', 'sigma', 'logE',
                  'trueRa', 'trueDec', 'trueE', 'ow')))
   #nu = self.nu
   nu.cz = -np.sin (nu.trueDec)
   w_aeff = 1 / (1e4 * np.log (10)) * nu.ow / nu.trueE / dOmega / dlogE

   h_aeff = histlite.hist (
       (nu.trueE, nu.cz), w_aeff,
       bins=(n_bins_E, n_bins_cz),
       range=((10**logEmin, 10**logEmax), (-1, 1)),
       log=(True, False),
   )

   fig = getfig (aspect=4/3., width=5)
   ax = fig.add_subplot (111)
   fig.subplots_adjust (bottom=.15, left=.15)
   result = histlite.plot2d (ax, h_aeff, cbar=True, log=True,
                       vmin=5e-6, vmax=1e4, zmin=5e-6)
   result['colorbar'].set_label (r'effective area $[\textrm{m}^2]$')
   ax.set_title('Aeff - IC{}'.format(str(year)))
   ax.set_xlabel ('neutrino energy [GeV]')
   ax.set_ylabel (r'$\cos(\textrm{zenith})$')

   plot_dir = misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/data')
   fig.savefig(plot_dir+'/aeff{}_2d.pdf'.format(str(year)))

def aeff1d (year):
   """
   Plot and tabulate Aeff.
   """

   #nu = self.nu ..... xn is the mc sample of a certain year. let's start with IC86I and work our way forward.
   if year == 86:
     full_xn = data_multi.MC86I()
     ngen = ngen_86I
   elif year ==79:
     full_xn = data_multi.MC79()
     ngen = ngen_79
   elif year ==59:
     full_xn = data_multi.MC59()
     ngen = ngen_59
   elif year ==40:
     full_xn = data_multi.MC40()
     ngen = ngen_40
   
   logEmin, logEmax = 1., 8.
   dlogE = 0.1
   n_bins_E = (logEmax - logEmin) / dlogE

   #I'll use the same curves as the IC86I page: https://wiki.icecube.wisc.edu/index.php/IC86_I_Point_Source_Analysis/Performance_Plots
   #Now we gotta do each of a number of ranges to match existing plots.
   decedges = np.array([-90.,-60.,-30.,0.,30.,60.,90.])
   
   aeff_hists = []

   for (dec_min, dec_max) in zip (decedges[:-1],decedges[1:]):
     dsd = np.sin(np.radians(dec_max - dec_min))
     dOmega = 2 * np.pi * dsd
     n_bins_sd = 2 / dsd
     print(dec_min,dec_max) 
     #first, define xn from the full set, before we pare it down by sindec.
     xn = full_xn    
     #Whichever year we have, the following occurs:
     mask = ((np.sin(np.radians(dec_min)) < np.sin(xn['trueDec'])) & (np.sin(xn['trueDec'])< np.sin(np.radians(dec_max))))
     xn = xn[mask]
     nu = arrays.Arrays (dict (
          (k,xn[k])
          for k in ('ra', 'sinDec', 'sigma', 'logE',
                  'trueRa', 'trueDec', 'trueE', 'ow')))
     #nu = self.nu
     #nu.cz = -np.sin (nu.trueDec)
     nu.sd = np.sin (nu.trueDec)
     w_aeff = 1 / (1e4 * np.log (10)) * nu.ow / nu.trueE / dOmega / dlogE
     h_aeff = histlite.hist (
       (nu.trueE), w_aeff,
       bins=(n_bins_E),
       range=(10**logEmin, 10**logEmax),
       log=(True)
     )
     aeff_hists.append(h_aeff)
   #that should give us all the hists we want. Now let's plot them on one axis.
   fig = plt.figure(figsize=(w, .75*w))
   ax = fig.add_subplot (111)
   fig.subplots_adjust (bottom=.15, left=.18)
   colors = ['red', 'green', 'orange', 'purple', 'blue', 'black']
   for i in range(len(decedges)-1):
     histlite.plot1d (ax, aeff_hists[i], log=True, color=colors[i],
                       vmin=5e-6, vmax=1e4, label =  r'${0:.1f} < \delta < {1:.1f}$'.format(decedges[i], decedges[i+1]))
   ax.set_title('Aeff - IC{}'.format(str(year)))
   ax.set_xlabel ('neutrino energy [GeV]')
   ax.set_ylabel (r'Effective Area $[\textrm{m}^2]$')
   ax.legend(loc='best',prop=propxsmall)
   ax.loglog()
   plot_dir = misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/data')
   fig.savefig(plot_dir+'/aeff{}_1d.pdf'.format(str(year)))



aeff2d(86)
aeff2d(79)
aeff2d(59)
aeff2d(40)

aeff1d(86)
aeff1d(79)
aeff1d(59)
aeff1d(40)

