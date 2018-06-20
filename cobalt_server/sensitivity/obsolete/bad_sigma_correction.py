import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
###This script imports the pickled data from stefan, as well as my data with correct and incorrect pull correction.###

datafolder = '/data/user/brelethford/Data/AGN_Core_Sample/pickle/'
stefanfolder = '/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/stefan_c/stefan_data/'

ben_exp = cache.load (datafolder + 'exp.pickle')
ben_mc = cache.load (datafolder + 'MC.pickle')

bad_exp = cache.load (datafolder + 'exp_bad_sigma.pickle')
bad_mc = cache.load (datafolder + 'MC_bad_sigma.pickle')

stefan_exp = np.genfromtxt(stefanfolder+'IC86-I_data.txt',names=True)
stefan_mc = np.genfromtxt(stefanfolder+'IC86-I_MC.txt',names=True)

#extract the sigmas

ben_exp_sigma =ben_exp['sigma']
ben_mc_sigma =ben_mc['sigma']

ben_hist_exp = histlite.hist(ben_exp_sigma,bins=100,log=True)
ben_hist_mc = histlite.hist(ben_mc_sigma,bins=100,log=True)

bad_exp_sigma =bad_exp['sigma']
bad_mc_sigma =bad_mc['sigma']

bad_hist_exp = histlite.hist(bad_exp_sigma,bins=100,log=True)
bad_hist_mc = histlite.hist(bad_mc_sigma,bins=100,log=True)

stefan_exp_sigma = stefan_exp['sigma']
stefan_mc_sigma = stefan_mc['sigma']

stefan_hist_exp = histlite.hist(stefan_exp_sigma,bins=100,log=True)
stefan_hist_mc = histlite.hist(stefan_mc_sigma,bins=100,log=True)

## Now to plot. ##
fig_sigma_exp = plt.figure()
ax=plt.gca()
histlite.plot1d(ax,ben_hist_exp,label = 'corrected', color = 'blue', linestyle = '-', errorbars=True) 
histlite.plot1d(ax,bad_hist_exp,label = 'original', color = 'red') 
histlite.plot1d(ax,stefan_hist_exp,label = 'stefan', color = 'k') 
ax.set_title(r'Sigma Distribution - Data')
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'density') 
#ax.set_xlim(-1.,1.)
#ax.set_ylim(1e-12,1e-10)
ax.loglog()
ax.grid() 
plt.legend(loc='upper right')
fig_sigma_exp.savefig('/home/brelethford/public_html/exp_sigmas.pdf')

fig_sigma_mc = plt.figure()
ax=plt.gca()
histlite.plot1d(ax,ben_hist_mc,label = 'corrected', color = 'blue',linestyle = '-', errorbars=True) 
histlite.plot1d(ax,bad_hist_mc,label = 'original', color = 'red') 
histlite.plot1d(ax,stefan_hist_mc,label = 'stefan', color = 'k') 
ax.set_title(r'Sigma Distribution - MC')
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$ensity') 
ax.set_xlim(1e-4,10)
#ax.set_ylim(1e-12,1e-10)
ax.loglog()
ax.grid() 
plt.legend(loc='upper right')
fig_sigma_mc.savefig('/home/brelethford/public_html/mc_sigmas.pdf')
