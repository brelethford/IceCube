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

projfolder='/data/user/brelethford/'
datafolder=projfolder+'Data/AGN_Core_Sample/'
filename_plots=projfolder + 'AGN_Core/Plots/tests/sigma_error/'

filename_pickle = datafolder+'pickle/'
if not os.path.exists(filename_pickle):
	os.mkdir(filename_pickle)

import data_multi

    ### Initial Plots ###
def plot_error(year):
    #The following currently devoted to making energy and angular error plots for a year of data.
    # init likelihood class
    if year == 40:
      llh = data_multi.init40(energy=True)
      mc = cache.load(filename_pickle+"IC40/mc.pickle")
      extra = cache.load(filename_pickle+"IC40/dpsi.pickle")
    if year == 59:
      llh = data_multi.init59(energy=True)
      mc = cache.load(filename_pickle+"IC59/mc.pickle")
      extra = cache.load(filename_pickle+"IC59/dpsi.pickle")
    if year == 79:
      llh = data_multi.init79(energy=True)
      mc = cache.load(filename_pickle+"IC79/mc.pickle")
      extra = cache.load(filename_pickle+"IC79/dpsi.pickle")
    elif year == 86:
      llh = data_multi.init86I(energy=True)
      mc = cache.load(filename_pickle+"IC86I/mc.pickle")
      extra = cache.load(filename_pickle+"IC86I/dpsi.pickle")
    dpsi = extra['dpsi']

    print(llh)

    # datatest
    #Currently don't need to remake these plots. but DONT DELETE

    colors=['b','g','y','r']
    gamma = np.linspace(1., 2.7, 4)

   # fig_energy, (ax1, ax2) = plt.subplots(ncols=2)
   # ax1.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Pseudo-Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # ax1.legend(loc="best")
   # ax1.set_title("Reconstructed Energy - IC{}".format(str(year)))
   # ax1.set_xlabel("logE")
   # ax1.set_ylabel("Relative Abundance")
   # ax1.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # ax2.set_title("Zoomed In")
   # ax2.set_xlabel("logE")
   # ax2.set_xlim(4,10)
   # ax2.set_ylim(1e-5,1)
   # ax2.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Pseudo-Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # fig_energy.savefig(filename_plots + 'energy_hists_IC{}.pdf'.format(str(year)))

    fig_angular_error, (ax1,ax2) = plt.subplots(ncols=2,figsize = (10,5))
    dec = np.arcsin(mc["sinDec"])
    angdist = np.degrees(dpsi) 

    ax1.hist([np.log10(np.degrees(mc["sigma"])) for i in gamma], label = [r"$\sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dashed',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    ax1.hist([np.log10(angdist) for i in gamma], label = [r"$\Delta \psi$ - $\gamma={0:.1f}$".format(g) for g in gamma],
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], linestyle = 'solid', color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)
    ax1.set_title("Reco MC Angular Error Check - IC{}".format(str(year)))
    ax1.set_xlabel(r"log$\sigma_{ang}$ (degrees)")
    ax1.set_ylabel("Relative Abundance")
    ax1.set_ylim(0,1.5)


    ax2.hist([(np.degrees(mc["sigma"])) for i in gamma], label = [r"$\sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dashed',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=1000, normed=True)

    ax2.hist([(angdist) for i in gamma], label = [r"$\Delta \psi$ - $\gamma={0:.1f}$".format(g) for g in gamma],
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], linestyle = 'solid', color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=1000, normed=True)
    ax2.legend(loc="upper right")
    ax2.set_xlim(0,5)
    ax2.set_ylim(0,3.5)
    ax2.set_xlabel(r"$\sigma_{ang}$ (degrees)")
    fig_angular_error.savefig(filename_plots + 'angular_error_hists_IC{}.pdf'.format(str(year)))

plot_error(40)
plot_error(59)
plot_error(79)
plot_error(86)


#Now let's also have plots for dpsi/sigma.
def plot_ratio(year):
    #The following currently devoted to making energy and angular error plots for a year of data.
    # init likelihood class
    if year == 40:
      llh = data_multi.init40(energy=True)
      mc = cache.load(filename_pickle+"IC40/mc.pickle")
      extra = cache.load(filename_pickle+"IC40/dpsi.pickle")
    if year == 59:
      llh = data_multi.init59(energy=True)
      mc = cache.load(filename_pickle+"IC59/mc.pickle")
      extra = cache.load(filename_pickle+"IC59/dpsi.pickle")
    if year == 79:
      llh = data_multi.init79(energy=True)
      mc = cache.load(filename_pickle+"IC79/mc.pickle")
      extra = cache.load(filename_pickle+"IC79/dpsi.pickle")
    elif year == 86:
      llh = data_multi.init86I(energy=True)
      mc = cache.load(filename_pickle+"IC86I/mc.pickle")
      extra = cache.load(filename_pickle+"IC86I/dpsi.pickle")
    dpsi = extra['dpsi']

    print(llh)

    # datatest
    #Currently don't need to remake these plots. but DONT DELETE

    colors=['b','g','y','r']
    gamma = np.linspace(1., 2.7, 4)

   # fig_energy, (ax1, ax2) = plt.subplots(ncols=2)
   # ax1.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Pseudo-Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # ax1.legend(loc="best")
   # ax1.set_title("Reconstructed Energy - IC{}".format(str(year)))
   # ax1.set_xlabel("logE")
   # ax1.set_ylabel("Relative Abundance")
   # ax1.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # ax2.set_title("Zoomed In")
   # ax2.set_xlabel("logE")
   # ax2.set_xlim(4,10)
   # ax2.set_ylim(1e-5,1)
   # ax2.hist([llh.exp["logE"]] + [mc["logE"] for i in gamma],
   #          weights=[np.ones(len(llh.exp))]
   #                   + [mc["ow"] * mc["trueE"]**(-g) for g in gamma],
   #          label=["Pseudo-Data"] + [r"$\gamma={0:.1f}$".format(g) for g in gamma], color= ['k'] + [colors[g] for g in range(len(gamma))],
   #          histtype="step", bins=100, log=True, normed=True, cumulative=-1)
   # fig_energy.savefig(filename_plots + 'energy_hists_IC{}.pdf'.format(str(year)))

    fig_ratio, (ax1,ax2) = plt.subplots(ncols=2,figsize=(10,5))
    dec = np.arcsin(mc["sinDec"])
    angdist = np.degrees(dpsi) 

    ax1.hist([np.log10(np.degrees(mc["sigma"])/(angdist)) for i in gamma], label = [r"$\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'solid',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    medians = [misc.weighted_median(np.log10(np.degrees(mc["sigma"])/(angdist)),mc["ow"] * mc["trueE"]**(-g)) for g in gamma]

    for g in range(len(gamma)):
      ax1.axvline(medians[g], color=colors[g], alpha = 0.3)

    ax1.set_title(r"Reco MC $\Delta \psi / \sigma$ Check - IC{}".format(str(year)))
    ax1.set_xlabel(r"log$\Delta \psi / \sigma_{ang}$")
    ax1.set_ylabel("Relative Abundance")
    ax1.set_ylim(0,1.5)
    ax1.set_xlim(-2.5,2.5)
    ax1.axvline(x=0, color='k')
    ax1.axvline(x=np.log10(1.1774), color='k')
    ax2.set_xlim(-5,5)


    ax2.hist([(np.degrees(mc["sigma"]))/(angdist) for i in gamma], label = [r"$\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'solid',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=1000, range = (0,5), normed=True)

    ax2.legend(loc="upper right")
    ax2.set_xlim(0,5)
    ax2.set_ylim(0,3.5)
    ax2.set_xlabel(r"$\Delta \psi / \sigma_{ang}$")
    fig_ratio.savefig(filename_plots + 'dpsi_sigma_ratio_IC{}.pdf'.format(str(year)))

plot_ratio(40)
plot_ratio(59)
plot_ratio(79)
plot_ratio(86)
