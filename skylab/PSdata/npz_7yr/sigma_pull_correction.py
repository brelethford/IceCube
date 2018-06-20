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
filename_pickle='/data/user/coenders/data/MultiYearPointSource/npz/'

filename_plots=projfolder + 'AGN_Core/Plots/tests/npz_sigma_error/'



    ### Initial Plots ###
def plot_error(year):
    #The following currently devoted to making energy and angular error plots for a year of data.
    # init mc
    mc = np.load(filename_pickle+"IC{}_MC.npy".format(year))
    dpsi=[astro.angular_distance(mc['trueRa'][i],mc['trueDec'][i],mc['ra'][i],np.arcsin(mc['sinDec'][i])) for i in range(len(mc))]

    mc_corr = np.load(filename_pickle+"IC{}_corrected_MC.npy".format(year))
    dpsi_corr=[astro.angular_distance(mc_corr['trueRa'][i],mc_corr['trueDec'][i],mc_corr['ra'][i],np.arcsin(mc_corr['sinDec'][i])) for i in range(len(mc_corr))]

    colors=['g']#['b','g','y','r']
    colors_corr=['r']#['b','g','y','r']
    gamma = np.array([2.])#np.linspace(1., 2.7, 4)

    fig_angular_error, (ax1,ax2) = plt.subplots(ncols=2,figsize = (10,5))
    dec = np.arcsin(mc["sinDec"])
    dec_corr = np.arcsin(mc_corr["sinDec"])
    angdist = np.degrees(dpsi) 
    angdist_corr = np.degrees(dpsi_corr) 

    ax1.hist([np.log10(np.degrees(mc["sigma"])) for i in gamma], label = [r"$\sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dashed',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    ax1.hist([np.log10(angdist) for i in gamma], label = [r"$\Delta \psi$ - $\gamma={0:.1f}$".format(g) for g in gamma],
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], linestyle = 'solid', color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    ax1.hist([np.log10(np.degrees(mc_corr["sigma"])) for i in gamma], label = [r"$\sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dashdot',
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    ax1.hist([np.log10(angdist_corr) for i in gamma], label = [r"'Corrected' $\Delta \psi$ - $\gamma={0:.1f}$".format(g) for g in gamma],
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], linestyle = 'dotted', color=[colors_corr[g] for g in range(len(gamma))],
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

    ax2.hist([(np.degrees(mc_corr["sigma"])) for i in gamma], label = [r"'corrected' $\sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dashdot',
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
             histtype="step", bins=1000, normed=True)

    ax2.hist([(angdist_corr) for i in gamma], label = [r"$\Delta \psi$ - $\gamma={0:.1f}$".format(g) for g in gamma],
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], linestyle = 'dotted', color=[colors_corr[g] for g in range(len(gamma))],
             histtype="step", bins=1000, normed=True)
    ax2.legend(loc="upper right")
    ax2.set_xlim(0,5)
    ax2.set_ylim(0,3.5)
    ax2.set_xlabel(r"$\sigma_{ang}$ (degrees)")
    fig_angular_error.savefig(filename_plots + 'angular_error_hists_IC{}.pdf'.format(str(year)))

#plot_error(40)
#plot_error(59)
#plot_error(79)
#plot_error(86)


#Now let's also have plots for dpsi/sigma.
def plot_ratio(year):
    #The following currently devoted to making energy and angular error plots for a year of data.
    # init likelihood class
    mc = np.load(filename_pickle+"IC{}_MC.npy".format(year))
    dpsi=[astro.angular_distance(mc['trueRa'][i],mc['trueDec'][i],mc['ra'][i],np.arcsin(mc['sinDec'][i])) for i in range(len(mc))]

    mc_corr = np.load(filename_pickle+"IC{}_corrected_MC.npy".format(year))
    dpsi_corr=[astro.angular_distance(mc_corr['trueRa'][i],mc_corr['trueDec'][i],mc_corr['ra'][i],np.arcsin(mc_corr['sinDec'][i])) for i in range(len(mc_corr))]

    colors=['g']#['b','g','y','r']
    colors_corr=['r']#['b','g','y','r']
    gamma = np.array([2.])#np.linspace(1., 2.7, 4)

    fig_ratio, (ax1,ax2) = plt.subplots(ncols=2,figsize=(10,5))
    dec = np.arcsin(mc["sinDec"])
    angdist = np.degrees(dpsi) 
    dec_corr = np.arcsin(mc_corr["sinDec"])
    angdist_corr = np.degrees(dpsi_corr) 

    ax1.hist([np.log10(np.degrees(mc["sigma"])/(angdist)) for i in gamma], label = [r"$\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'solid',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    medians = [misc.weighted_median(np.log10(np.degrees(mc["sigma"])/(angdist)),mc["ow"] * mc["trueE"]**(-g)) for g in gamma]

    for g in range(len(gamma)):
      ax1.axvline(medians[g], color=colors[g], linestyle = 'dotted', alpha = 0.3)

    ax1.hist([np.log10(np.degrees(mc_corr["sigma"])/(angdist_corr)) for i in gamma], label = [r"'Corrected' $\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dotted',
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
             histtype="step", bins=100, normed=True)

    medians = [misc.weighted_median(np.log10(np.degrees(mc_corr["sigma"])/(angdist_corr)),mc["ow"] * mc["trueE"]**(-g)) for g in gamma]

    for g in range(len(gamma)):
      ax1.axvline(medians[g], color=colors_corr[g], marker='+',alpha = 0.3)

    ax1.set_title(r"Reco MC $\Delta \psi / \sigma$ Check - IC{}".format(str(year)))
    ax1.set_xlabel(r"log$\Delta \psi / \sigma_{ang}$")
    ax1.set_ylabel("Relative Abundance")
    ax1.set_ylim(0,1.5)
    ax1.set_xlim(-2.5,2.5)
    #ax1.axvline(x=0, color='k')
    ax1.axvline(x=np.log10(1./1.1774), color='k')
    ax2.set_xlim(-5,5)


    ax2.hist([(np.degrees(mc["sigma"]))/(angdist) for i in gamma], label = [r"$\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'solid',
             weights=[mc["ow"] * mc["trueE"]**(-g) for g in gamma], color=[colors[g] for g in range(len(gamma))],
             histtype="step", bins=1000, range = (0,5), normed=True)

    ax2.hist([(np.degrees(mc_corr["sigma"]))/(angdist_corr) for i in gamma], label = [r"'Corrected' $\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dotted',
             weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
             histtype="step", bins=1000, range = (0,5), normed=True)

    ax2.legend(loc="upper right")
    ax2.set_xlim(0,5)
    ax2.set_ylim(0,3.5)
    ax2.set_xlabel(r"$\Delta \psi / \sigma_{ang}$")
    fig_ratio.savefig(filename_plots + 'dpsi_sigma_ratio_IC{}.pdf'.format(str(year)))

#plot_ratio(40)
#plot_ratio(59)
#plot_ratio(79)
#plot_ratio(86)

#Finally let's also have 2d pull plots.
def plotpull2d(year):
    #The following currently devoted to making energy and angular error plots for a year of data.
    # init likelihood class
    mc = np.load(filename_pickle+"IC{}_MC.npy".format(year))
    dpsi=[astro.angular_distance(mc['trueRa'][i],mc['trueDec'][i],mc['ra'][i],np.arcsin(mc['sinDec'][i])) for i in range(len(mc))]

    mc_corr = np.load(filename_pickle+"IC{}_corrected_MC.npy".format(year))
    dpsi_corr=[astro.angular_distance(mc_corr['trueRa'][i],mc_corr['trueDec'][i],mc_corr['ra'][i],np.arcsin(mc_corr['sinDec'][i])) for i in range(len(mc_corr))]

    colors=['g']#['b','g','y','r']
    colors_corr=['r']#['b','g','y','r']
    gamma = np.array([2.])#np.linspace(1., 2.7, 4)

    fig_pull, (ax) = plt.subplots(ncols=1,figsize=(10,5))
    dec = np.arcsin(mc["sinDec"])
    angdist = np.degrees(dpsi) 
    dec_corr = np.arcsin(mc_corr["sinDec"])
    angdist_corr = np.degrees(dpsi_corr) 

    #medians = [misc.weighted_median(np.log10(np.degrees(mc["sigma"])/(angdist)),mc["ow"] * mc["trueE"]**(-g)) for g in gamma]

    #for g in range(len(gamma)):
    #  ax1.axvline(medians[g], color=colors[g], linestyle = 'dotted', alpha = 0.3)

    #ax1.hist([np.log10(np.degrees(mc_corr["sigma"])/(angdist_corr)) for i in gamma], label = [r"'Corrected' $\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dotted',
    #         weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
    #         histtype="step", bins=100, normed=True)

    #medians = [misc.weighted_median(np.log10(np.degrees(mc_corr["sigma"])/(angdist_corr)),mc["ow"] * mc["trueE"]**(-g)) for g in gamma]

    #for g in range(len(gamma)):
    #  ax1.axvline(medians[g], color=colors_corr[g], marker='+',alpha = 0.3)

    #ax1.set_title(r"Reco MC $\Delta \psi / \sigma$ Check - IC{}".format(str(year)))
    #ax1.set_xlabel(r"log$\Delta \psi / \sigma_{ang}$")
    #ax1.set_ylabel("Relative Abundance")
    #ax1.set_ylim(0,1.5)
    #ax1.set_xlim(-2.5,2.5)
    #ax1.axvline(x=0, color='k')
    #ax1.axvline(x=np.log10(1./1.1774), color='k')
    #set gamma
    g=2.0
    pull = np.log10(np.degrees(mc["sigma"])/(angdist))
    #Need to calc the median of pull for each energyrange:
    pullrange = (-4,4)
    erange = (3,8)
    ebins = 50
    pullbins = 200
    eedges = np.linspace(erange[0],erange[1],ebins+1)
    ezones=zip(eedges[:-1],eedges[1:])
    emid = [np.median(ezone) for ezone in ezones]
    emasks = [(np.log10(mc["trueE"]) > ezone[0]) & (np.log10(mc["trueE"]) < ezone[1]) for ezone in ezones]

    medians = [misc.weighted_median(pull[emask],mc["ow"][emask] * mc["trueE"][emask]**(-g)) for emask in emasks]
    #pull_corr = [(np.degrees(mc_corr["sigma"]))/(angdist_corr) for i in gamma]
    hpull = histlite.hist((np.log10(mc['trueE']),pull),
                           weights = mc["ow"] * mc["trueE"]**(-g),
                           bins=(ebins,pullbins), range = (erange,pullrange), log = (False,False))

    histlite.plot2d(ax, hpull.normalize([-1]), label = [r"$\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g)], cbar=True, cmap='jet', color=colors[0])
    
    ax.scatter(emid,medians, color = 'white')             
    ax.axhline(y=np.log10(1.1774), color='white', linestyle = 'dashed')

    #ax.hist([(np.degrees(mc_corr["sigma"]))/(angdist_corr) for i in gamma], label = [r"'Corrected' $\Delta \psi / \sigma$ - $\gamma={0:.1f}$".format(g) for g in gamma], linestyle = 'dotted',
    #         weights=[mc_corr["ow"] * mc_corr["trueE"]**(-g) for g in gamma], color=[colors_corr[g] for g in range(len(gamma))],
    #         histtype="step", bins=1000, range = (0,5), normed=True)
    ax.set_title("Pull for IC{}".format(str(year)))
    ax.legend(loc="upper right")
    ax.set_xlim(3,8)
    ax.set_xlabel("log(trueE[GeV])")
    ax.set_ylim(-2,2)
    ax.set_ylabel(r"log($\Delta \psi / \sigma_{ang}$)")
    fig_pull.savefig(filename_plots + 'opp_pull_IC{}.pdf'.format(str(year)))

#plotpull2d(40)
#plotpull2d(59)
#plotpull2d(79)
plotpull2d(86)
