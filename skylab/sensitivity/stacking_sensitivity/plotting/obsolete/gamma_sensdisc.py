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
###This script imports the sensitivities from the submitter and plots them.###
#datafolder = '/data/user/brelethford/Output/all_sky_sensitivity/results/box/' ##Remember to change this depending on which results you're looking at!##
###########################################
#This function should pick out the various sensitivities for different gammas.
###########################################

misc.tex_mpl_rc()
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')
propxxsmall = mpl.font_manager.FontProperties (size='xx-small')
w=4

#Okay, so we need to make sure that our fluxes are at 100 TeV reference energy. Currently they're at 1 TeV, so we'll have to convert them based on gamma.

def refChange(flux,gamma):
    return np.multiply(flux,100**(2-gamma))

def plotSensitivity(datafolder):
    sensefolder = datafolder+'sensitivity/'
    stack_files = [cache.load(sensefolder+file) for file in os.listdir(misc.ensure_dir(sensefolder)) if file.startswith('gamma')]    
    ##Let's find the gammas used first.
    gammas=[float(filename[5:-6]) for filename in  os.listdir(misc.ensure_dir(sensefolder)) if filename.startswith('gamma')]
    flux_stack = [stack_file['flux'] for stack_file in stack_files]
    inj_stack = [stack_file['mu'] for stack_file in stack_files]
    results_stack = sorted(zip(gammas,inj_stack,flux_stack))
    gammas,inj,sensitivity = zip(*results_stack)
    sens100 = [refChange(sensitivity[i],gammas[i]) for i in range(len(sensitivity))]
    return gammas,inj,sens100

def plotDisc(datafolder):
    sensefolder = datafolder+'disc/'
    stack_files = [cache.load(sensefolder+file) for file in os.listdir(misc.ensure_dir(sensefolder)) if file.startswith('gamma')]    
    ##Let's find the gammas used first.
    gammas=[float(filename[5:-6]) for filename in  os.listdir(misc.ensure_dir(sensefolder)) if filename.startswith('gamma')]
    flux_stack = [stack_file['flux'] for stack_file in stack_files]
    inj_stack = [stack_file['mu'] for stack_file in stack_files]
    results_stack = sorted(zip(gammas,inj_stack,flux_stack))
    gammas,inj,disc = zip(*results_stack)
    disc100 = [refChange(disc[i],gammas[i]) for i in range(len(disc))]
    return gammas,inj,disc100

#I'll plot both the sens and disc now.
gamma_uniform, inj_uniform, sensitivity_uniform = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/')
gamma_flux, inj_flux, sensitivity_flux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/')
gamma_redshift, inj_redshift, sensitivity_redshift = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/')
'''
## Sensitivity plot. ##
fig_sensitivity = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_uniform,sensitivity_uniform, label='Uniform', color = 'black')
plt.plot(gamma_flux,sensitivity_flux, label='Flux', color = 'green')
plt.plot(gamma_redshift,sensitivity_redshift, label='Redshift', color = 'red')
ax.set_title(r'Stacked Sensitivity')
ax.set_xlabel(r'$\gamma$')
#ax.set_ylabel(r'$E^2 \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
ax.set_ylabel(r'$E^{\gamma} \frac{d\phi}{dE} [\frac{TeV^{\gamma-1}}{cm^2 s}]$') 
ax.set_xlim(1.,4.)
#ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid()
plt.legend(loc='lower right',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_sensitivity.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensitivity.pdf')

#Injected Events Plot

fig_injections = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_uniform,inj_uniform,label='Uniform', color = 'black')
plt.plot(gamma_flux,inj_flux, label='Flux', color = 'green')
plt.plot(gamma_redshift,inj_redshift, label='Redshift', color = 'red')
ax.set_title(r'Stacked Injections')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$n_{inj}$') 
ax.set_xlim(1.,4.)
#ax.set_ylim(1e-12,1e-10)
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_injections.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_injections.pdf')
#######################################
'''

##Okay, now let's do all the same stuff, but for discovery potentials.
gamma_uniform_disc, inj_uniform_disc, disc_uniform = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/')
gamma_flux_disc, inj_flux_disc, disc_flux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/')
gamma_redshift_disc, inj_redshift_disc, disc_redshift = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/')
'''
fig_disc = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_uniform_disc,disc_uniform, label='Uniform', color = 'black')
plt.plot(gamma_flux_disc,disc_flux, label='Flux', color = 'green')
plt.plot(gamma_redshift_disc,disc_redshift, label='Redshift', color = 'red')
ax.set_title(r'Stacked Discovery Potential')
ax.set_xlabel(r'$\gamma$')
#ax.set_ylabel(r'$E^2 \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
ax.set_ylabel(r'$E^{\gamma} \frac{d\phi}{dE} [\frac{TeV^{\gamma-1}}{cm^2 s}]$') 
ax.set_xlim(1.,4.)
#ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_disc.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_discovery.pdf')
'''

fig_sensdisc = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_uniform,sensitivity_uniform, '--', label='Sens - Equal', color = 'black')
plt.plot(gamma_uniform_disc,disc_uniform, label='Disc - Equal', color = 'black')
plt.plot(gamma_flux,sensitivity_flux, '--', label='Sens - Flux', color = 'green')
plt.plot(gamma_flux_disc,disc_flux, label='Disc - Flux', color = 'green')
plt.plot(gamma_redshift,sensitivity_redshift, '--', label='Sens - Redshift', color = 'red')
plt.plot(gamma_redshift_disc,disc_redshift, label='Disc - Redshift', color = 'red')
ax.set_title(r'Stacked sensitivities/discovery')
ax.set_xlabel(r'$\gamma$')
#ax.set_ylabel(r'$E^2 \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
#ax.set_ylabel(r'$E^{\gamma} \frac{d\phi}{dE} [\frac{TeV^{\gamma-1}}{cm^2 s}]$') 
ax.set_ylabel(r'$E^2 \cdot \left( \frac{E}{100 \textrm{\small{TeV}}} \right)  ^{\gamma-2} \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$')
ax.set_xlim(1.,4.)
#ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid() 
plt.legend(loc='lower center',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_sensdisc.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc.pdf')
fig_sensdisc.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc.png')
