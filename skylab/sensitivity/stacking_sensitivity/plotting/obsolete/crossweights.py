import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import scipy as sp
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
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
    mu_stack = [stack_file['n_inj'] for stack_file in stack_files]
    results_stack = sorted(zip(gammas,mu_stack,flux_stack))
    gammas,mu,sensitivity = zip(*results_stack)
    sens100 = [refChange(sensitivity[i],gammas[i]) for i in range(len(sensitivity))]
    return gammas,mu,sens100

def plotDisc(datafolder):
    sensefolder = datafolder+'disc/'
    stack_files = [cache.load(sensefolder+file) for file in os.listdir(misc.ensure_dir(sensefolder)) if file.startswith('gamma')]    
    ##Let's find the gammas used first.
    gammas=[float(filename[5:-6]) for filename in  os.listdir(misc.ensure_dir(sensefolder)) if filename.startswith('gamma')]
    flux_stack = [stack_file['flux'] for stack_file in stack_files]
    mu_stack = [stack_file['n_inj'] for stack_file in stack_files]
    results_stack = sorted(zip(gammas,mu_stack,flux_stack))
    gammas,mu,disc = zip(*results_stack)
    disc100 = [refChange(disc[i],gammas[i]) for i in range(len(disc))]
    return gammas,mu,disc100

#I'll plot both the sens and disc now. I'll need to do this for each of two injection models.
gamma_sens_llheq_injflux, mu_sens_llheq_injflux, sensitivity_llheq_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/flux_inj/')
gamma_sens_llheq_injz, mu_sens_llheq_injz, sensitivity_llheq_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/redshift_inj/')

gamma_sens_llhflux_injflux, mu_sens_llhflux_injflux, sensitivity_llhflux_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/flux_inj/')
gamma_sens_llhflux_injz, mu_sens_llhflux_injz, sensitivity_llhflux_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/redshift_inj/')

gamma_sens_llhz_injflux, mu_sens_llhz_injflux, sensitivity_llhz_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/flux_inj/')
gamma_sens_llhz_injz, mu_sens_llhz_injz, sensitivity_llhz_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/redshift_inj/')


##Okay, now let's do all the same stuff, but for discovery potentials.
gamma_disc_llheq_injflux, mu_disc_llheq_injflux, disc_llheq_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/flux_inj/')
gamma_disc_llheq_injz, mu_disc_llheq_injz, disc_llheq_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/uniform/redshift_inj/')

gamma_disc_llhflux_injflux, mu_disc_llhflux_injflux, disc_llhflux_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/flux_inj/')
gamma_disc_llhflux_injz, mu_disc_llhflux_injz, disc_llhflux_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/redshift_inj/')

gamma_disc_llhz_injflux, mu_disc_llhz_injflux, disc_llhz_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/flux_inj/')
gamma_disc_llhz_injz, mu_disc_llhz_injz, disc_llhz_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/redshift/redshift_inj/')

fig_sensdisc = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injflux,sensitivity_llheq_injflux, '--', label='llhequal, injflux', color = 'black')
plt.plot(gamma_sens_llheq_injz,sensitivity_llheq_injz, '--', label='llhequal, injz', color = 'dimgrey')

plt.plot(gamma_sens_llhflux_injflux,sensitivity_llhflux_injflux, '--', label='llhflux, injflux', color = 'green')
plt.plot(gamma_sens_llhflux_injz,sensitivity_llhflux_injz, '--', label='llhflux, injz', color = 'lightgreen')

plt.plot(gamma_sens_llhz_injflux,sensitivity_llhz_injflux, '--', label='llhz, injflux', color = 'red')
plt.plot(gamma_sens_llhz_injz,sensitivity_llhz_injz, '--', label='llhz, injz', color = 'orangered')

plt.plot(gamma_disc_llheq_injflux,disc_llheq_injflux, color = 'black')
plt.plot(gamma_disc_llheq_injz,disc_llheq_injz, color = 'dimgrey')

plt.plot(gamma_disc_llhflux_injflux,disc_llhflux_injflux, color = 'green')
plt.plot(gamma_disc_llhflux_injz,disc_llhflux_injz, color = 'lightgreen')

plt.plot(gamma_disc_llhz_injflux,disc_llhz_injflux, color = 'red')
plt.plot(gamma_disc_llhz_injz,disc_llhz_injz, color = 'orangered')

ax.set_title(r'Stacked sensitivities/discovery - no cut')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$E^2 \cdot \left( \frac{E}{100 \textrm{\small{TeV}}} \right)  ^{\gamma-2} \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$')
ax.set_xlim(1.,4.)
ax.semilogy()
ax.grid() 
plt.legend(loc='lower center',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_sensdisc.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc.pdf')
fig_sensdisc.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc.png')

########################################
#Okay, Now let's plot injected events:

fig_mu = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injflux,mu_sens_llheq_injflux, '--', label='llhequal, injflux', color = 'black')
plt.plot(gamma_sens_llheq_injz,mu_sens_llheq_injz, '--', label='llhequal, injz', color = 'dimgrey')

plt.plot(gamma_sens_llhflux_injflux,mu_sens_llhflux_injflux, '--', label='llhflux, injflux', color = 'green')
plt.plot(gamma_sens_llhflux_injz,mu_sens_llhflux_injz, '--', label='llhflux, injz', color = 'lightgreen')

plt.plot(gamma_sens_llhz_injflux,mu_sens_llhz_injflux, '--', label='llhz, injflux', color = 'red')
plt.plot(gamma_sens_llhz_injz,mu_sens_llhz_injz, '--', label='llhz, injz', color = 'orangered')

plt.plot(gamma_disc_llheq_injflux,mu_disc_llheq_injflux, color = 'black')
plt.plot(gamma_disc_llheq_injz,mu_disc_llheq_injz, color = 'dimgrey')

plt.plot(gamma_disc_llhflux_injflux,mu_disc_llhflux_injflux, color = 'green')
plt.plot(gamma_disc_llhflux_injz,mu_disc_llhflux_injz, color = 'lightgreen')

plt.plot(gamma_disc_llhz_injflux,mu_disc_llhz_injflux, color = 'red')
plt.plot(gamma_disc_llhz_injz,mu_disc_llhz_injz, color = 'orangered')

ax.set_title(r'Stacked injected events - no cut')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$n_{inj}$')
ax.grid() 
plt.legend(loc='lower center',prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2) 
fig_mu.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_mu.pdf')
fig_mu.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_mu.png')

def printSens(gammalist,sensitivitylist,mulist):
  sens=0
  for i in range(len(gammalist)):
    if gammalist[i] == 2.0:
      sens = sensitivitylist[i]
      mu = mulist[i]
  return 'sens: '+str(sens)+', mu: '+  str(mu)
#Let's show how the sensitivities look for no energy restriction.

print(' ')
print ('sensitivities: no restriction on energy')

print ('llheq, injflux: =' + printSens(gamma_sens_llheq_injflux,sensitivity_llheq_injflux,mu_sens_llheq_injflux))
print ('llheq, injz: =' + printSens(gamma_sens_llheq_injz,sensitivity_llheq_injz,mu_sens_llheq_injz))
print ('llhflux, injflux: =' + printSens(gamma_sens_llhflux_injflux,sensitivity_llhflux_injflux,mu_sens_llhflux_injflux))
print ('llhflux, injz: =' + printSens(gamma_sens_llhflux_injz,sensitivity_llhflux_injz,mu_sens_llhflux_injz))
print ('llhz, injflux: =' + printSens(gamma_sens_llhz_injflux,sensitivity_llhz_injflux,mu_sens_llhz_injflux))
print ('llhz, injz: =' + printSens(gamma_sens_llhz_injz,sensitivity_llhz_injz,mu_sens_llhz_injz))

#Okay, now we'll do the same thing, but for a dataset with the injection range restricted



##########################################
pad = 0.26
def icprelim (fig, x=pad + .12, y=1 - pad + .1, **kw):
   """Mark a figure as preliminary."""
   if 'color' not in kw:
       kw['color'] = 'red'
   if 'weight' not in kw:
       kw['weight'] = 'bold'
   fig.text (x, y, 'IceCube Preliminary', **kw)

gamma_sens_llheq_injflux, mu_sens_llheq_injflux, sensitivity_llheq_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/uniform/flux_inj/')
gamma_sens_llheq_injz, mu_sens_llheq_injz, sensitivity_llheq_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/uniform/redshift_inj/')

gamma_sens_llhflux_injflux, mu_sens_llhflux_injflux, sensitivity_llhflux_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/flux/flux_inj/')
gamma_sens_llhflux_injz, mu_sens_llhflux_injz, sensitivity_llhflux_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/flux/redshift_inj/')

gamma_sens_llhz_injflux, mu_sens_llhz_injflux, sensitivity_llhz_injflux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/redshift/flux_inj/')
gamma_sens_llhz_injz, mu_sens_llhz_injz, sensitivity_llhz_injz = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/redshift/redshift_inj/')


##Okay, now let's do all the same stuff, but for discovery potentials.
gamma_disc_llheq_injflux, mu_disc_llheq_injflux, disc_llheq_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/uniform/flux_inj/')
gamma_disc_llheq_injz, mu_disc_llheq_injz, disc_llheq_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/uniform/redshift_inj/')

gamma_disc_llhflux_injflux, mu_disc_llhflux_injflux, disc_llhflux_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/flux/flux_inj/')
gamma_disc_llhflux_injz, mu_disc_llhflux_injz, disc_llhflux_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/flux/redshift_inj/')

gamma_disc_llhz_injflux, mu_disc_llhz_injflux, disc_llhz_injflux = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/redshift/flux_inj/')
gamma_disc_llhz_injz, mu_disc_llhz_injz, disc_llhz_injz = plotDisc( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/redshift/redshift_inj/')


#flux injected plots
fig_sensdisc_cut_flux = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injflux,sensitivity_llheq_injflux, '--', label='equal', color = 'black')

plt.plot(gamma_sens_llhflux_injflux,sensitivity_llhflux_injflux, '--', label='flux', color = 'green')

plt.plot(gamma_sens_llhz_injflux,sensitivity_llhz_injflux, '--', label='redshift', color = 'red')

plt.plot(gamma_disc_llheq_injflux,disc_llheq_injflux, color = 'black')

plt.plot(gamma_disc_llhflux_injflux,disc_llhflux_injflux, color = 'green')

plt.plot(gamma_disc_llhz_injflux,disc_llhz_injflux, color = 'red')

ax.set_title(r'flux injection')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$E^2 \cdot \left( \frac{E}{100 \textrm{\small{TeV}}} \right)  ^{\gamma-2} \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$')
ax.set_xlim(1.,4.)
ax.semilogy()
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall, title=r'\small{TS weighting model}')
plt.subplots_adjust (left=.2, bottom=.2) 
icprelim(fig_sensdisc_cut_flux,pad)

fig_sensdisc_cut_flux.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc_flux.pdf')
fig_sensdisc_cut_flux.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc_flux.png')
#Now for redshift injected plots
fig_sensdisc_cut_redshift = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injz,sensitivity_llheq_injz, '--', label='equal', color = 'black')

plt.plot(gamma_sens_llhflux_injz,sensitivity_llhflux_injz, '--', label='flux', color = 'green')

plt.plot(gamma_sens_llhz_injz,sensitivity_llhz_injz, '--', label='redshift', color = 'red')

plt.plot(gamma_disc_llheq_injz,disc_llheq_injz, color = 'black')

plt.plot(gamma_disc_llhflux_injz,disc_llhflux_injz, color = 'green')

plt.plot(gamma_disc_llhz_injz,disc_llhz_injz, color = 'red')


ax.set_title(r'redshift injection')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$E^2 \cdot \left( \frac{E}{100 \textrm{\small{TeV}}} \right)  ^{\gamma-2} \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$')
ax.set_xlim(1.,4.)
ax.semilogy()
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall, title=r'\small{TS weighting model}')
plt.subplots_adjust (left=.2, bottom=.2) 

icprelim(fig_sensdisc_cut_redshift,pad)

fig_sensdisc_cut_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc_redshift.pdf')
fig_sensdisc_cut_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/stacked_sensdisc_redshift.png')


#And finally, injected event plots
fig_mu_cut_flux = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injflux,mu_sens_llheq_injflux, '--', label='equal', color = 'black')

plt.plot(gamma_sens_llhflux_injflux,mu_sens_llhflux_injflux, '--', label='flux', color = 'green')

plt.plot(gamma_sens_llhz_injflux,mu_sens_llhz_injflux, '--', label='redshift', color = 'red')

plt.plot(gamma_disc_llheq_injflux,mu_disc_llheq_injflux, color = 'black')

plt.plot(gamma_disc_llhflux_injflux,mu_disc_llhflux_injflux, color = 'green')

plt.plot(gamma_disc_llhz_injflux,mu_disc_llhz_injflux, color = 'red')

ax.set_title(r'flux $n_{inj}$')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$n_{inj}$')
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall, title=r'\small{TS weighting model}')
plt.subplots_adjust (left=.2, bottom=.2) 
icprelim(fig_mu_cut_flux,pad)

fig_mu_cut_flux.savefig('/data/user/brelethford/AGN_Core/Plots/mu_cut_flux.pdf')
fig_sensdisc_cut_flux.savefig('/data/user/brelethford/AGN_Core/Plots/mu_cut_flux.png')
#Now for redshift injected plots
fig_mu_cut_redshift = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(gamma_sens_llheq_injz,mu_sens_llheq_injz, '--', label='equal', color = 'black')

plt.plot(gamma_sens_llhflux_injz,mu_sens_llhflux_injz, '--', label='flux', color = 'green')

plt.plot(gamma_sens_llhz_injz,mu_sens_llhz_injz, '--', label='redshift', color = 'red')

plt.plot(gamma_disc_llheq_injz,mu_disc_llheq_injz, color = 'black')

plt.plot(gamma_disc_llhflux_injz,mu_disc_llhflux_injz, color = 'green')

plt.plot(gamma_disc_llhz_injz,mu_disc_llhz_injz, color = 'red')


ax.set_title(r'redshift $n_{inj}$')
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$n_{inj}$')
ax.grid() 
plt.legend(loc='lower right',prop=propxxsmall, title=r'\small{TS weighting model}')
plt.subplots_adjust (left=.2, bottom=.2) 

icprelim(fig_mu_cut_redshift,pad)

fig_mu_cut_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/mu_cut_redshift.pdf')
fig_mu_cut_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/mu_cut_redshift.png')









print (' ')
print ('sensitivities: injection energy restriction')

print ('llheq, injflux: =' + printSens(gamma_sens_llheq_injflux,sensitivity_llheq_injflux,mu_sens_llheq_injflux))
print ('llheq, injz: =' + printSens(gamma_sens_llheq_injz,sensitivity_llheq_injz,mu_sens_llheq_injz))
print ('llhflux, injflux: =' + printSens(gamma_sens_llhflux_injflux,sensitivity_llhflux_injflux,mu_sens_llhflux_injflux))
print ('llhflux, injz: =' + printSens(gamma_sens_llhflux_injz,sensitivity_llhflux_injz,mu_sens_llhflux_injz))
print ('llhz, injflux: =' + printSens(gamma_sens_llhz_injflux,sensitivity_llhz_injflux,mu_sens_llhz_injflux))
print ('llhz, injz: =' + printSens(gamma_sens_llhz_injz,sensitivity_llhz_injz,mu_sens_llhz_injz))


