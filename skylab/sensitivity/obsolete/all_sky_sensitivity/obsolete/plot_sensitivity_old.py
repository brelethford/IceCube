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
#################################

misc.tex_mpl_rc()
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')
propxxsmall = mpl.font_manager.FontProperties (size='xx-small')
w=4

stefanfolder = '/data/user/brelethford/Output/all_sky_sensitivity/stefan_c/IC86-I/'
stefanfiles = [cache.load(stefanfolder+file) for file in os.listdir(stefanfolder) if file.endswith('.pkl')]

stefan_flux = [file[2.0]['flux'][0] for file in stefanfiles]
stefan_inj = [file[2.0]['mu'][0] for file in stefanfiles]
stefan_sindec = [np.sin(file['dec']) for file in stefanfiles]

results_stefan = sorted(zip(stefan_sindec,stefan_inj,stefan_flux))
sindecrange_script,inj_script,flux_script = zip(*results_stefan)




################################
## Read in sensitivity calculations yielded from submitter ##
#files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]

## Extract flux and declination (which I get from the filenames) ##
#flux = [file['flux'][0] for file in files]
#inj = [file['mu'][0] for file in files]
#dec = [np.float(line.split('dec_',1)[-1].split('.array')[0]) for line in os.listdir(datafolder) if line.endswith('.array')]
#sindecrange=np.sin(np.radians(dec))

## I want them in the right order for plotting (to make connecting them easier and better looking). To do this, I'll sort a zipped list of sindecrange and flux, and sort it by sindec. ##

#results = sorted(zip(sindecrange,flux))
#sindecrange,flux = zip(*results)

##Now I'll load in the results that I got from mhuber's stacking code, but for single source sensitivities. Currently only one is finished, I'll have to come back to do this later. ##

#Here's a definition for getting sensitivities in the right format from the specified directories.

def plotSensitivity(datafolder_stack_bin):
    stack_files = [[cache.load(datafolder_stack_bin+decfolder+'/sensitivity/'+file) for file in os.listdir(misc.ensure_dir(datafolder_stack_bin+decfolder+'/sensitivity/')) if file.endswith('.array')] for decfolder in os.listdir(datafolder_stack_bin) if decfolder.startswith('dec')]

    ##Let's find the declinations used first.
    decfolders = [file for file in os.listdir(datafolder_stack_bin) if file.startswith('dec')]
    dec_stack =  [np.float(decfolder[3:]) for decfolder in decfolders]

    #We might not have gotten every value possible, so remove the blanks.

    for i in reversed(range(len(stack_files))):
      if stack_files[i] == []:
        stack_files.remove([])
        dec_stack.pop(i)


    flux_stack = [stack_file[0]['flux'] for stack_file in stack_files]
    inj_stack = [stack_file[0]['mu'] for stack_file in stack_files]
    sindecrange_stack=np.sin(np.radians(dec_stack))

    results_stack = sorted(zip(sindecrange_stack,inj_stack,flux_stack))
    sindecrange_stack,inj_stack,flux_stack = zip(*results_stack)

    return sindecrange_stack, inj_stack, flux_stack

#I'll plot regular stacking and no-energy stacking for now.
sindec_stack, inj_stack, flux_stack = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/')
sindec_no_energy, inj_no_energy, flux_no_energy = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/no_energy/')
sindec_double, inj_double, flux_double = plotSensitivity('/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/double_source/')

#And now for gamma=1,3 for regular sensitivity calculation.
sindec_gamma1, inj_gamma1, flux_gamma1 = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/gamma/1.0/')

sindec_gamma3, inj_gamma3, flux_gamma3 = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/gamma/3.0/')

##Also need to compare these results to jfeintzeig's IC86I sensitivities.##

filename4 = open('/data/user/brelethford/Data/jfeintzeig_sensitivity/4yr_allsky.csv','r')
rawdata4=filename4.read()
filename4.close()
table = [map(str, row.split()) for row in rawdata4.strip().split("\n")]
x4= [table[i][0].split(',')[0] for i in range(len(table))] #sindec
y4= [table[i][1].split(',')[0] for i in range(len(table))] #sensitivity

filename1 = open('/data/user/brelethford/Data/jfeintzeig_sensitivity/IC86I_upgoing.csv','r')
rawdata1=filename1.read()
filename1.close()
table = [map(str, row.split()) for row in rawdata1.strip().split("\n")]
x1= [table[i][0].split(',')[0] for i in range(len(table))] #sindec
y1= [table[i][1].split(',')[0] for i in range(len(table))] #sensitivity

#Let's include stefan's sensitivities#

datafolder_stefan = '/data/user/brelethford/Data/stefan_c_sensitivity/'
sindec_stefan, flux_stefan = np.genfromtxt (datafolder_stefan + 'orig_sens.txt').T
#Currently in GeV - change to TeV
flux_stefan /= 1e3

sindecmask=np.array(np.logical_and((sindec_stack>sindec_stefan[0]), (sindec_stack < sindec_stefan[-1])))
interp_flux=np.interp(np.array(sindec_stack)[sindecmask], sindec_stefan, flux_stefan)
#interp_flux is a list, where flux_stack is a tuple of lists of one element. Before plotting the ratio, I'll change it into a regular list.
## Sensitivity plot. ##
fig_sensitivity = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_stack,flux_stack, ':',label='Space+Energy llh-Model', color = 'red', linewidth=3)
plt.plot(sindec_no_energy,flux_no_energy, ':', label='Space Only llh-Model', color = 'green', linewidth = 3)
plt.plot(sindec_double,flux_double, ':', label='double source', color = 'blue', linewidth =3)
plt.plot(sindecrange_script,flux_script, '--', label='stefan - from script', color = 'k')
plt.plot(sindec_stefan,flux_stefan,label='stefan IC86 sensitivity', color='k')
ax.set_title(r'Single Source Sensitivity')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$E^2 \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
ax.set_xlim(-1.,1.)
ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid() 
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_sensitivity.savefig('/data/user/brelethford/AGN_Core/Plots/single_source_sensitivity.pdf')

#flux_stack in wrong format - change to a list, like interp.
flux_reformat=[]
for i in flux_stack:
  flux_reformat.append(i[0])

#Ratio Plot
fig_ratio = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(np.array(sindec_stack)[sindecmask],np.array(flux_reformat)[sindecmask]/interp_flux)
ax.set_title(r'Me vs. Stefan single source sensitivity ratio')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylim(.8,1.2)
ax.set_xlim(-1.,1.)
ax.grid()
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_ratio.savefig('/data/user/brelethford/AGN_Core/Plots/single_source_ratio.pdf')

#Injected Events Plot
fig_injected_events = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_stack,inj_stack, ':',label='Space+Energy llh-Model', color = 'red', linewidth = 3)
plt.plot(sindec_no_energy,inj_no_energy, ':', label='Space Only llh-Model', color = 'green', linewidth = 3)
plt.plot(sindec_double,inj_double, ':', label='double source', color = 'blue', linewidth = 3)
plt.plot(sindecrange_script,inj_script, '--', label='stefan - from script', color = 'k')
#plt.plot(x4,y4,label='jfeintzeig 4yr sensitivity', color='blue')
#Need to get stefan's injected events from somewhere...
#plt.plot(sindec_stefan,inj_stefan,label='stefan IC86 sensitivity', color='red')
ax.set_title(r'Single Source Injected Events')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$\mu_{inj}$')
ax.grid()
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_injected_events.savefig('/data/user/brelethford/AGN_Core/Plots/single_source_injected_events.pdf')

#Here I'm checking being able to do calculations with a different spectral index.
#Different Gamma Sensitivities:
fig_gamma = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_gamma1,flux_gamma1, '-',label=r'$\gamma=1$', color = 'k', alpha=.3)
plt.plot(sindec_stack,flux_stack, '-',label=r'$\gamma=2$', color = 'k', alpha=.6)
plt.plot(sindec_gamma3,flux_gamma3, '-',label=r'$\gamma=3$', color = 'k', alpha=1)
#plt.plot(sindec_stefan,flux_stefan,label='stefan IC86 sensitivity', color='k')
ax.set_title(r'$\gamma$-dependent Sensitivity')
ax.set_xlabel(r'$\sin{(\delta)}$')
#ax.set_ylabel(r'$E^2 \frac{E}{100TeV}^{\gamma-2} \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
ax.set_ylabel(r'$E^{ \gamma} \frac{d\phi}{dE} [\frac{TeV^{\gamma-1}}{cm^2 s}]$')
ax.set_xlim(-1.,1.)
#ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid() 
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_gamma.savefig('/data/user/brelethford/AGN_Core/Plots/gamma_sensitivity.pdf')



#######################################


##Okay, now let's do all the same stuff, but for discovery potentials.

def plotDisc(datafolder_stack_bin):
    stack_files = [[cache.load(datafolder_stack_bin+decfolder+'/disc/'+file) for file in os.listdir(misc.ensure_dir(datafolder_stack_bin+decfolder+'/disc/')) if file.endswith('.array')] for decfolder in os.listdir(datafolder_stack_bin) if decfolder.startswith('dec')]

    ##Let's find the declinations used first.
    decfolders = [file for file in os.listdir(datafolder_stack_bin) if file.startswith('dec')]
    dec_stack =  [np.float(decfolder[3:]) for decfolder in decfolders]

    #We might not have gotten every value possible, so remove the blanks.

    for i in reversed(range(len(stack_files))):
      if stack_files[i] == []:
        stack_files.remove([])
        dec_stack.pop(i)


    flux_stack = [stack_file[0]['flux'] for stack_file in stack_files]
    inj_stack = [stack_file[0]['mu'] for stack_file in stack_files]
    sindecrange_stack=np.sin(np.radians(dec_stack))

    results_stack = sorted(zip(sindecrange_stack,inj_stack,flux_stack))
    sindecrange_stack,inj_stack,flux_stack = zip(*results_stack)

    return sindecrange_stack, inj_stack, flux_stack

######### NOTE - I'm reassigning these values, just because it's easier. If I need them both on the same plot I'll do it again later.

sindec_stack, inj_stack, flux_stack = plotDisc( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/')
sindec_no_energy, inj_no_energy, flux_no_energy = plotDisc( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/no_energy/')
sindec_double, inj_double, flux_double = plotDisc('/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/double_source/')

#Plot discovery potential
fig_disc = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_stack,flux_stack, ':',label='Space+Energy llh-Model', color = 'red', linewidth=3)
plt.plot(sindec_no_energy,flux_no_energy, ':', label='Space Only llh-Model', color = 'green', linewidth = 3)
plt.plot(sindec_double,flux_double, ':', label='double source', color = 'blue', linewidth =3)
#plt.plot(sindecrange_script,flux_script, '--', label='stefan - from script', color = 'k')
#plt.plot(sindec_stefan,flux_stefan,label='stefan IC86 sensitivity', color='k')
ax.set_title(r'Single Source Discovery Potential')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$E^2 \frac{dN}{dE} [\frac{TeV}{cm^2 s}]$') 
ax.set_xlim(-1.,1.)
#ax.set_ylim(1e-12,1e-10)
ax.semilogy()
ax.grid() 
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_disc.savefig('/data/user/brelethford/AGN_Core/Plots/single_source_discovery_potential.pdf')
