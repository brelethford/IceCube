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


#Here's a definition for getting sensitivities in the right format from the specified directories.

#def plotSensitivity(datafolder_stack_bin):
#    stack_files = [[cache.load(datafolder_stack_bin+decfolder+'/sensitivity/'+file) for file in os.listdir(misc.ensure_dir(datafolder_stack_bin+decfolder+'/sensitivity/')) if file.endswith('.array')] for decfolder in os.listdir(datafolder_stack_bin) if decfolder.startswith('dec')]
#
#    ##Let's find the declinations used first.
#    decfolders = [file for file in os.listdir(datafolder_stack_bin) if file.startswith('dec')]
#    dec_stack =  [np.float(decfolder[3:]) for decfolder in decfolders]
#
#    #We might not have gotten every value possible, so remove the blanks.
#
#    for i in reversed(range(len(stack_files))):
#      if stack_files[i] == []:
#        stack_files.remove([])
#        dec_stack.pop(i)
#
#
#    flux_stack = [stack_file[0]['flux'] for stack_file in stack_files]
#    inj_stack = [stack_file[0]['mu'] for stack_file in stack_files]
#    sindecrange_stack=np.sin(np.radians(dec_stack))
#    TS_stack = [stack_file[0]['TSval'][0] for stack_file in stack_files]
#
#    results_stack = sorted(zip(sindecrange_stack,inj_stack,flux_stack, TS_stack))
#    sindecrange_stack,inj_stack,flux_stack, TS_stack = zip(*results_stack)
#
#    return sindecrange_stack, inj_stack, flux_stack, TS_stack

#I'll plot regular stacking and no-energy stacking for now.
#sindec_stack, inj_stack, flux_stack, TS_stack = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/')
#sindec_no_energy, inj_no_energy, flux_no_energy, TS_no_energy = plotSensitivity( '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/no_energy/')
#sindec_double, inj_double, flux_double, TS_double = plotSensitivity('/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/double_source/')
#I grouped these ones differently, so let's redo a function for them alone.

def plotyear(datafolder):
    stack_files = [[cache.load(datafolder+decfolder+'/'+file) for file in os.listdir(datafolder+decfolder) if file.startswith('sens')] for decfolder in os.listdir(datafolder) if decfolder.startswith('dec')]
    ##Let's remove the files that didn't give back sensitivities first.
    decfolders = [file for file in os.listdir(datafolder) if file.startswith('dec')]
    dec_stack =  [np.float(decfolder[3:]) for decfolder in decfolders]

    for i in reversed(range(len(stack_files))):
      if stack_files[i] == []:
        stack_files.remove([])
        dec_stack.pop(i)

    flux_stack = [stack_files[index][0]['flux'] for index in range(len(stack_files))]
    inj_stack = [stack_files[index][0]['mu'] for index in range(len(stack_files))]
    sindecrange_stack=np.sin(np.radians(dec_stack))
    TS_stack = [stack_files[index][0]['TSval'] for index in range(len(stack_files))]
    results_stack = sorted(zip(sindecrange_stack,inj_stack,flux_stack,TS_stack))
    sindecrange_stack,inj_stack,flux_stack,TS_stack = zip(*results_stack)
    return sindecrange_stack, inj_stack, flux_stack, TS_stack
#    stack_files = [[cache.load(datafolder_stack_bin+decfolder+'/sensitivity/'+file) for file in os.listdir(misc.ensure_dir(datafolder_stack_bin+decfolder+'/sensitivity/')) if file.endswith('.array')] for decfolder in os.listdir(datafolder_stack_bin) if decfolder.startswith('dec')]

datastorage = '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/single_year/'
datafolder86 = datastorage + 'IC86/'
datafolder79 = datastorage + 'IC79/'
datafolder59 = datastorage + 'IC59/'
datafolder40 = datastorage + 'IC40/'
sindec_86, inj_86, flux_86, TS_86 = plotyear(datafolder86)
sindec_79, inj_79, flux_79, TS_79 = plotyear(datafolder79)
sindec_59, inj_59, flux_59, TS_59 = plotyear(datafolder59)
sindec_40, inj_40, flux_40, TS_40 = plotyear(datafolder40)

#interp_flux is a list, where flux_stack is a tuple of lists of one element. Before plotting the ratio, I'll change it into a regular list.
## Sensitivity plot. ##
fig_sensitivity = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
#plt.plot(sindec_stack,flux_stack, '--',label='Space+Energy llh-Model', color = 'red')
#plt.plot(sindec_no_energy,flux_no_energy, '--', label='Space Only llh-Model', color = 'green')
#plt.plot(sindec_double,flux_double, '--', label='Double Source', color = 'blue')
#plt.plot(sindecrange_script,flux_script, '--', label='stefan - from script', color = 'k')
plt.plot(sindec_86,flux_86, label = '86I',color='red')
plt.plot(sindec_79,flux_79, label = '79',color='blue')
plt.plot(sindec_59,flux_59, label = '59',color='green')
plt.plot(sindec_40,flux_40, label = '40',color='purple')
ax.set_title(r'Single Source Sensitivity')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$E^2 \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$') 
ax.set_xlim(-1.,1.)
ax.set_ylim(1e-12,1e-8)
ax.semilogy()
ax.grid() 
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_sensitivity.savefig('/data/user/brelethford/AGN_Core/Plots/allsky_check/sensitivity_1yr.pdf')

fig_injected_events = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
#plt.plot(sindec_stack,inj_stack, ':',label='Space+Energy llh-Model', color = 'red', linewidth = 3)
#plt.plot(sindec_no_energy,inj_no_energy, ':', label='Space Only llh-Model', color = 'green', linewidth = 3)
#plt.plot(sindec_double,inj_double, ':', label='double source', color = 'blue', linewidth = 3)
plt.plot(sindec_86,inj_86, '--', label='86I', color = 'red')
plt.plot(sindec_79,inj_79, '--', label='79', color = 'blue')
plt.plot(sindec_59,inj_59, '--', label='59', color = 'green')
plt.plot(sindec_40,inj_40, '--', label='40', color = 'purple')
ax.set_title(r'Single Source Injected Events')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$\mu_{inj}$')
ax.grid()
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_injected_events.savefig('/data/user/brelethford/AGN_Core/Plots/allsky_check/injected_events_1yr.pdf')

#fig_TS_med = plt.figure(figsize=(w, .75*w))
#ax=plt.gca()
#plt.plot(sindec_stack,TS_stack, ':',label='Space+Energy llh-Model', color = 'red', linewidth = 3)
#plt.plot(sindec_no_energy,TS_no_energy, ':', label='Space Only llh-Model', color = 'green', linewidth = 3)
#plt.plot(sindec_double,TS_double, ':', label='double source', color = 'blue', linewidth = 3)
#plt.plot(sindec_retry,TS_retry, '--', label='retry', color = 'k')
#ax.set_title(r'Single Source background TS median')
#ax.set_xlabel(r'$\sin{(\delta)}$')
#ax.set_ylabel(r'TS')
#ax.grid()
#plt.legend(loc='upper left', prop=propxxsmall)
#plt.subplots_adjust (left=.2, bottom=.2)
#fig_TS_med.savefig('/data/user/brelethford/AGN_Core/Plots/TS_median.pdf')
