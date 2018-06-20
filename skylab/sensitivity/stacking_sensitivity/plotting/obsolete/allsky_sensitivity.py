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
plt.plot(sindec_86,flux_86, label = '86I',color='red')
plt.plot(sindec_79,flux_79, label = '79',color='blue')
plt.plot(sindec_59,flux_59, label = '59',color='green')
plt.plot(sindec_40,flux_40, label = '40',color='purple')
ax.set_title(r'Single Year')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$E^2 \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$') 
ax.set_xlim(-1.,1.)
ax.set_ylim(1.5e-13,1e-9)
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

##Now let's do the same thing, this time for multiple years (1,2,3,4). Start with IC40, add one each time.

multistorage = '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/multi_year/'

#datafolder40 is the same as before, but for consistancy's sake let's call it 1year

datafolder1yr = datastorage + 'IC40/'
datafolder2yr = multistorage + '2/'
datafolder3yr = multistorage + '3/'
datafolder4yr = multistorage + '4/'

sindec_1yr, inj_1yr, flux_1yr, TS_1yr = plotyear(datafolder1yr)
sindec_2yr, inj_2yr, flux_2yr, TS_2yr = plotyear(datafolder2yr)
sindec_3yr, inj_3yr, flux_3yr, TS_3yr = plotyear(datafolder3yr)
sindec_4yr, inj_4yr, flux_4yr, TS_4yr = plotyear(datafolder4yr)

#interp_flux is a list, where flux_stack is a tuple of lists of one element. Before plotting the ratio, I'll change it into a regular list.
## Sensitivity plot. ##
fig_sensitivity = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_1yr,flux_1yr, label = '1yr: IC40',color='purple')
plt.plot(sindec_2yr,flux_2yr, label = '2yr: IC40-59',color='green')
plt.plot(sindec_3yr,flux_3yr, label = '3yr: IC40-79',color='blue')
plt.plot(sindec_4yr,flux_4yr, label = '4yr: IC40-86I',color='red')
ax.set_title(r'Cumulative')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$E^2 \frac{\partial N}{\partial E} [\frac{\textrm{\small{TeV}}}{\textrm{\small{cm}}^2 \textrm{\small{s}}}]$') 
ax.set_xlim(-1.,1.)
ax.set_ylim(1.5e-13,1e-9)
ax.semilogy()
ax.grid() 
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_sensitivity.savefig('/data/user/brelethford/AGN_Core/Plots/allsky_check/sensitivity_multiyear.pdf')

fig_injected_events = plt.figure(figsize=(w, .75*w))
ax=plt.gca()
plt.plot(sindec_1yr,inj_1yr, '--', label = '1yr: IC40',color='purple')
plt.plot(sindec_2yr,inj_2yr, '--', label = '2yr: IC40-59',color='green')
plt.plot(sindec_3yr,inj_3yr, '--', label = '3yr: IC40-79',color='blue')
plt.plot(sindec_4yr,inj_4yr, '--', label = '4yr: IC40-86I',color='red')

ax.set_title(r'Single Source Injected Events')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$\mu_{inj}$')
ax.grid()
plt.legend(loc='upper right', prop=propxxsmall)
plt.subplots_adjust (left=.2, bottom=.2)
fig_injected_events.savefig('/data/user/brelethford/AGN_Core/Plots/allsky_check/injected_events_multiyear.pdf')
