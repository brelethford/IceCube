import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import sys
from icecube.umdtools import cache
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
###This script imports the sensitivities from the submitter and plots them.###

datafolder = '/data/user/brelethford/Output/all_sky_sensitivity/results/n_iter/' ##Remember to change this depending on which results you're looking at!##

## Read in sensitivity calculations yielded from submitter ##
files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]

## Extract flux and declination (which I get from the filenames) ##
flux = [file['flux'][0] for file in files]
dec = [np.float(line.split('dec_',1)[-1].split('.array')[0]) for line in os.listdir(datafolder) if line.endswith('.array')]
sindecrange=np.sin(np.radians(dec))

## I want them in the right order for plotting (to make connecting them easier and better looking). To do this, I'll sort a zipped list of sindecrange and flux, and sort it by sindec. ##

results = sorted(zip(sindecrange,flux))
sindecrange,flux = zip(*results)

#Also need to compare these results to jfeintzeig's IC86I sensitivities.##

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


## Now to plot. ##
fig_sensitivity = plt.figure()
ax=plt.gca()
plt.plot(sindecrange,flux, label='sensitivity')
plt.plot(x4,y4, label='jfeintzeig 4yr sensitivity')
plt.plot(x1,y1, label='jfeintzeig IC86 sensitivity', color='red')
ax.set_title(r'Single Source Sensitivity - n_iter')
ax.set_xlabel(r'$\sin{(\delta)}$')
ax.set_ylabel(r'$\Phi_{0}[\frac{1}{GeV cm^2 s}]$') 
ax.semilogy() 
plt.legend(loc='upper right')
fig_sensitivity.savefig('/data/user/brelethford/AGN_Core/Plots/single_source_sensitivity_n_iter.pdf')

