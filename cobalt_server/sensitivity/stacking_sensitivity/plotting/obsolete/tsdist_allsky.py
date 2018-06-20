import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import chi2
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite

from skylab import statistics
fitfun = statistics.delta_chi2

##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

###This script imports the sensitivities from the submitter and plots them.###
#picklefolder = '/data/user/brelethford/Data/SwiftBAT70m/pickle/'

## Define fcn to read in background trials previously logged ##

def getBckg(datafolder):
    files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
    n_inj=[]
    nsources=[]
    TS=[]
    beta=(0.5) #For background ts
    TS_beta=[] #Calculated from the total TS median after we get all the TS.
    beta_err=[]
    gamma=[]
    for file in files:
      for item in range(len(file['n_inj'])):
        n_inj.append(file['n_inj'][item])
        nsources.append(file['nsources'][item])
        TS.append(file['TS'][item])
        gamma.append(file['gamma'][item])

    TSs=TS
    TS_beta = np.percentile(TSs, 100.*(1. - beta))
    m=np.count_nonzero(np.asarray(TSs) > (TS_beta))
    i = len(TSs)
    fraction = float(m)/float(i)
    beta_err = (np.sqrt(fraction * (1. - fraction) / float(i)) if 0 < beta < 1 else 1.)
    bckg_trials = {'n_inj':n_inj,'nsources':np.asarray(nsources), 'TS':np.asarray(TS), 'beta':beta, 'beta_err':beta_err, 'TS_beta':TS_beta, 'gamma':np.asarray(gamma)}
    return bckg_trials

datafolder_79 = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/background_trials/'


bckg_uniform = getBckg(datafolder_uniform)

TSsens = np.median(bckg_uniform['TS'])
TSdisc = np.percentile(bckg_uniform['TS'], (1-(2.867e-7))*100)
##Now we make hists of the test statistics ##
bins = 80
range = (0.0,30.0)
uniform_hist =histlite.Hist.normalize(histlite.hist(bckg_uniform['TS'],bins=bins,range=range))

#I'll do the same thing for the injected events to achieve sensitivity and discovery potential for uniform llh uniform inj.

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()


histlite.plot1d(ax,uniform_hist,histtype='step',label='background',color='blue')
histlite.plot1d(ax,senshist,histtype='step',label='sensitivity: n inj = {0:.2f}'.format(sensinj[0]) ,color='red')
histlite.plot1d(ax,dischist,histtype='step',label='discovery potential: n inj = {0:.2f}'.format(discinj[0]) ,color='green')

plt.axvline(TSsens, linestyle='--', color = 'k')
plt.axvline(TSdisc, linestyle='--', color = 'k')
ax.set_title(r'TS distribution - example')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,30)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/ts_example_ninj.pdf')
