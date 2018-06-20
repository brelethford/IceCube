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
weibfun = statistics.weib
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

datafolder_2LAC = '/data/user/brelethford/Output/stacking_sensitivity/2LAC/flux_3yr/background_trials/'

bckg_2LAC = getBckg(datafolder_2LAC)

print (bckg_2LAC['TS_beta'])


##Now we make hists of the test statistics ##
bins = 100
range = (0.0,20.0)
h_2LAC = histlite.hist(bckg_2LAC['TS'],bins=bins,range=range)

##I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chi2fit_2LAC = fitfun(bckg_2LAC['TS'],df=2., floc=0., fscale=1.)

#Now for a weibfit on the distribution.
weib_flux_2LAC = weibfun(bckg_2LAC['TS'],df=2., floc=0., fscale=1.)

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

weib_fits = [weib_flux_2LAC]
chi2_fits = [chi2fit_2LAC]
colors=['blue']
for chi2_fit,weib_fit,color in zip(chi2_fits,weib_fits,colors):
  x = np.linspace(0,20,100)
  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color, label=r'$\tilde{\chi}^2$')#: df='+str(round(chi2_fit.par[0],2)))
  ax.plot(x, weib_fit.pdf(x), linestyle='--', color=color, label = 'weibull')

histlite.plot1d(ax,h_2LAC.normalize(integrate=True),histtype='step',label='flux-weighted',color='blue')
ax.set_title(r'2LAC Background TS - 3yr (59-86I)')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall, ncol=2)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_2LAC_3yr.pdf')
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_2LAC_3yr.png')
