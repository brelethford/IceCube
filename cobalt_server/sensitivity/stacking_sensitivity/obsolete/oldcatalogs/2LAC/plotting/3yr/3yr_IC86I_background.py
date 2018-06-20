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

datafolder_flux = '/data/user/brelethford/Output/stacking_sensitivity/2LAC/flux_3yr/background_trials/'

bckg_flux = getBckg(datafolder_flux)

print (bckg_flux['TS_beta'])


#bckg_uniform = cache.load(picklefolder+'bckg_trials_uniform.pickle')
#bckg_redshift = cache.load(picklefolder+'bckg_trials_redshift.pickle')
#bckg_flux = cache.load(picklefolder+'bckg_trials_flux.pickle')

##Now we make hists of the test statistics ##
bins = 80
range = (0.0,20.0)
flux_hist =histlite.Hist.normalize(histlite.hist(bckg_flux['TS'],bins=bins,range=range))

##I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chifit_flux = fitfun(bckg_flux['TS'],df=2., floc=0., fscale=1.).par[0]

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()


chi_degs = [1,2,chifit_flux]
colors=['black','gray','green']
for df,color in zip(chi_degs,colors):
  x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99999,df), 100)
  rv = chi2(df)
  chi_dist=rv.pdf(x)
  ax.plot(x, chi_dist/sum(chi_dist), linestyle=':',color=color, label=r'$\tilde{\chi}^2$: df='+str(round(df,2)))

histlite.plot1d(ax,flux_hist,histtype='step',label='flux',color='green')
ax.set_title(r'Background TS for 2LAC Survey - 3yr')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/2LAC/bckg_TS_2LAC_3yr.pdf')
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/2LAC/bckg_TS_2LAC_3yr.png')
