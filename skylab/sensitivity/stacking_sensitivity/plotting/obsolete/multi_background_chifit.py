import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
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

datafolder_uniform = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/uniform/background_trials/'
datafolder_flux = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/flux/background_trials/'
datafolder_redshift = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/redshift/background_trials/'

bckg_uniform = getBckg(datafolder_uniform)
bckg_flux = getBckg(datafolder_flux)
bckg_redshift = getBckg(datafolder_redshift)

print (bckg_uniform['TS_beta'])
print (bckg_flux['TS_beta'])
print (bckg_redshift['TS_beta'])


#bckg_uniform = cache.load(picklefolder+'bckg_trials_uniform.pickle')
#bckg_redshift = cache.load(picklefolder+'bckg_trials_redshift.pickle')
#bckg_flux = cache.load(picklefolder+'bckg_trials_flux.pickle')

##Now we make hists of the test statistics ##
bins = 1000
range = (0.0,20.0)
uniform_hist =histlite.hist(bckg_uniform['TS'],bins=bins,range=range)
redshift_hist =histlite.hist(bckg_redshift['TS'],bins=bins,range=range)
flux_hist =histlite.hist(bckg_flux['TS'],bins=bins,range=range)

#I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chi2fit_uniform = fitfun(bckg_uniform['TS'],df=2., floc=0., fscale=1.)
chi2fit_redshift = fitfun(bckg_redshift['TS'],df=2., floc=0., fscale=1.)
chi2fit_flux = fitfun(bckg_flux['TS'],df=2., floc=0., fscale=1.)

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

chi2_fits = [chi2fit_uniform,chi2fit_redshift,chi2fit_flux]
colors=['blue','red','green']
for chi2_fit,color in zip(chi2_fits,colors):
  x = np.linspace(0,20,100)
  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color, label=r'$\tilde{\chi}^2$: df='+str(round(chi2_fit.par[0],2)))

histlite.plot1d(ax,uniform_hist.normalize(integrate=True),histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,redshift_hist.normalize(integrate=True),histtype='step',label='redshift',color='red')
histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',label='flux',color='green')
ax.set_title(r'70m Background TS - 3yr (59-86I)')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall, ncol=2)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr.pdf')
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr.png')
