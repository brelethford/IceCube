import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
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

##Now we make hists of the test statistics ##
#The choice in range is to sidestep the issue that most of the bins are at TS=0, which screws up our fit.
bins = 100-1
range = (20./(bins+1),20.0)

uniform_hist =histlite.hist(bckg_uniform['TS'],bins=bins,range=range)
redshift_hist =histlite.hist(bckg_redshift['TS'],bins=bins,range=range)
flux_hist =histlite.hist(bckg_flux['TS'],bins=bins,range=range)
#bincenters will be the same for all of these hists
binedges = uniform_hist.bins[0]
x = [np.mean([binedges[i+1],binedges[i]]) for i in np.arange(bins)]

#I'll try epinat's method of doing this.
class my_weibull(stats.rv_continuous):
  def _cdf(self,x,a,b):
    return np.where(x >0,1-np.exp(-(x/a)**b),0)

my_fittingfunc = my_weibull()

#Now let's histogram the cumulative
##Let's normalize first - I'm going to have to see the cumulative distribution separately.
uniform_norm = uniform_hist.normalize(integrate=False)
redshift_norm = redshift_hist.normalize(integrate=False)
flux_norm = flux_hist.normalize(integrate=False)

###Let's try fitting the cumulative sum to a weibull distribution
#first, define the cumulative sum for these hists.

uniform_cumsum = uniform_norm.cumsum()
flux_cumsum = flux_norm.cumsum()
redshift_cumsum = redshift_norm.cumsum()

#Now let's fit them to an exponential weibull function.

def weibull(x,a,b):
    return 1 - np.exp(-(x/a)**b)

pars_all=[]

## Now to plot. ##
fig_cumsum = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

popt_unif, pcov_unif = scipy.optimize.curve_fit(weibull, x, uniform_cumsum.values)
popt_flux, pcov_flux = scipy.optimize.curve_fit(weibull, x, flux_cumsum.values)
popt_redshift, pcov_redshift = scipy.optimize.curve_fit(weibull, x, redshift_cumsum.values)

my_fittingfunc = my_weibull()   

ax.plot(x, my_fittingfunc.cdf(x,popt_unif[0], popt_unif[1]), color="black")
ax.plot(x, my_fittingfunc.cdf(x,popt_flux[0], popt_flux[1]), color="black")
ax.plot(x, my_fittingfunc.cdf(x,popt_redshift[0], popt_redshift[1]), color="black")
ax.plot(x, weibull(x,*popt_unif), "--", color = "magenta")
ax.plot(x, weibull(x,*popt_flux), "--", color = "magenta")
ax.plot(x, weibull(x,*popt_redshift), "--", color = "magenta")

histlite.plot1d(ax,uniform_norm.cumsum(),histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,redshift_norm.cumsum(),histtype='step',label='redshift',color='red')
histlite.plot1d(ax,flux_norm.cumsum(),histtype='step',label='flux',color='green')
ax.set_title(r'70m Bckg TS $\textgreater$ 0 - 3yr (59-86I) - cumsum')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_xlim(0,10)
ax.set_ylim(0,1)
plt.legend(loc='upper right', prop=propxsmall, ncol=1)
fig_cumsum.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr_cumsum.pdf')
fig_cumsum.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr_cumsum.png')

##Okay, this looks pretty good. Now let's plot the actual pdf using the fits from the cumsum.

#Same as before, but now we don't worry about TS = 0 (we already have our fit for the tail, so let's just follow that.)
'''
bins = 100
range = (0.,20.0)

uniform_hist =histlite.hist(bckg_uniform['TS'],bins=bins,range=range)
redshift_hist =histlite.hist(bckg_redshift['TS'],bins=bins,range=range)
flux_hist =histlite.hist(bckg_flux['TS'],bins=bins,range=range)
#bincenters will be the same for all of these hists
binedges = uniform_hist.bins[0]
x = [np.mean([binedges[i+1],binedges[i]]) for i in np.arange(bins)]

uniform_norm = uniform_hist.normalize(integrate=False)
redshift_norm = redshift_hist.normalize(integrate=False)
flux_norm = flux_hist.normalize(integrate=False)
'''
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

ax.plot(x, my_fittingfunc.pdf(x,popt_unif[0], popt_unif[1]), color="blue", linestyle = '--')
ax.plot(x, my_fittingfunc.pdf(x,popt_flux[0], popt_flux[1]), color="green", linestyle = '--')
ax.plot(x, my_fittingfunc.pdf(x,popt_redshift[0], popt_redshift[1]), color="red", linestyle = '--')
#ax.plot(x, weibull(x,*popt_unif), "--", color = "magenta")
#ax.plot(x, weibull(x,*popt_flux), "--", color = "magenta")
#ax.plot(x, weibull(x,*popt_redshift), "--", color = "magenta")

histlite.plot1d(ax,uniform_hist.normalize(integrate=True),histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',label='flux',color='green')
histlite.plot1d(ax,redshift_hist.normalize(integrate=True),histtype='step',label='redshift',color='red')


#histlite.plot1d(ax,redshift_norm.cumsum(),histtype='step',label='redshift',color='red')
#histlite.plot1d(ax,flux_norm.cumsum(),histtype='step',label='flux',color='green')
ax.set_title(r'70m Bckg TS - 3yr (59-86I)')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy()
plt.legend(loc='upper right', prop=propxsmall, ncol=1)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr_weib.pdf')
#fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bckgTS_70m_3yr_weib.png')





'''


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
'''
