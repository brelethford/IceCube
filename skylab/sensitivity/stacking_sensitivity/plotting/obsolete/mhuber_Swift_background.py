import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from scipy.stats import chi2, norm, exponweib
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite

from skylab import statistics
fitfun = statistics.delta_chi2
weibfun=statistics.weib

#Note: this script is for plotting the old TS information - before bstacking. Treat accordingly.

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

def get_inj(folder):
    results = cache.load(folder+'sensitivity/gamma2.0.array')
    return results['TS']

datafolder = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/4yr/'

datafolder_uniform = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/4yr/uniform/'
datafolder_flux = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/4yr/flux/'
datafolder_redshift = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/4yr/redshift/'

background_uniform = datafolder_uniform+'background_trials/'
background_flux = datafolder_flux+'background_trials/'
background_redshift = datafolder_redshift+'background_trials/'

sensitivity_uniform = datafolder_uniform+'uniform_inj/'
sensitivity_flux = datafolder_flux+'flux_inj/'
sensitivity_redshift = datafolder_redshift+'redshift_inj/'


#Load in bckg and inj TSs
bckg_uniform = getBckg(background_uniform)
bckg_flux = getBckg(background_flux)
bckg_redshift = getBckg(background_redshift)

inj_uniform = get_inj(sensitivity_uniform)
inj_flux = get_inj(sensitivity_flux)
inj_redshift = get_inj(sensitivity_redshift)

print ('background-only trials: median TS')
print (bckg_uniform['TS_beta'])
print (bckg_flux['TS_beta'])
print (bckg_redshift['TS_beta'])

print ('background-only trials: max TS')
print max(bckg_uniform['TS'])
print max(bckg_flux['TS'])
print max(bckg_redshift['TS'])

print ('background-only trials: percentage of TS=0')
print str(1-np.float(np.count_nonzero(bckg_uniform['TS']))/np.float(len(bckg_uniform['TS'])))
print str(1-np.float(np.count_nonzero(bckg_flux['TS']))/np.float(len(bckg_flux['TS'])))
print str(1-np.float(np.count_nonzero(bckg_redshift['TS']))/np.float(len(bckg_redshift['TS'])))


print ('with signal trials: median TS')
print (np.median(inj_uniform))
print (np.median(inj_flux))
print (np.median(inj_redshift))

print ('with signal trials: max TS')
print max(inj_uniform)
print max(inj_flux)
print max(inj_redshift)

print ('with signal trials: percentage of TS=0')
print str(1-np.float(np.count_nonzero(inj_uniform))/np.float(len(inj_uniform)))
print str(1-np.float(np.count_nonzero(inj_flux))/np.float(len(inj_flux)))
print str(1-np.float(np.count_nonzero(inj_redshift))/np.float(len(inj_redshift)))

##Now we make hists of the test statistics ##
bins = 50
range = (0.0,40.0)
uniform_hist =histlite.hist(bckg_uniform['TS'],bins=bins,range=range)
redshift_hist =histlite.hist(bckg_redshift['TS'],bins=bins,range=range)
flux_hist =histlite.hist(bckg_flux['TS'],bins=bins,range=range)

uniform_inj_hist =histlite.hist(inj_uniform,bins=bins,range=range)
redshift_inj_hist =histlite.hist(inj_redshift,bins=bins,range=range)
flux_inj_hist =histlite.hist(inj_flux,bins=bins,range=range)

#I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chi2fit_uniform = fitfun(bckg_uniform['TS'],df=2., floc=0., fscale=1.)
chi2fit_redshift = fitfun(bckg_redshift['TS'],df=2., floc=0., fscale=1.)
chi2fit_flux = fitfun(bckg_flux['TS'],df=2., floc=0., fscale=1.)

#Now for a weibfit on the distribution.
weib_uniform = weibfun(bckg_uniform['TS'],df=2., floc=0., fscale=1.)
weib_redshift = weibfun(bckg_redshift['TS'],df=2., floc=0., fscale=1.)
weib_flux = weibfun(bckg_flux['TS'],df=2., floc=0., fscale=1.)

#In addition, let's get the completed with-signal distributions.


## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

chi2_fits = [chi2fit_uniform,chi2fit_redshift,chi2fit_flux]
weib_fits = [weib_uniform,weib_redshift,weib_flux]
colors=['blue','red','green']
#for chi2_fit, weib_fit, color in zip(chi2_fits, weib_fits, colors):
#  x = np.linspace(0,20,100)
#  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color, label=r'$\tilde{\chi}^2$')#: df='+str(round(chi2_fit.par[0],2)))

for chi2_fit, weib_fit, color in zip(chi2_fits, weib_fits, colors):
  #x = np.linspace(0,20,100)
  #ax.plot(x, weib_fit.pdf(x), linestyle='--', color=color, label = 'weibull')
  plt.axvline(chi2_fit.isf(norm.sf(0)), color = color, linestyle = ':')
  print (chi2_fit)
  #plt.axvline(weib_fit.isf(norm.sf(5)), color = color, linestyle = '--')
  
#Finally plot the background TS
histlite.plot1d(ax,uniform_hist.normalize(integrate=True),histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,redshift_hist.normalize(integrate=True),histtype='step',label='redshift',color='red')
histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',label='flux',color='green')

#And the signal TS
histlite.plot1d(ax,uniform_inj_hist.normalize(integrate=True),histtype='step', alpha = 0.5, color='blue')
histlite.plot1d(ax,redshift_inj_hist.normalize(integrate=True),histtype='step',alpha = 0.5, color='red')
histlite.plot1d(ax,flux_inj_hist.normalize(integrate=True),histtype='step', alpha = 0.5, color='green')

#And include the preliminary thing.

def ic_prelim(fig, x = 0.75, y = 0.8, **kw):
    """Marks maps and plots as preliminary"""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text(x, y, "IceCube Preliminary", **kw)

#ax.set_title(r'70m TS distributions - 4yr (40-86I)')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ic_prelim(fig_bckg, x = 0.5, y=0.6)
ax.set_ylabel(r'Probability Density') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,40)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall, ncol=1)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/mhuberTS_70m_4yr.pdf')
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/mhuberTS_70m_4yr.png')

