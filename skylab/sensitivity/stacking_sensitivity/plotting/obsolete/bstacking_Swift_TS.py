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

##let's define weibfun here:


##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

###This script imports the sensitivities from the submitter and plots them.###
#picklefolder = '/data/user/brelethford/Data/SwiftBAT70m/pickle/'

## Define fcn to read in background trials previously logged ##




datafolder_uniform = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/uniform/background_trials/'
datafolder_flux = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/flux/background_trials/'
datafolder_redshift = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/redshift/background_trials/'

def get_TS(datafolder):
  files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
  TS=[]
  for file in files:
    for entry in file:
      if entry[0]==0:
        TS.append(entry[1])
  TSs=TS
  chi2fit = fitfun(TSs, df=2., floc=0., fscale=1.)
  weib = weibfun(TSs,df=2., floc=0., fscale=1.)
  return TSs, chi2fit, weib

TS_uniform,chi2fit_uniform,weib_uniform = get_TS(datafolder_uniform)
TS_flux,chi2fit_flux,weib_flux = get_TS(datafolder_flux)
TS_redshift,chi2fit_redshift,weib_redshift = get_TS(datafolder_redshift)

weighteach = ['uniform','flux','redshift']
TSeach = [TS_uniform,TS_flux,TS_redshift]
chi_each = [chi2fit_uniform,chi2fit_flux,chi2fit_redshift]
weib_each = [weib_uniform,weib_flux,weib_redshift]
for weight, TS, chi2fit in zip(weighteach,TSeach,chi_each):
  print('For ' + str(weight) +', the median TS is: ' + str(np.asscalar(chi2fit.isf(0.5))))
  print('               max TS is: ' + str(max(TS)))
  print('               percentage of TS = 0 is: ') + str(1.0-np.float(np.count_nonzero(TS))/np.float(len(TS)))

##Now we make hists of the test statistics ##
bins = 50
range = (0.0,40.0)
uniform_hist =histlite.hist(TS_uniform,bins=bins,range=range)
flux_hist =histlite.hist(TS_flux,bins=bins,range=range)
redshift_hist =histlite.hist(TS_redshift,bins=bins,range=range)

fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

colors=['blue','green','red']
for chi2_fit, weib_fit, color in zip(chi_each, weib_each, colors):
  x = np.linspace(0,20,100)
  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color, label=r'$\tilde{\chi}^2$')#: df='+str(round(chi2_fit.par[0],2)))
  plt.axvline(np.asscalar(chi2_fit.isf(0.5)),color=color)
#for chi2_fit, weib_fit, color in zip(chi_each, weib_each, colors):
#  ax.plot(x, weib_fit.pdf(x), linestyle='--', color=color, label = 'weibull')
#  plt.axvline(chi2_fit.isf(norm.sf(5)), color = color, linestyle = ':')
#  plt.axvline(weib_fit.isf(norm.sf(5)), color = color, linestyle = '--')

#Now we have to get the TS from the injected trials

inj_uniform = cache.load('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/uniform/uniform_inj/sensitivity/gamma2.0.array')
inj_flux = cache.load('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/flux/flux_inj/sensitivity/gamma2.0.array')
inj_redshift = cache.load('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/bstacking_4yr/redshift/redshift_inj/sensitivity/gamma2.0.array')

def getTSinj(data):
    TS = []
    for trial in data['trials']:
        TS.append(trial[1])
    return TS

uniform_inj_hist =histlite.hist(getTSinj(inj_uniform),bins=bins,range=range)
flux_inj_hist =histlite.hist(getTSinj(inj_flux),bins=bins,range=range)
redshift_inj_hist =histlite.hist(getTSinj(inj_redshift),bins=bins,range=range)
#background
histlite.plot1d(ax,uniform_hist.normalize(integrate=True),histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',label='flux',color='green')
histlite.plot1d(ax,redshift_hist.normalize(integrate=True),histtype='step',label='redshift',color='red')
#And with signal
histlite.plot1d(ax,uniform_inj_hist.normalize(integrate=True),histtype='step', alpha = 0.5, color='blue')
histlite.plot1d(ax,redshift_inj_hist.normalize(integrate=True),histtype='step',alpha = 0.5, color='red')
histlite.plot1d(ax,flux_inj_hist.normalize(integrate=True),histtype='step', alpha = 0.5, color='green')
              
def ic_prelim(fig, x = 0.75, y = 0.8, **kw):
    """Marks maps and plots as preliminary"""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text(x, y, "IceCube Preliminary", **kw)                                                             
#ax.set_title(r'70m Background TS - 4yr (40-86I)')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Probability Density')
ax.set_ylim(8e-5,1.0)
ic_prelim(fig_bckg, x = 0.5, y=0.6)
ax.set_xlim(0,40)
ax.semilogy()
plt.legend(loc='upper right', prop=propxsmall, ncol=2)

fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bstacking_bckgTS_70m_4yr.pdf')
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/backgroundTS/bstacking_bckgTS_70m_4yr.png')

