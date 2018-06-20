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

alldata = '/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/bstacking_3yr/flux/'

bckg_sirin = alldata + 'background_trials/'
bckg_schatto = alldata + 'schatto79_background_trials/'
sens_sirin = alldata + 'flux_inj/sensitivity/'
sens_schatto = alldata + 'schatto79_flux_inj/sensitivity/'
disc_sirin = alldata + 'flux_inj/discovery/'
disc_schatto = alldata + 'schatto79_flux_inj/discovery/'

bins = 40
range = (0.0,40.0)

def makehist(datafolder):
  files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
  TS=[]
  for file in files:
    for entry in file:
      if entry[0]==0:
        if entry[1]<0:
          TS.append(0)
        else:
          TS.append(entry[1])
  TSs=TS
  chi2fit_flux = fitfun(TSs, df=1., floc=0., fscale=1.)
  print('background only: median TS = {}'.format(str(np.asscalar(chi2fit_flux.isf(0.5)))))
  print ('max TS = {}'.format(str(max(TSs))))
  print ('Percentage of TS=0 is: {}'.format(str(1.0-np.float(np.count_nonzero(TSs))/np.float(len(TSs)))))
  ##Now we make hists of the test statistics ##
  flux_hist =histlite.hist(TSs,bins=bins,range=range)
  return flux_hist,chi2fit_flux

hist_sirin,chi_sirin = makehist(bckg_sirin)
hist_schatto,chi_schatto = makehist(bckg_schatto)
def makehist_inj(datafolder):
  files = cache.load(datafolder+'gamma2.0.array')

  TSs=files['trials']['TS']

  chi2fit_flux = fitfun(TSs, df=1., floc=0., fscale=1.)
  weib_flux = weibfun(TSs,df=2., floc=0., fscale=1.)
  print('background plus injected: median TS = {}'.format(str(np.asscalar(chi2fit_flux.isf(0.5)))))
  print ('max TS = {}'.format(str(max(TSs))))
  print ('Percentage of TS=0 is: {}'.format(str(1.0-np.float(np.count_nonzero(TSs))/np.float(len(TSs)))))
  
  ##Now we make hists of the test statistics ##
  flux_hist =histlite.hist(TSs,bins=bins,range=range)
  return flux_hist,chi2fit_flux

def printflux(datafolder):
  files = cache.load(datafolder+'gamma2.0.array')
  flux= files['flux']
  mu = files['mu']
  TSval = files['TSval']
  return (flux, mu, TSval)

#extract hists for injected trials
hist_sirin_sens,chi_sirin_sens = makehist_inj(sens_sirin)
hist_schatto_sens,chi_schatto_sens = makehist_inj(sens_schatto)
hist_sirin_disc,chi_sirin_disc = makehist_inj(disc_sirin)
hist_schatto_disc,chi_schatto_disc = makehist_inj(disc_schatto)

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

labels = ['sirin', 'schatto']
flux_hists = [hist_sirin,hist_schatto]
chi2_fits = [chi_sirin,chi_schatto]

sens_hists = [hist_sirin_sens,hist_schatto_sens]
sens_fits = [chi_sirin_sens,chi_schatto_sens]
disc_hists = [hist_sirin_disc,hist_schatto_disc]
disc_fits = [chi_sirin_disc,chi_schatto_disc]

colors=['green', 'blue']
#First for bckg
for flux_hist, chi2_fit, color,label in zip(flux_hists,chi2_fits, colors,labels):
  x = np.linspace(0,40,100)
  histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',color=color,label='{} - '.format(label) + r'$\tilde{\chi}^2$' + ': df={}'.format(str(round(chi2_fit.par[0],2))))
  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color)#, label='{} - '.format(label) + r'$\tilde{\chi}^2$' + ': df={}'.format(str(round(chi2_fit.par[0],2))))
#now for sens, then disc 
for flux_hist, chi2_fit, color in zip(sens_hists,sens_fits, colors):
  histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',color=color)
  #plt.axvline(chi2_fit.isf(0.1), color = color, linestyle = ':')

for flux_hist, chi2_fit, color in zip(disc_hists,disc_fits, colors):
  histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',color=color)
  #plt.axvline(chi2_fit.isf(0.5), color = color, linestyle = ':')

###NOTE - don't use the following two lines, they don't' give a good chi2 fit for some reason. I think we'd specifically need to use the deltachi2 function. 
#y = np.linspace(1e-5,20,100) #nonzero start to prevent infinity
#ax.plot(y, chi2.pdf(y,1), linestyle=':',color='k', label=r'$\tilde{\chi}^2$: df=1')
#for chi2_fit, weib_fit, color in zip(chi2_fits, weib_fits, colors):
#  x = np.linspace(0,20,100)
#  ax.plot(x, weib_fit.pdf(x), linestyle='--', color=color, label = 'weibull')
#  plt.axvline(chi2_fit.isf(norm.sf(5)), color = color, linestyle = ':')
#  plt.axvline(weib_fit.isf(norm.sf(5)), color = color, linestyle = '--')

#Now let's get the signal

#inj_flux = cache.load('/data/user/brelethford/Output/stacking_sensitivity/4yr_Starburst/4yr/flux_mhuber_git/flux_inj/sensitivity/gamma2.0.array')

#def getTSinj(data):
#    TS = []
#    for trial in data['trials']:
#        TS.append(trial[1])
#    return TS

#flux_inj_hist =histlite.hist(getTSinj(inj_flux),bins=bins,range=range)

#histlite.plot1d(ax,flux_inj_hist.normalize(integrate=True),histtype='step', alpha = 0.5, color='green')

def ic_prelim(fig, x = 0.75, y = 0.8, **kw):
    """Marks maps and plots as preliminary"""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text(x, y, "IceCube Preliminary", **kw)


#ax.set_title(r'4yr Starburst TS - updated')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ic_prelim(fig_bckg, x = 0.5, y=0.3)
ax.set_ylabel(r'Probability Density') 
ax.set_ylim(8e-5,1.2)
ax.set_xlim(0,40)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall, ncol=1)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/starburst/bstacking_starburst3yr_bckgTS.pdf')

print ( 'sens - sirin: flux, mu, TSval - ' + str(printflux(sens_sirin)))
print ( 'sens - schatto: flux, mu, TSval - ' + str(printflux(sens_schatto)))
print( 'disc - sirin: flux, mu, TSval - ' +str(printflux(disc_sirin)))
print( 'disc - schatto: flux, mu, TSval - ' +str(printflux(disc_schatto)))

