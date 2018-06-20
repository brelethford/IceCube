import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite
##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

###This script imports the sensitivities from the submitter and plots them.###
picklefolder = '/data/user/brelethford/Data/SwiftBAT70m/pickle/'

## Read in background trials previously logged ##

bckg_uniform = cache.load(picklefolder+'bckg_trials_uniform.pickle')
bckg_redshift = cache.load(picklefolder+'bckg_trials_redshift.pickle')
bckg_flux = cache.load(picklefolder+'bckg_trials_flux.pickle')

##Now we make hists of the test statistics ##
bins = 80
range = (0.0,20.0)
uniform_hist =histlite.Hist.normalize(histlite.hist(bckg_uniform['TS'],bins=bins,range=range))
redshift_hist =histlite.Hist.normalize(histlite.hist(bckg_redshift['TS'],bins=bins,range=range))
flux_hist =histlite.Hist.normalize(histlite.hist(bckg_flux['TS'],bins=bins,range=range))

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

##I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chifit_uniform = chi2.fit(bckg_uniform['TS'])[0]
chifit_redshift = chi2.fit(bckg_redshift['TS'])[0]
chifit_flux = chi2.fit(bckg_flux['TS'])[0]

chi_degs = [1,2,chifit_uniform,chifit_redshift,chifit_flux]
colors=['black','gray','blue','red','green']
for df,color in zip(chi_degs,colors):
  x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99999,df), 100)
  rv = chi2(df)
  chi_dist=rv.pdf(x)
  ax.plot(x, chi_dist/sum(chi_dist), linestyle=':',color=color, label=r'$\tilde{\chi}^2$: df='+str(round(df,2)))

histlite.plot1d(ax,uniform_hist,histtype='step',label='uniform',color='blue')
histlite.plot1d(ax,redshift_hist,histtype='step',label='redshift',color='red')
histlite.plot1d(ax,flux_hist,histtype='step',label='flux',color='green')
ax.set_title(r'Background TS for 70m Survey')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/bckg_TS_70m.pdf')

