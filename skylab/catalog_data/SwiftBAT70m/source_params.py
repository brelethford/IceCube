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

###In this script I'll show some characteristics of the catalogue.###
picklefolder = '/data/user/brelethford/Data/SwiftBAT70m/pickle/'

params=cache.load(picklefolder+'params.pickle')
src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

##We'll plot redshift, flux, and 1/distance squared. since removing the symb sources, there are none with z=0, so we can do this easily now.##

#fig_redshift= plt.figure (figsize=(w, .75*w))
#ax=plt.gca()
#redshift_hist =histlite.hist(redshift,bins=100, log=True)
#histlite.plot1d(ax,redshift_hist,histtype='step')
#ax.set_title(r'Redshift Distribution for 70m Survey')
#ax.set_xlabel(r'z')
#plt.subplots_adjust (left=.2, bottom=.2)
#ax.set_ylabel(r'Number of Sources')
#ax.loglog()
#plt.legend(loc='upper right', prop=propsmall)
#fig_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/70m_redshift.pdf')
#
#fig_flux= plt.figure (figsize=(w, .75*w))
#ax=plt.gca()
#flux_hist =histlite.hist(flux,bins=100,log=True)
#histlite.plot1d(ax,flux_hist,histtype='step')
#ax.set_title(r'Flux Distribution for 70m Survey')
#ax.set_xlabel(r'flux $[10^{-12} \textrm{ergs}/\textrm{s}/\textrm{cm}^2]$')
#plt.subplots_adjust (left=.2, bottom=.2)
#ax.set_ylabel(r'Number of Sources')
#ax.loglog()
#ax.set_ylim(0.8)
#plt.legend(loc='upper right', prop=propsmall)
#fig_flux.savefig('/data/user/brelethford/AGN_Core/Plots/70m_flux.pdf')
#fig_flux.savefig('/data/user/brelethford/AGN_Core/Plots/70m_flux.png')
#
#fig_dist= plt.figure (figsize=(w, .75*w))
#ax=plt.gca()
#dist_hist =histlite.hist((np.power(redshift,-2.)),bins=100, log=True)
#histlite.plot1d(ax,dist_hist,histtype='step')
#ax.set_title(r'$\frac{1}{r^2}$ Distribution for 70m Survey')
#ax.set_xlabel(r'$\frac{1}{r^2}$ in arbitrary units')
#plt.subplots_adjust (left=.2, bottom=.2)
#ax.set_ylabel(r'Number of Sources')
#ax.loglog()
#ax.set_ylim(0.8)
#plt.legend(loc='upper right', prop=propsmall)
#fig_dist.savefig('/data/user/brelethford/AGN_Core/Plots/70m_dist.pdf')
#fig_dist.savefig('/data/user/brelethford/AGN_Core/Plots/70m_dist.png')
#
#fig_redshiftxflux= plt.figure (figsize=(w, .75*w))
#ax=plt.gca()
#plt.scatter(flux,redshift,s=1)
#ax.set_title(r'Redshift vs. Flux')
#ax.set_xlabel(r'flux $[10^{-12} ergs/s/cm^2]$')
#plt.subplots_adjust (left=.2, bottom=.2)
#ax.set_ylabel(r'z')
#plt.legend(loc='upper right', prop=propsmall)
#fig_redshiftxflux.savefig('/data/user/brelethford/AGN_Core/Plots/flux_vs_z.pdf')
#
#print ('Number of sources = ' + str(len(src_ra)))

#It's also probably worth having a sindec hist of the sources, but I want to weight the sources to each of the three schemes. First let's get in sindec...
sindec = np.sin(src_dec)

#Now I have to establish three hists - one for each weighting scheme...
bins = 100
range = (-1.0,1.0)
unif_hist =histlite.Hist.normalize(histlite.hist(sindec, bins=bins, range=range),integrate=True)
flux_hist =histlite.Hist.normalize(histlite.hist(sindec, weights = np.array(flux), bins=bins, range=range),integrate=True)
redshift_hist =histlite.Hist.normalize(histlite.hist(sindec, weights = np.power(redshift,-2), bins=bins, range=range),integrate=True)
 
#Now plot.

fig_sindec= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
histlite.plot1d(ax,unif_hist,histtype='step', color = 'k', label = 'equal')
histlite.plot1d(ax,flux_hist,histtype='step', color = 'green', label = 'flux')
histlite.plot1d(ax,redshift_hist,histtype='step', color = 'red', label = 'redshift')
ax.set_title(r'Weighted Source Declinations')
ax.set_xlabel(r'$\sin{(\delta)}$')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Probability Density')
plt.legend(loc='upper left', prop=propsmall)
fig_sindec.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/SwiftBAT70m/weighted_decs.pdf')
fig_sindec.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/SwiftBAT70m/weighted_decs.png')

#Now for gamma - I want to see the spectral index for the catalog.
bins = 30
range = (1,4.0)
gamma_hist = histlite.hist(gamma, bins=bins, range=range)
 
#Now plot.

fig_gamma= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
histlite.plot1d(ax,gamma_hist,histtype='step', color = 'k')
ax.set_title(r'X-ray Spectral index')
ax.set_xlabel(r'${\gamma}$')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Counts')
plt.legend(loc='upper left', prop=propsmall)
fig_gamma.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/SwiftBAT70m/gamma.pdf')

#scatterplots

fig_sindec= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
plt.scatter(np.sin(src_dec), np.ones_like(src_dec)/len(src_dec), color = 'k', s=1, label = 'equal')
plt.scatter(np.sin(src_dec), np.array(flux) / sum(flux), color = 'green', s=1, label = 'flux')
redshift_flux = np.power(redshift,-2)
plt.scatter(np.sin(src_dec), redshift_flux/sum(redshift_flux), color = 'red', s = 1, label = 'redshift')
ax.set_title(r'Source Weights')
ax.set_xlabel(r'$\sin{(\delta)}$')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Weight')
ax.set_xlim(-1.,1.)
plt.legend(loc='upper left', prop=propsmall)
fig_sindec.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/SwiftBAT70m/scatter_weights.pdf')
fig_sindec.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/SwiftBAT70m/scatter_weights.png')
