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
picklefolder = '/data/user/brelethford/Data/starburst/pickle/'

params=cache.load(picklefolder+'params.pickle')
src_ra, src_dec, redshift, Dl, flux =  params['ra'], params['dec'], params['z'], params['DL_Gpc'], params['S60m']

##We'll plot redshift,  1/distance squared.

fig_redshift= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
redshift_hist =histlite.hist(redshift,bins=40, log=True)
histlite.plot1d(ax,redshift_hist,histtype='step')
ax.set_title(r'Redshift Distribution for Starburst Catalog')
ax.set_xlabel(r'z')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Number of Sources')
ax.loglog()
plt.legend(loc='upper right', prop=propsmall)
fig_redshift.savefig('/data/user/brelethford/AGN_Core/Plots/starburst/redshift.pdf')

fig_dist= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
dist_hist =histlite.hist((np.power(redshift,-2.)),bins=40, log=True)
histlite.plot1d(ax,dist_hist,histtype='step')
ax.set_title(r'$\frac{1}{r^2}$ Distribution for Starburst Catalog')
ax.set_xlabel(r'$\frac{1}{r^2}$ in arbitrary units')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Number of Sources')
ax.loglog()
ax.set_ylim(0.8)
plt.legend(loc='upper right', prop=propsmall)
fig_dist.savefig('/data/user/brelethford/AGN_Core/Plots/starburst/dist.pdf')
#Note: we can use luminosity distance for this!

print ('Number of sources = ' + str(len(src_ra)))

#It's also probably worth having a sindec hist of the sources, but I want to weight the sources to each of the three schemes. First let's get in sindec...

#Now I have to establish three hists - one for each weighting scheme...
bins = 50
range = (-90.0,90.0)
unif_hist =histlite.hist(src_dec, bins=bins, range=range)
range = (0.5,3.5)
flux_hist =histlite.hist(np.log10(flux), bins=bins, range=range)
#redshift_hist =histlite.Hist.normalize(histlite.hist(sindec, weights = np.power(redshift,-2), bins=bins, range=range))
 
#Now plot.
'''
fig_sindec= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
histlite.plot1d(ax,unif_hist,histtype='step', color = 'k', label = 'equal')
#histlite.plot1d(ax,flux_hist,histtype='step', color = 'green', label = 'flux')
#histlite.plot1d(ax,redshift_hist,histtype='step', color = 'red', label = 'redshift')
ax.set_title(r'Source declinations')
ax.set_xlabel(r'Declination')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Dec distribution')
ax.set_xlim(-90.,90.)
ax.set_ylim(0,20)
plt.legend(loc='upper left', prop=propsmall)
fig_sindec.savefig('/data/user/brelethford/AGN_Core/Plots/starburst/dec_hist.pdf')
'''
fig_flux= plt.figure (figsize=(w, .75*w))
ax=plt.gca()
histlite.plot1d(ax,flux_hist,histtype='step')
ax.set_title(r'Flux Distribution for Starburst Catalog')
ax.set_xlabel(r'log10 (Flux$_{FIR}$/Jy)')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Number of Sources')
ax.set_ylim(0,13)
ax.set_xlim(0.5,3.5)
plt.legend(loc='upper right', prop=propsmall)
fig_flux.savefig('/data/user/brelethford/AGN_Core/Plots/starburst/flux_hist.pdf')
