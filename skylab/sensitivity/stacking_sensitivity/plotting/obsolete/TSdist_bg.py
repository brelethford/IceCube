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
from optparse import OptionParser
import argparse
fitfun = statistics.delta_chi2
weibfun=statistics.weib

##This script should be able to plot the bckg TS, along with the TS for sensdisc, given any catalog and number of years - we'll ask cat, n as args.

##arguments:

parser = OptionParser (usage = '%prog [options]')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

opts, args = parser.parse_args ()
catalog = opts.catalog

llhweight = opts.llhweight
injweight = opts.injweight
if not llhweight:
  llhweight = 'equal'
if not injweight:
  injweight = llhweight
##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

#for now:
years = 4

#Should mstack vs bstack be an option? yeah, might as well.
alldata = '/data/user/brelethford/Output/stacking_sensitivity/{0}/{1}yr/{2}/'.format(catalog,str(years),llhweight)

bckg = alldata + 'background_trials/'

bins = 40
range = (0.0,10.0)

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
  print ('Number of trials is: {}'.format(str(len(TSs))))
  eta = np.float(np.count_nonzero(TSs))/np.float(len(TSs))
  print ('Percentage of TS = 0 is: {}'.format(str(1 -eta )))
  ##Now we make hists of the test statistics ##
  flux_hist =histlite.hist(TSs,bins=bins,range=range)
  return flux_hist,chi2fit_flux, eta

hist_bckg,chi_bckg, eta = makehist(bckg) 
#extract hists for injected trials
## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

labels = ['{0}yr'.format(str(years))]
flux_hists = [hist_bckg]
chi2_fits = [chi_bckg]

colors=['k']
#First for bckg
for flux_hist, chi2_fit, color,label in zip(flux_hists,chi2_fits, colors,labels):
  x = np.linspace(0,10,100)
  histlite.plot1d(ax,flux_hist.normalize(integrate=True),histtype='step',color=color,label='{} - '.format(label) + r'$\tilde{\chi}^2$' + ': df={}'.format(str(round(chi2_fit.par[0],2)))+ r', $\eta = $' + str(np.round(eta,2)) + ')')
  ax.plot(x, chi2_fit.pdf(x), linestyle=':',color=color)#, label='{} - '.format(label) + r'$\tilde{\chi}^2$' + ': df={}'.format(str(round(chi2_fit.par[0],2))))

def ic_prelim(fig, x = 0.75, y = 0.8, **kw):
    """Marks maps and plots as preliminary"""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text(x, y, "IceCube Preliminary", **kw)


#ax.set_title(r'{0} weighted'.format(llhweight,catalog))
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ic_prelim(fig_bckg, x = 0.5, y=0.3)
ax.set_ylabel(r'Probability Density') 
ax.set_ylim(8e-3,1.2)
ax.set_xlim(0,10)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall, ncol=1)
misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/catalog/{0}/TSdist/'.format(catalog))
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/catalog/{0}/TSdist/TSdist_{1}yr_{2}weight.pdf'.format(catalog,str(years),llhweight))
misc.ensure_dir ('/home/brelethford/public_html/catalogs/{0}/plots/TSdist/'.format(catalog))
fig_bckg.savefig('/home/brelethford/public_html/catalogs/{0}/plots/TSdist/TSdist_{1}yr_{2}weight.pdf'.format(catalog,str(years),llhweight))


