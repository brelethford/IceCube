import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
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
datafolder = '/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/'
## Read in background trials previously logged ##
files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')] 
bckg_single=[]
for file in files:
  bckg_single.append(list(file['TS']))


##Now we make hists of the test statistics ##
bins = 80
range = (0.0,20.0)
single_hist =histlite.Hist.normalize(histlite.hist(bckg_single[0],bins=bins,range=range))

## Now to plot. ##
fig_bckg = plt.figure (figsize=(w, .75*w))
ax=plt.gca()

##I'll include a chi squared distribution w/ DOF=1 (and 2, just because). I'll also show the best fitting chi2 dist for each weighting scheme.##
chifit_single = chi2.fit(bckg_single[0])[0]

chi_degs = [1,2,chifit_single]
colors=['black','gray','blue']
for df,color in zip(chi_degs,colors):
  x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99999,df), 100)
  rv = chi2(df)
  chi_dist=rv.pdf(x)
  ax.plot(x, chi_dist/sum(chi_dist), linestyle=':',color=color, label=r'$\tilde{\chi}^2$: df='+str(round(df,2)))

histlite.plot1d(ax,single_hist,histtype='step',label='single source',color='blue')
ax.set_title(r'Background TS for 70m Survey')
ax.set_xlabel(r'TS')
plt.subplots_adjust (left=.2, bottom=.2)
ax.set_ylabel(r'Normalized Counts') 
ax.set_ylim(8e-5,1)
ax.set_xlim(0,16)
ax.semilogy() 
plt.legend(loc='upper right', prop=propxsmall)
fig_bckg.savefig('/data/user/brelethford/AGN_Core/Plots/bckg_TS_single_source.pdf')

