#!/usr/bin/env python
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


#data - mike
miketest = {'flux':[4.11814226e-09,4.37604463e-09,4.04030662e-09,3.50666124e-09,2.92795573e-09,2.31680668e-09,1.79871132e-09,1.33476691e-09,8.20183760e-10,2.81398718e-10,1.51897087e-10,1.56410475e-10,1.71821319e-10,1.84804601e-10,2.02581141e-10,2.19322079e-10,2.37784516e-10,2.56127821e-10,2.77857740e-10,3.40234607e-10,4.27790687e-10], 'dec':np.linspace(-1.,1.,21)}

#data - me
bentest = cache.load ('/data/user/brelethford/Output/stacking_sensitivity/testing/onesource.pickle')

#Now plot.
plt.clf()
fig_mu2flux= plt.figure (figsize=(w, w))
gs = mpl.gridspec.GridSpec (2,1, height_ratios=[3,1], hspace=0.0)
ax1 = plt.subplot (gs[0])
ax2 = plt.subplot (gs[1])
ax2.yaxis.set_ticks(np.arange(-0.6,1.4,0.2))

fig_mu2flux.suptitle(r'IC40 through 86I flux yielding 1 event')
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.plot(bentest['dec'],bentest['flux'], color = 'red', label = 'ben')
ax1.plot(miketest['dec'],miketest['flux'], color = 'blue', label = 'mike')
ax2.plot(bentest['dec'], [bentest['flux'][i] / miketest['flux'][i] for i in range(len(bentest['flux']))], color = 'k')
ax2.set_xlabel(r'$\sin{(\delta})$')
ax1.set_ylabel(r'flux [$E^2 \frac{dN}{dE} \frac{1}{GeV cm^2 s}$]')
ax2.set_ylabel(r'ratio (ben/mike)')
plt.subplots_adjust (left=.2, bottom=.2)
ax1.semilogy()
ax1.set_ylim(1e-10,1e-8)
ax2.set_ylim(0.6,1.4)
ax1.legend(loc='upper right', prop=propsmall)

gs.tight_layout(fig_mu2flux, rect=[0, 0.03, 1, 0.97])
#fig_IC79compare.subplots_adjust(top=0.80)
fig_mu2flux.savefig('/data/user/brelethford/AGN_Core/Plots/mu2flux_test/onesource.pdf')

