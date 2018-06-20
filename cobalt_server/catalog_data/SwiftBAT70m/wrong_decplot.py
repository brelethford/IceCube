#!/usr/bin/env python

'''Make sure to source this before running: /opt/i3shared/meta-projects/icerec/V04-10-00/build/env-shell.sh'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import itertools
import tables
import healpy as hp
import pylab as py
from scipy import optimize as opt
from scipy.stats import chi2, norm
from scipy import stats as st
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, neutrinoflux, histlite, NewNuFlux
import re
misc.tex_mpl_rc()
propsmall = mpl.font_manager.FontProperties (size='small')
w=4
#The following script is used for plotting the SwiftBAT catalog in healpy.


filename=open('/data/user/brelethford/Data/SwiftBAT70m/SeyfertPure.csv','r')
rawdata=filename.read()
filename.close()
table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
#table[0] gives each column title. Data gives each row.
data=table[1:]
for i in range(len(data)):
  for j in range(len(data[i])):
    data[i][j]=data[i][j].split(',')

#need to reshape every three rows into one row.
catalog=[list(itertools.chain.from_iterable(data[i])) for i in range(len(data))]



points=[(float(catalog[i][2]),float(catalog[i][3])) for i in range(len(catalog))]
ra,dec=zip(*points)
###### And here is what I had been previously doing to my source catalogue before stacking sensitivities.
dec = np.sin(dec)
######
redshift=[float(catalog[i][15]) for i in range(len(catalog))]
gamma=[float(catalog[i][11]) for i in range(len(catalog))]
flux = [float(catalog[i][7]) for i in range(len(catalog))]
#Let's plot the flux and redshift distributions:
'''
figredshift = plt.figure (figsize=(w, .75*w))
ax=plt.gca()
ax.semilogy()
histplot=histlite.hist(redshift,bins=14,range=(0,.35))
histlite.plot1d(ax,histplot)
ax.set_xlabel('redshift')
ax.set_ylabel('sources')
ax.legend(loc='upper left', prop=propsmall)
ax.set_title ('redshift distribution')
plt.subplots_adjust (left=.2, bottom=.2)
figredshift.savefig('/home/relethford/Desktop/Link to Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/redshift.pdf')

figgamma = plt.figure (figsize=(w, .75*w))
ax=plt.gca()
ax.semilogy()
histplot=histlite.hist(gamma,bins=20,range=(0,4))
histlite.plot1d(ax,histplot)
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel('sources')
ax.legend(loc='upper left', prop=propsmall)
ax.set_title (r'$\gamma$ distribution')
plt.subplots_adjust (left=.2, bottom=.2)
figgamma.savefig('/home/relethford/Desktop/Link to Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/gamma.pdf')

figflux = plt.figure (figsize=(w, .75*w))
ax=plt.gca()
histplot=histlite.hist(flux,bins=15,range=(0,300))
histlite.plot1d(ax,histplot)
ax.semilogy()
ax.set_xlabel('flux')
ax.set_ylabel('sources')
ax.legend(loc='upper left', prop=propsmall)
ax.set_title ('flux distribution')
plt.subplots_adjust (left=.2, bottom=.2)
figflux.savefig('/home/relethford/Desktop/Link to Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/flux.pdf')
'''
#Now to make the flux evidenced in the healpy maps. To do this, I want to get the fluxes normalized so that the biggest one is... say... 100 (30). I'll also take the log10 of each point so that the amplitudes are logarithmically scaled.

fluxplot=[(i*100)/max(flux) for i in flux]
#fluxplot=[(np.log10(i)*30)/np.log10(max(flux)) for i in flux]

#Rot = 180 shifts the axis from the center to the edge
#centered by default
hp.mollview(title = 'Equatorial Map of SwiftBAT AGNs', cbar = False,
            rot = 180, notext = True, cmap=None, coord='C')
hp.graticule(coord='C', color='DimGrey')
py.title(r"wrongdec - Equatorial", fontsize = 25, fontweight='bold')
#hp.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
hp.projscatter(ra,dec,coord='C',lonlat=True, s = fluxplot)
hp.projtext(185, 2,'180', lonlat=True, fontweight = 'bold')
hp.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
hp.projtext(275, 2, '270', lonlat=True, fontweight = 'bold')
hp.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
hp.projtext(359, 2, '360', lonlat=True, fontweight = 'bold')
hp.projtext(193, -8, 'RA (deg)', lonlat=True, fontweight = 'bold')
hp.projtext(350, 30, '30', lonlat=True, fontweight = 'bold')
hp.projtext(340, 60, '60', lonlat=True, fontweight = 'bold')
#hp.projtext(5, -5, 'Dec (deg)', lonlat=True)
hp.projtext(358, -33.5, '-30', lonlat=True, fontweight = 'bold')
hp.projtext(358, -63.5, '-60', lonlat=True, fontweight = 'bold')
ra_gplane = np.arange(0.,361.,1.)
dec_gplane = np.zeros(len(ra_gplane))

hp.projplot(ra_gplane, dec_gplane, coord='G', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)

py.savefig('/data/user/brelethford/AGN_Core/Plots/wrongdec_eq.png')

#Galactic plane requires no rotation to put it in the form we want.
hp.mollview(title = 'Galactic Map of SwiftBAT AGNs', cbar = False,
            cmap = None, notext = 'True', coord = 'G')
hp.graticule(coord='G',color='DimGrey')
py.title(r"wrongdec - Galactic", fontsize = 25, fontweight = 'bold')
#hp.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
hp.projscatter(ra,dec, coord='C', lonlat=True, s = fluxplot)
hp.projtext(175, 2,'180', lonlat=True, fontweight = 'bold')
hp.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
hp.projtext(275, 2, '-90', lonlat=True, fontweight = 'bold')
hp.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
hp.projtext(-160, 2, '-180', lonlat=True, fontweight = 'bold')
hp.projtext(0, -8, 'Latitude (deg)', lonlat=True, fontweight = 'bold')
hp.projtext(175, 30, '30', lonlat=True, fontweight = 'bold')
hp.projtext(165, 60, '60', lonlat=True, fontweight = 'bold')
#hp.projtext(5, -5, 'Longitude (deg)', lonlat=True)
hp.projtext(178, -33.5, '-30', lonlat=True, fontweight = 'bold')
hp.projtext(175, -63.5, '-60', lonlat=True, fontweight = 'bold')
ra_cplane = np.arange(0.,361.,1.)
dec_cplane = np.zeros(len(ra_cplane))
hp.projplot(ra_cplane, dec_cplane, coord='C', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)

py.savefig('/data/user/brelethford/AGN_Core/Plots/wrongdec_gal.png')

print ( 'number of sources = ' + str(len(ra)))

