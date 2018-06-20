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
from icecube.umdtools import cache
from icecube import icetray, dataclasses, neutrinoflux, histlite, NewNuFlux

#The following script is used for plotting the SwiftBAT catalog in healpy.

filename=open('/home/relethford/Documents/IceCube_Research/Scripts/AGN_Core/Swift-BAT_70m/seyfert0.csv','r')
rawdata=filename.read()
filename.close()
table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
#table[0] gives each column title. Data gives each row.
data=table[1:]
# data[i][j] = [[data[i][j].split(',') for j in range(len(data[i]))] for i in range(len(data))] doesn't work right now for some reason...
for i in range(len(data)):
  for j in range(len(data[i])):
    data[i][j]=data[i][j].split(',')

#need to reshape every three rows into one row.
catalog=[list(itertools.chain.from_iterable(data[i])) for i in range(len(data))]
for i in range(len(catalog)):
  del catalog[i][7]

points=[(float(catalog[i][2]),float(catalog[i][3])) for i in range(len(catalog))]
ra,dec=zip(*points)

#Okay, I want to display these points in both equatorial and galactic coordinates. For some reason I can't rotate the whole plot, so instead I while simply rotate every point by 180 degrees so that 0 starts on the right, not the center. To do this I simply add 180 to every RA. Will says my Galactic plane is still off, but adding 180 to each RA for the galactic plane doesn't seem to do anything.

#First equatorial...
plt.clf()
hp.graticule(coord='C',color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Equatorial")
hp.projscatter(np.add(ra,180),dec,coord='C',lonlat=True)
horizon = 90.
hp.projtext(np.pi/2, np.pi+ 0.15, s='0$^\circ$', color='DimGrey')
hp.projtext(np.pi/2, 3*np.pi-0.05, color='DimGrey',s='360$^\circ$')
ra_gplane = np.add(np.arange(0.,361.,1.),180)
dec_gplane = np.ones(len(ra_gplane))*(180.-horizon)
hp.projplot(dec_gplane,ra_gplane, color='DimGrey', lonlat=True, linewidth=1., alpha=0.5, coord='G')
#plt.show()
py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/equatorialmap.png')


#Now Galactic. I assume I don't need to rotate this one.
plt.clf()
hp.graticule(coord='G',color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Galactic")
#hp.projscatter(ra,dec,coord='G', lonlat=True)
hp.projscatter(0,30,coord='C', lonlat=True)
hp.projtext(np.pi/2-0.025, -np.pi+0.4, s='-180$^\circ$', color='DimGrey',coord='G')
hp.projtext(np.pi/2-0.025, np.pi-0.05, color='DimGrey', coord='G', s='180$^\circ$')
hp.projplot(dec_gplane,ra_gplane, color='DimGrey', linewidth=1., alpha=0.5,coord='C')
plt.show()
py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/galacticmap.png')

#Okay, I must be calculating the planes wrong. For Galactic plane on the equatorial map, I should take the horizon and plot it in 'G' coordinates, right? Vice Versa for Equatorial plane on galacitc map?

#MAYBE it's because I'm rotating everything else in the plane, but not the galactic plane itself? But that doesn't account for why the equatorial is backwards, too.

#following: without rotation.
'''
#First equatorial...
hp.graticule(coord='C',color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Equatorial")
hp.projscatter(ra,dec,coord='C',lonlat=True)
horizon = 90.
hp.projtext(np.pi/2, 0.15, s='0$^\circ$', color='DimGrey')
hp.projtext(np.pi/2, 2*np.pi-0.05, color='DimGrey',s='360$^\circ$')
ras = np.arange(0.,361.,1.)*np.pi/180.
decls_1 = np.ones(len(ras))*(180.-horizon) *np.pi/180.
hp.projplot(decls_1,ras, color='DimGrey', linewidth=1., alpha=0.5, coord='G')

py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/equatorialmap.png')
plt.clf()

#Now Galactic.
hp.graticule(coord='G',color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Galactic")
hp.projscatter(ra,dec,coord='G', lonlat=True)
hp.projtext(np.pi/2-0.025, -np.pi+0.4, s='-180$^\circ$', color='DimGrey',coord='G')
hp.projtext(np.pi/2-0.025, np.pi-0.05, color='DimGrey', coord='G', s='180$^\circ$')
hp.projplot(decls_1,ras, color='DimGrey', linewidth=1., alpha=0.5,coord='C',)
py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/galacticmap.png')
plt.clf()
'''
