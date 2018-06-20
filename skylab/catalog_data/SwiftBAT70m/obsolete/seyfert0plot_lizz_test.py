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

#Okay, I want to display these points in both equatorial and galactic coordinates. I'm going to shamelessly steal Lizz' code for plotting, just to see how my galactic plane lines up in correspondance to Will's code.

ra_gtest=(5,5,-5,-5)
dec_gtest=(5,-5,5,-5)

ra_ctest=(265,265,255,255)
dec_ctest=(-35,-25,-35,-25)

#This is the actual plotting part
#If you don't specify, it's in equatorial
#Rot = 180 shifts the axis from the center to the edge
#centered by default
hp.mollview(title = 'Equatorial Map of SwiftBAT AGNs', cbar = False,
            rot = 180, notext = False, cmap=None,coord='C')
hp.graticule(coord='C', color='DimGrey')
#hp.graticule(coord='G', color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Equatorial")
#hp.projscatter(ra,dec,coord='C',lonlat=True)
hp.projscatter(ra_ctest,dec_ctest, coord='C',lonlat=True)
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
dec_gplane = np.ones(len(ra_gplane))*(90)
hp.projplot(dec_gplane,ra_gplane, rot=180, coord='G', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)
py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/equatorialmap.png')

#This is the galactic one, which required no rotation BUT
#it required a coordinate transform
#which is so much easier with healpy
#coord = ['C', 'G'] basically says take equatorial map (C) and make Galactic (G)
hp.mollview(title = 'Galactic Map of SwiftBAT AGNs', cbar = False,
            cmap = None, notext = 'True', coord = 'G')
hp.graticule(coord='G',color='DimGrey')
#hp.graticule(coord='C',color='DimGrey')
py.title("SwiftBAT 70 month Seyfert-0 Catalogue - Galactic")
#hp.projscatter(ra,dec, rot=180, coord='C', lonlat=True)
hp.projscatter(ra_ctest,dec_ctest, coord='C',lonlat=True)
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
ra_eplane = np.arange(0.,361.,1.)
dec_eplane = np.ones(len(ra_eplane))*(90.)
hp.projplot(dec_eplane,ra_eplane, rot=180, coord='C', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)
py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/galacticmap.png')


#Problem - My equatorial map is showing the wrong plane, I think?

#CURRENTLY: I'm only getting the correct galactic plane on my equatorial plot by rotating BOTH the mollview projection AND the galactic plane. Is this correct?
