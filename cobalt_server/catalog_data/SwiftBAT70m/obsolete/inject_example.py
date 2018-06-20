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
import re

#The following script is used to show a random skyplot, then one with injected signal.

points=[(np.random.uniform(0,360),np.random.uniform(-90,90)) for i in range(300)]
ra,dec=zip(*points)

points_sig=[(np.random.uniform(130,150),np.random.uniform(20,40)) for i in range(20)]
ra_sig,dec_sig=zip(*points_sig)

#Now to make the flux evidenced in the healpy maps. To do this, I want to get the fluxes normalized so that the biggest one is... say... 100 (30). I'll also take the log10 of each point so that the amplitudes are logarithmically scaled.

#Rot = 180 shifts the axis from the center to the edge
#centered by default
hp.mollview(title = 'Equatorial Map of SwiftBAT AGNs', cbar = False,
            rot = 180, notext = True, cmap=None, coord='C')
hp.graticule(coord='C', color='DimGrey')
py.title("Isotropic Distribution")
#hp.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
hp.projscatter(ra,dec,coord='C',lonlat=True, color='b')
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

py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/background_example.png')


hp.mollview(title = 'Equatorial Map of SwiftBAT AGNs', cbar = False,
            rot = 180, notext = True, cmap=None, coord='C')
hp.graticule(coord='C', color='DimGrey')
py.title("Isotropic Distribution with Injected Signal")
#hp.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
hp.projscatter(ra,dec,coord='C',lonlat=True, color='b')
hp.projscatter(ra_sig,dec_sig,coord='C',lonlat=True, color = 'r')
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

py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/signal_example.png')
