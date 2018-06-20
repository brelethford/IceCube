#!/usr/bin/env python

'''Make sure to source this before running: /opt/i3shared/meta-projects/icerec/V04-10-00/build/env-shell.sh'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import itertools
import tables
import healpy as hp
import pylab as pl
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
ra,dec=zip(*coord)

#Okay, I want to display these points in both equatorial and galactic coordinates.
print len(catalog)

def plotEq():
	hp.mollview(title="SwiftBAT 70 month Seyfert-0 Catalogue")
	hp.graticule(coord='C',color='DimGrey')
	hp.projscatter(ra,dec,coord='C')
	horizon = 90.
	hp.projtext(np.pi/2, 0.15, s='0$^\circ$', color='DimGrey')
	hp.projtext(np.pi/2, 2*np.pi-0.05, color='DimGrey',s='360$^\circ$')
	ras = np.arange(0.,361.,1.)*np.pi/180.
	decls_1 = np.ones(len(ras))*(180.-horizon) *np.pi/180.
	hp.projplot(decls_1,ras, color='DimGrey', linewidth=1., alpha=0.5, coord='G')
	hp.projscatter(0.0, 0.0, color='DimGrey', marker='s',coord='G', lonlat=True)
	
pl.show()

