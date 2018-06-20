#!/usr/bin/env python

import numpy as np
from icecube.umdtools import cache

#make sources
N=50

dec_limits = np.radians([-85.0,85.0])
ra_limits = np.radians([0.0,360.0])
decs = np.linspace(dec_limits[0],dec_limits[1],N)
ras = np.linspace(ra_limits[0],2*ra_limits[1],N)%(2*np.pi)
ns = np.arange(1,N+1)

params = {'n':ns,'dec':decs,'ra':ras}

cache.save(params,'params.pickle')
#Also save it as pickle and a txtfile for the svn site
outfolder = "/data/i3home/brelethford/csky/stacktest/teststack50/"
outfile = open(outfolder+"params.txt","w")
outfile.write('N     Dec (radians)     RA (radians)\n')
for n,dec,ra in zip(ns,decs,ras):
  outfile.write('{0}     {1:f}     {2:f}\n'.format(n,dec,ra))



