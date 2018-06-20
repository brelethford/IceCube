#!/usr/bin/env python

import numpy as np
import icecube.astro as astro
from icecube.umdtools import cache

rawparams = np.genfromtxt('blackhole233.csv',skip_header=True, delimiter=',').T

l,b,z,mass,flux2micron,d1Mpc,M_by_R2 = rawparams[1:]

ra, dec = astro.gal_to_equa(l,b)

params = {'ra':list(ra),'dec':list(dec),'z':list(z),'mass':list(mass),'flux2m':list(flux2micron),'distMpc':list(d1Mpc),'M_by_R2':list(M_by_R2)}

cache.save(params,'/data/user/brelethford/Data/blackhole/pickle/params.pickle')
