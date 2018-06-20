#!/usr/bin/env python

import numpy as np
import icecube.astro as astro
from icecube.umdtools import cache

rawparams = np.genfromtxt('10PWN.csv',skip_header=True, delimiter=',').T

ra,dec,weight = rawparams

params = {'ra':list(np.radians(ra)),'dec':list(np.radians(dec)),'weight':list(weight)}

cache.save(params,'/data/user/brelethford/Data/10PWN/pickle/params.pickle')
