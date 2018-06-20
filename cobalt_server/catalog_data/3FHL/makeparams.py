import numpy as np
from icecube.umdtools import misc, cache

cat = np.load('3FHL_All.npy')

ra, dec, flux, gamma, src_type = cat['ra'], cat['dec'], cat['flux50_2000'], cat['gamma'], cat['type']

#Finally let's save these values.

params = {'ra':list(ra),'dec':list(dec),'flux':list(flux), 'gamma':list(gamma), 'src_type':list(src_type)}

cache.save(params, '/data/user/brelethford/Data/3FHL_All/pickle/params.pickle')



