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

###In this script I'll show some characteristics of the catalogue.###
picklefolder = '/data/user/brelethford/Data/SwiftBAT70m/pickle/'

params=cache.load(picklefolder+'params.pickle')
src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']

N = len(src_ra)

arr= np.empty((N, ), dtype=[("ra", np.float), ("dec", np.float),
                                 ("z", np.float), ("gamma", np.float),
                                 ('flux', np.float), ('lum', np.float)])
arr["ra"] = src_ra
arr["dec"] = src_dec
arr["gamma"] = gamma
arr["z"] = redshift
arr["flux"] = flux
arr["lum"] = lum

np.save(picklefolder+'seyferts.npy', arr)
