#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from skylab.psLLH_stack import PointSourceLLH #What I'm using - uses energy info
from skylab.ps_model_stack import ClassicLLH  #calculates llh from spatial information only
from scipy.stats import chi2
import healpy as hp
import itertools
from scipy.signal import convolve2d
from skylab.ps_injector_stack import PointSourceInjector
from skylab.psLLH_stack import MultiPointSourceLLH
from skylab.utils import poisson_weight
from optparse import OptionParser
import argparse

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'

## Input my commandline parameters ##

parser = OptionParser ( usage = '%prog [options]')
parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'GAMMA',
                help = 'spectral index for injection')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'uniform', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model for point source searches.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = 'uniform', metavar = 'INJWEIGHT',
                help = 'Sets the weighting used for signal injection for point source searches.')

parser.add_option ('--longrun', dest = 'longrun', 
                default = False, action = 'store_true',
                help = 'Allows longer runtime for sensitivity calculation (>40 hours)')

parser.add_option ('--n', dest = 'n', type = int, 
                default = 4, metavar = 'N',
                help = 'Number of years of data')


opts,args = parser.parse_args ()
Gamma=opts.gamma

llhweight = opts.llhweight
injweight = opts.injweight
n = opts.n

def get_SwiftBAT_params():
    filename=open(datafolder + 'SwiftBAT70m/SeyfertPure.csv','r')
    rawdata=filename.read()
    filename.close()
    table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
    data=table[1:]
    for i in range(len(data)):
      for j in range(len(data[i])):
        data[i][j]=data[i][j].split(',')
    catalog=[list(itertools.chain.from_iterable(data[i])) for i in range(len(data))]
    points=[(float(catalog[i][2]),float(catalog[i][3])) for i in range(len(catalog))]
    ra,dec=zip(*points)
    ra=list(ra)
    dec=list(dec)
    redshift=[float(catalog[i][15]) for i in range(len(catalog))]
    lum =[float(catalog[i][16]) for i in range(len(catalog))]
    gamma=[float(catalog[i][11]) for i in range(len(catalog))]
    flux = [float(catalog[i][7]) for i in range(len(catalog))]
    params={'ra':ra,'dec':dec,'redshift':redshift,'gamma':gamma,'flux':flux,'lum':lum}
    return params

params = cache.get (datafolder + 'SwiftBAT70m/pickle/params.pickle', get_SwiftBAT_params)


###Now I'll read in the background TS stuff.###
datafolder='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}yr/{1}/background_trials/'.format(str(n),llhweight)

files = [cache.load(datafolder+file) for file in os.listdir(datafolder) if file.endswith('.array')]
n_inj=[]
nsources=[]
TS=[]
beta=(0.5) #For background ts
TS_beta=[] #Calculated from the total TS median after we get all the TS.
beta_err=[]
gamma=[]

for file in files:
  for item in range(len(file['n_inj'])):
    n_inj.append(file['n_inj'][item])
    nsources.append(file['nsources'][item])
    TS.append(file['TS'][item])
    gamma.append(file['gamma'][item])

TSs=TS
TS_beta = np.percentile(TSs, 100.*(1. - beta))
m=np.count_nonzero(np.asarray(TSs) > (TS_beta))
i = len(TSs)
fraction = float(m)/float(i)
beta_err = (np.sqrt(fraction * (1. - fraction) / float(i)) if 0 < beta < 1 else 1.)

##Now we have all the pieces of the original dictionary. Time to glue bckg_trials back in place, in their proper file type.##
bckg_trials = {'n_inj':n_inj,'nsources':np.asarray(nsources), 'TS':np.asarray(TS), 'beta':beta, 'beta_err':beta_err, 'TS_beta':TS_beta, 'gamma':np.asarray(gamma)}

src_ra, src_dec, redshift, gamma, flux, lum =  params['ra'], params['dec'], params['redshift'], params['gamma'], params['flux'], params['lum']
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)

## Time to define my modelweights in a dictionary. ##


## Now to import my llh model framework. ##
import data_multi

## Like in the background trials, we have to define which llhmodel to use.
llh79 = data_multi.init79(energy=True,mode = 'box')
llh86I= data_multi.init86I(energy=True,mode = 'box')
llh59= data_multi.init59(energy=True,mode = 'box')
llh40= data_multi.init40(energy=True,mode = 'box')


#We've loaded in the appropriate llh samples, now let's put them both in the blender (not sure about weighting)
samples = [llh40,llh59,llh79,llh86I]

llhmodel = data_multi.multi_init(samples,energy=True)

mc = llhmodel.mc

total_mu = 0
for i in range(len(src_dec)):
    inj = PointSourceInjector(gamma = 2.0, sinDec_bandwidth = 0.05, src_dec = src_dec[i], seed=0)
    inj.fill(src_ra[i],src_dec[i], mc)
    total_mu += inj.flux2mu(5.61e-12/519.)

print(total_mu)
