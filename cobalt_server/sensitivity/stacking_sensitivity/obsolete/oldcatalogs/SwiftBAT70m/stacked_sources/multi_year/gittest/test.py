#!/usr/bin/env python
import numpy as np
import scipy as sp
import tables
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
from bstacking.psLLH import PointSourceLLH #What I'm using - uses energy info
from bstacking.ps_model import ClassicLLH  #calculates llh from spatial information only
import itertools
from bstacking.ps_injector import PointSourceInjector
from bstacking.psLLH import MultiPointSourceLLH
from bstacking.utils import poisson_weight

sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core")

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/'

filename_plots=projfolder+'Plots/AGNCore/Stacking/'

## Input my commandline parameters ##
gamma = 2.0
Gamma = 2.0

llhweight = 'uniform'
injweight = 'uniform'
n=4

params = cache.load (datafolder + 'SwiftBAT70m/pickle/params.pickle')

projfolder='/home/brelethford/Documents/IceCube_Research/'
datafolder='/data/user/brelethford/Data/AGN_Core_Sample/'
filename_pickle = datafolder+'pickle/'

livetime_IC86I = 332.61

from bstacking.psLLH import PointSourceLLH, MultiPointSourceLLH, StackingPointSourceLLH
from bstacking.ps_model import  ClassicLLH, EnergyLLH, PowerLawLLH

def init86I(energy=False, mode = 'all', **kwargs): #weighting=None, mode='all', **kwargs):
    arr_exp = cache.load(filename_pickle+"IC86I/exp.pickle")
    arr_mc = cache.load(filename_pickle+"IC86I/mc.pickle")
    #This hem... not sure about what it means, but stefan uses it in load.py.
    #Obviously it's for the hemisphere line - but why doesn't he use it for every year?
    hem = np.sin(np.radians(-5.))
    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.2, 10 + 1),
                         np.linspace(-0.2, hem, 4 + 1),
                         np.linspace(hem, 0.2, 5 + 1),
                         np.linspace(0.2, 1., 10),
                         ]))
    #check this later, BEN!
    energy_bins = [np.linspace(1., 10., 67 + 1), dec_bins]

    if energy:

        llh_model = EnergyLLH(energy_bins,
                          sinDec_bins=dec_bins)
    else:

        llh_model = ClassicLLH(#["logE"], min(50, Nexp // 50),
                                #twodim_range=[0.9 * arr_mc["logE"].min(),
                                 #             1.1 * arr_mc["logE"].max()],
                                sinDec_bins=dec_bins,
                                sinDec_range=[-1., 1.])

    llh = StackingPointSourceLLH(arr_exp, arr_mc, livetime=livetime_IC86I, llh_model=llh_model,
                         mode=mode, #hemisphere=dict(Full=[-np.inf, np.inf]),
                         #nsource=Nexp / 100.,
                         #nsource_bounds=(-Nexp / 2., Nexp / 2.)
                         #               if not energy else (0., Nexp / 2.),
                         seed=np.random.randint(2**32),
                         **kwargs)

    return llh


llh86I= init86I(energy=True)

print (llh86I)

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

## Time to define my modelweights in a dictionary. ##

modelweights = {'flux':flux, 'redshift': list(np.power(redshift,-2)), 'uniform':list(np.ones_like(redshift))}

##Remember, weighted sensitivity requires src dec in radians.#
src_dec= np.radians(src_dec)
src_ra = np.radians(src_ra)

##Now, I'll input my injection weighting scheme from the commandline. ##
if injweight == 'uniform':
  inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, seed=0)
else:
  inj = PointSourceInjector(Gamma, sinDec_bandwidth=.05, src_dec= src_dec, theo_weight = modelweights['{}'.format(injweight)], seed=0) 

sensitivity = StackingPointSourceLLH.weighted_sensitivity(llh86I,src_ra=src_ra,src_dec=src_dec,alpha=.5,beta=.9,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},longrun=True,bckg_trials=bckg_trials,eps=0.02,n_iter=250, w_theo=None)
print sensitivity

#discovery = PointSourceLLH.weighted_sensitivity(llhmodel,src_ra=src_ra,src_dec=src_dec,alpha=2.867e-7,beta=.5,inj=inj,trials={'n_inj':[],'TS':[],'nsources':[],'gamma':[]},bckg_trials=bckg_trials,eps=0.01,n_iter=250)
#print discovery

#choose an output dir, and make sure it exists
this_dir = os.path.dirname(os.path.abspath (__file__))
sens_dir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}yr/{1}/{2}_inj/sensitivity/'.format(str(n),llhweight, injweight))

# save the output
outfile_sens = sens_dir + 'gamma{}.array'.format(Gamma)

print 'Saving', outfile_sens, '...'
cache.save(sensitivity, outfile_sens)
#cache.save(discovery, outfile_disc)

