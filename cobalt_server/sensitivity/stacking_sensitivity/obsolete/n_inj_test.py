import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import scipy as sp
import os
import sys
from icecube.umdtools import cache, misc
from scipy import optimize as opt
from icecube import icetray, dataclasses, histlite, astro
###########################################
#This function should pick out the various sensitivities for different gammas.
###########################################

misc.tex_mpl_rc()
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')
propxxsmall = mpl.font_manager.FontProperties (size='xx-small')
w=4

#Okay, so we need to make sure that our fluxes are at 100 TeV reference energy. Currently they're at 1 TeV, so we'll have to convert them based on gamma.

def refChange(flux,gamma):
    return np.multiply(flux,100**(2-gamma))

def plotSensitivity(sensefolder):
    stack_files = [cache.load(sensefolder+file) for file in os.listdir(misc.ensure_dir(sensefolder)) if file.startswith('gamma2.0')]    
    ##Let's find the gammas used first.
    gammas=[float(filename[5:-6]) for filename in  os.listdir(misc.ensure_dir(sensefolder)) if filename.startswith('gamma')]
    flux_stack = [stack_file['flux'] for stack_file in stack_files]
    mu_stack = [stack_file['mu'] for stack_file in stack_files]
    return flux_stack,mu_stack

##All I'm interested in is seeing the number of n_inj for each sample.
flux_stack, mu_stack = plotSensitivity('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/flux/flux_inj/sensitivity/')

flux_stack_cut, mu_stack_cut = plotSensitivity('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/flux/flux_inj/sensitivity/')


flux1,mu1 = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/onesource/')
flux2,mu2 = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/twosource/')

flux1_cut,mu1_cut = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/onesource/')
flux2_cut,mu2_cut = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m_mc_cut/twosource/')

print ('stacked n_inj = ')
print (mu_stack,flux_stack)

print ('stacked cut n_inj = ')
print (mu_stack_cut,flux_stack_cut)

print ('One Source n_inj = ')
print (mu1,flux1)

print ('Two Source n_inj = ')
print (mu2,flux2)


print ('One Source Cut n_inj = ')
print (mu1_cut,flux1_cut)

print ('Two Source Cut n_inj = ')
print (mu2_cut,flux2_cut)

params = cache.load ( '/data/user/brelethford/Data/SwiftBAT70m/pickle/params.pickle' )

print ( 'First ten sources in catalog: ')
print (params['dec'][0:10])

print ('check if I can staple params together:')
print (params['dec'][0:1]+params['dec'][5:8])

