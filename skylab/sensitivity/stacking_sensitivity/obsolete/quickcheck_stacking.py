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
# For this script I'm only interested in checking out the stacked sensitivities for gamma = 2, for each weighting scheme, as well as n_inj. No plotting from this script.
###########################################


#Okay, so we need to make sure that our fluxes are at 100 TeV reference energy. Currently they're at 1 TeV, so we'll have to convert them based on gamma.

def refChange(flux,gamma):
    return np.multiply(flux,100**(2-gamma))

def plotSensitivity(datafolder): 
    sensefolder = datafolder+'sensitivity/'
    stack_files = [cache.load(sensefolder+file) for file in os.listdir(misc.ensure_dir(sensefolder)) if file.startswith('gamma')]    
    ##Let's find the gammas used first.
    gammas=[float(filename[5:-6]) for filename in  os.listdir(misc.ensure_dir(sensefolder)) if filename.startswith('gamma')]
    flux_stack = [stack_file['flux'] for stack_file in stack_files]
    mu_stack = [stack_file['mu'] for stack_file in stack_files]
    results_stack = sorted(zip(gammas,mu_stack,flux_stack))
    gammas,mu,sensitivity = zip(*results_stack)
    sens100 = [refChange(sensitivity[i],gammas[i]) for i in range(len(sensitivity))]
    return gammas,mu,sens100


#I'll plot both the sens and disc now. I'll need to do this for each of two injection models.
gamma_sens_equal, mu_sens_equal, sensitivity_equal = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/uniform/')

gamma_sens_flux, mu_sens_flux, sensitivity_flux = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/flux/')

gamma_sens_redshift, mu_sens_redshift, sensitivity_redshift = plotSensitivity( '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/3yr/redshift/')

def printSens(gammalist,sensitivitylist,mulist):
  sens=0
  for i in range(len(gammalist)):
    if gammalist[i] == 2.0:
      sens = sensitivitylist[i]
      mu = mulist[i]
  return 'sens: '+str(sens)+', mu: '+  str(mu)

print(' ')

print ('equal: =' + printSens(gamma_sens_equal,sensitivity_equal,mu_sens_equal))
print ('flux: =' + printSens(gamma_sens_flux,sensitivity_flux,mu_sens_flux))
print ('redshift: =' + printSens(gamma_sens_redshift,sensitivity_redshift,mu_sens_redshift))










