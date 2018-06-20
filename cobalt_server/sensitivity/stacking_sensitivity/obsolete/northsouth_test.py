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

def plotSensitivity(filename):
    stack_file = cache.load(filename)
    return stack_file['flux'], stack_file['mu']

##All I'm interested in is seeing the number of n_inj for each sample.
datafolder = '/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/'

flux_onenorth, mu_onenorth = plotSensitivity(datafolder+'onenorth.array')
flux_twonorth, mu_twonorth = plotSensitivity(datafolder+'twonorth.array')
flux_onesouth, mu_onesouth = plotSensitivity(datafolder+'onesouth.array')
flux_twosouth, mu_twosouth = plotSensitivity(datafolder+'twosouth.array')
flux_both, mu_both = plotSensitivity(datafolder+'both.array')

##Let's load in the sources so we know where they are:
params = cache.load('/data/user/brelethford/Data/SwiftBAT70m/pickle/params.pickle')

onenorthdec = [params['dec'][2]]
twonorthdec = [params['dec'][2],params['dec'][6]]
onesouthdec = [params['dec'][0]]
twosouthdec = [params['dec'][0],params['dec'][7]]
bothdec = [params['dec'][0],params['dec'][6]]

onenorthweight = [params['flux'][2]]
twonorthweight = [params['flux'][2],params['flux'][6]]
onesouthweight = [params['flux'][0]]
twosouthweight = [params['flux'][0],params['flux'][7]]
bothweight = [params['flux'][0],params['flux'][6]]

print ('one northern source at ' + str(onenorthdec)+ ', flux = '+str(onenorthweight))
print (mu_onenorth,flux_onenorth)

print ('two northern sources at ' + str(twonorthdec)+', flux = '+str(twonorthweight))
print (mu_twonorth,flux_twonorth)

print ('one southern source at ' + str(onesouthdec)+', flux = '+str(onesouthweight))
print (mu_onesouth,flux_onesouth)

print ('two southern sources at ' + str(twosouthdec)+', flux = '+str(twosouthweight))
print (mu_twosouth,flux_twosouth)

print ('one north, one south at ' + str(bothdec)+', flux = '+str(bothweight))
print (mu_both,flux_both)


##Something's wrong - how do the background TS values differ?.
def backgroundTS(filename):
  files = [cache.load(filename+file) for file in os.listdir(filename) if file.endswith('.array')]
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
  print(i) 
  ##Now we have all the pieces of the original dictionary. Time to glue bckg_trials back in place, in their proper file type.##
  bckg_trials = {'n_inj':n_inj,'nsources':np.asarray(nsources), 'TS':np.asarray(TS), 'beta':beta, 'beta_err':beta_err, 'TS_beta':TS_beta, 'gamma':np.asarray(gamma)}
  return bckg_trials['TS_beta']

bckg_onenorth = backgroundTS('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/onenorth/background_trials/')
print (bckg_onenorth)
bckg_twonorth = backgroundTS('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/twonorth/background_trials/')
print (bckg_twonorth)
bckg_onesouth = backgroundTS('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/onesouth/background_trials/')
print (bckg_onesouth)
bckg_twosouth = backgroundTS('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/twosouth/background_trials/')
print (bckg_twosouth)
bckg_both = backgroundTS('/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/both/background_trials/')
print (bckg_both)
