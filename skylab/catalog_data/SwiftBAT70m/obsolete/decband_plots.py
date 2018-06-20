'''Make sure to source this before running: /opt/i3shared/meta-projects/icerec/V04-11-02/build/env-shell.sh'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import tables
import os
from icecube.umdtools import cache
from icecube import icetray, dataclasses, neutrinoflux, histlite, NewNuFlux, astro


#Okay. Let's start by getting the data from all the years onto here. We'll put the SwiftBAT catalogue in later - first we want to see the PS data, then we'll match it up with the SwiftBAT catalogue.
projfolder='/home/relethford/Documents/IceCube_Research/'
datafolder=projfolder+'Data/AGN_Core_Sample/'

filename_pickle = projfolder+'Scripts/AGN_Core/pickle/'
if not os.path.exists(filename_pickle):
	os.mkdir(filename_pickle)

#NOTE - Have to get a better way of combining these datasets (I don't want to write everything twice!)

#Some of the following data is not yet necessary, or I just don't know how to access those types of files. Right now IC86I will be sufficient.
#f40= gz files
#f79_MESE=tables.openFile (datafolder+'IC79_MESE/MergedDataIC79.hd5')
#f79_SplinePS= root or C files
#f86_1_MESE=tables.openFile (datafolder+'IC86_2011_MESE/MergedMESEDataIC86.hd5')

#For the datasets with up and down separated, I'm going to use a dictionary to load them in, for ease of access when combining them through the getBoth function later.

fsim_86_1={'up':tables.openFile (datafolder+'IC86_2011_PS/MonteCarlo/PrunedNugenUpgoing4866.hd5'),'down':tables.openFile (datafolder+'IC86_2011_PS/MonteCarlo/PrunedNugenDowngoing4866.hd5')}

fdata_86_1={'up':tables.openFile (datafolder+'IC86_2011_PS/Data/PrunedDataUpgoing.hd5'),'down':tables.openFile (datafolder+'IC86_2011_PS/Data/PrunedDataDowngoing.hd5')}
#f_86_2=tables.openFile (datafolder+'IC86_2012_MESE/MergedMESEIC86_II.hd5')

#Let me define a function which takes a datafile (such as fsim_86_1) and applies a given string (such as SplineMPE.cols.zenith) to get both values (zenith for upgoing and zenith for downgoing)

def getBoth(nodeName, colName, datafile):
  file_up=datafile['up']
  file_down=datafile['down']
  value_up = file_up.getNode('/',nodeName).col(colName)[:]
  value_down = file_down.getNode('/',nodeName).col(colName)[:]
  combined = np.array(list(value_up)+list(value_down))
  return combined

#for right now, I'm just gonna try to plot the two datasets onto a skyplot, since where the events come from will become relevant later. So: find ra and dec (90-zen) for each event. Using the icecube astro package, this is easy enough.  

#We need zen and azi in radians in order to calculate dec and ra. luckily they are given in this form.
zen86I_sim = getBoth('SplineMPE','zenith',fsim_86_1)
zen86I_data = getBoth('SplineMPE','zenith',fdata_86_1)

azi86I_sim = getBoth('SplineMPE','azimuth',fsim_86_1)
azi86I_data = getBoth('SplineMPE','azimuth',fdata_86_1)

#Get the reco energy while we're at it...

muex86I_sim =getBoth('SplineMPEMuEXDifferential','energy',fsim_86_1)
muex86I_data = getBoth('SplineMPEMuEXDifferential','energy',fdata_86_1)

#Doing ra and dec calculations take a while - pickle?
mjd86I_sim=getBoth('SplineMPE','time',fsim_86_1)
mjd86I_data = getBoth('SplineMPE','time',fdata_86_1)

#I'll recalculate the ras, since we need to unblind before using the actual ras.
def get_coords():
	ra86I_sim,dec86I_sim=astro.dir_to_equa(zen86I_sim,azi86I_sim,mjd86I_sim)
	ra86I_data,dec86I_data=astro.dir_to_equa(zen86I_data,azi86I_data,mjd86I_data) 
	#scramble RA! We can't know this until unblinding. (Should I fix these in place, or renew them each time?)
	ra86I_sim = np.random.random(len(ra86I_sim))*2*np.pi
	ra86I_data = np.random.random(len(ra86I_data))*2*np.pi
	return ra86I_sim, dec86I_sim , ra86I_data, dec86I_data

ra86I_sim, dec86I_sim , ra86I_data, dec86I_data = cache.get (filename_pickle+"coords.pickle", get_coords)

#Let's now make a histogram of the sin(decs). In order to hist a MC set, we should first weight it with...

honda = NewNuFlux.makeFlux('honda2006')
honda.knee_reweighting_model = 'gaisserH3a_elbert' 
flux = honda.getFlux

#Weighting fcn for MC. Later on we'll need to limit the data to a certain dec band, so let's make the masking optional. To pass it a mask, we'll need to specify the dec band (i.e., pass the function [bandmask_sim[bandnum]]).

def simweight(fsim,gamma,mask=None):
	if mask == None:
	  mask=np.cast[np.bool](np.ones_like(getBoth('I3MCWeightDict','OneWeight',fsim)))
	OneWeight=getBoth('I3MCWeightDict','OneWeight',fsim)[:][mask]
	zen_sim = getBoth('MCPrimary1','zenith',fsim)[mask]
	energy_sim=getBoth('MCPrimary1','energy',fsim)[mask]
	type_sim=getBoth('MCPrimary1','type',fsim)[mask]
	type_sim[type_sim==68]=14 #neutrinos
	type_sim[type_sim==69]=-14 #antineutrinos
	n_gen = np.unique (getBoth('I3MCWeightDict','NEvents',fsim)[mask]) * len (np.unique (getBoth('I3EventHeader','Run',fsim)[mask]))
	fluence = flux (type_sim, energy_sim, np.cos(zen_sim))
	factor=2
	atmosweight = (OneWeight / n_gen) * factor * fluence
	astroweight = (OneWeight / n_gen)*(1/1e18)*(np.power((energy_sim/1e5),-1*gamma))
	return atmosweight #+ astroweight


#We'll need each MC set's livetime in order to weight properly.... not sure why?
liveTime86I=332.61*86400

#dec Hist (using gamma=2)
#The following just shows the populations per dec band of MC and Data.

fig_dec_compare=plt.figure()
ax=plt.gca()
histdec_data=histlite.plot1d(ax,histlite.hist(np.sin(dec86I_data),bins=30), color='r', label='data')
histdec_sim=histlite.plot1d(ax,histlite.hist(np.sin(dec86I_sim), bins=30, weights=liveTime86I*simweight(fsim_86_1,2)),color='b', label='MC')
ax.set_xlabel('sin(dec)')
ax.set_ylabel('events')
ax.legend(loc='upper right')
ax.set_ylim(ymin=0)
ax.set_title('event population per dec band in IC86-I - data')
fig_dec_compare.savefig(projfolder+'Plots/AGNCore/Stacking/IC86I_populations_by_dec')

#Now that I have the dec distribution, I want to divide into dec bands (let's say 6 bands of 30 deg each), and plot the energy (MuEx) of both data and MC for a given spectral parameter (gamma = 1,2,2.5,3)

#The following two lines are masks corresponding to dec bands - the first bandmask is a mask to retrieve data which lies in the lowest decband, the second bandmask the second lowest dec band, etc.
def decdiv(divisions):
	return np.arange(-np.pi/2.,np.pi/2.,np.pi/float(divisions))

decbands = decdiv(7)
bandmask_sim=[np.ma.getmask(np.ma.masked_inside(dec86I_sim, decbands[i], decbands[i+1])) for i in range(len(decbands)-1)]
bandmask_data=[np.ma.getmask(np.ma.masked_inside(dec86I_data, decbands[i], decbands[i+1])) for i in range(len(decbands)-1)]


hargs1=dict(bins=(24),range=(10,1e7),log=(True))

#Let's define a function to calculate weighting for MC MuEx, given a gamma. For now, I'll use 1 and 10^-18 for atmosnorm and astronorm, respectively.

def energyhist_sim(fsim,gamma,bandnum):
	weight=simweight(fsim,gamma,mask=bandmask_sim[bandnum])
	return histlite.hist (muex86I_sim[bandmask_sim[bandnum]],weights=liveTime86I*weight,**hargs1)

#Let's plot the lowest and highest bands. Later, we'll make a function to plot any dec band. With 7 divisions, each dec band will be a 30 degree slice:

decrange=[r'-90$^{\circ}$ to -60$^{\circ}$',r'-60$^{\circ}$ to -30$^{\circ}$',r'-30$^{\circ}$ to 0$^{\circ}$',r'0$^{\circ}$ to 30$^{\circ}$',r'30$^{\circ}$ to 60$^{\circ}$',r'60$^{\circ}$ to 90$^{\circ}$']

def plotDecBand(bandnum):
  fig_band=plt.figure()
  ax=plt.gca()
  histdec_data=histlite.plot1d(ax,histlite.hist(muex86I_data[bandmask_data[bandnum]],**hargs1), color='r', label='data')
#use gamma=2,bandnum=0
  histdec_sim=histlite.plot1d(ax,energyhist_sim(fsim_86_1,2,bandnum), color='b', label='MC')
  ax.loglog()
  ax.set_xlabel('MuEx')
  ax.set_ylabel('events')
  ax.legend(loc='upper right')
  ax.set_ylim(ymin=0)
  ax.set_title('energy pdf in IC86-I, gamma=2: '+decrange[bandnum])
  fig_band.savefig(projfolder+'Plots/AGNCore/Stacking/IC86I_energydist: band '+str(bandnum))

for i in range(6):
  plotDecBand(i)
