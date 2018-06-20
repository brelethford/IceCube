#!/usr/bin/env python

'''Make sure to source this before running: /opt/i3shared/meta-projects/icerec/V04-10-00/build/env-shell.sh'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
from scipy import optimize as opt
from scipy.stats import chi2, norm
from scipy import stats as st
import tables
from icecube.umdtools import cache
from icecube import icetray, dataclasses, neutrinoflux, histlite, NewNuFlux

#The following script is used for the purposes of categorizing the Swift BAT 70M catalogue.

filename=open('/home/relethford/Documents/Papers/X-Ray Catalogs/BAT70m_basics','r')
rawdata=filename.read()
table = [map(str, row.split()) for row in rawdata.strip().split("\n")]
data=table[1:]
datatype=[]
for i in range(len(data)):
	try:
		datatype.append(data[i][2].split(',')[-1])
	except:
		pass
'''
#This bit gives us the numbers associated with an assoc. of 0. To get a catalogue with ONLY assoc. = 0s, we should do a mask on the table/rawdata before finding the seyfert galaxies.
data_assoc=[]
for i in range(len(data)):
	try:
		data_assoc.append(data[i][2].split(',')[-3])
	except:
		pass
'''

datatype=filter(lambda a: a!= 'SWIFT', datatype)
uniques = np.unique(datatype)
typecount=[datatype.count(string) for string in uniques]
typedict_all=dict(zip(uniques,typecount))
#The next line is added to get rid of any objects which appear a substantial number of times.
typedict=dict((k,v) for (k,v) in typedict_all.iteritems() if v>10)
typedict['Unknown']=typedict.pop('')
#Above is the data of each source, but I'm sure some of these 'Types' aren't types at all. We'll figure it out soon... In any case, now I need to create a list for each type, with a number corresponding to the number of sources that are that type. I'll do this by looping through the sources...

paramhargs = {'ha': 'left',
			 'va': 'center',
			 'bbox': dict(facecolor='b', alpha=.1)}

figtype = plt.figure (1)
ypos=np.arange(len(typedict.values()))
total=sum(typedict.values())
ax=plt.gca()
ax.set_xlabel('Source Type')
ax.set_ylabel('Frequency in Catalogue')
ax.set_title('Swift BAT 70-Month Catalogue - Most Frequent Sources')
plt.xticks(ypos,typedict.keys(), rotation='vertical')
plt.bar(ypos, typedict.values(), align="center")
ax.text (4, 140, 'Total # of sources = ' + str(total), **paramhargs)

plt.subplots_adjust(bottom=0.2)
figtype.savefig("/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/sourcetype.pdf")

figseyfert = plt.figure (2)

typedict_sy=dict((k,v) for (k,v) in typedict_all.iteritems() if 'Sy' in k)
ypos=np.arange(len(typedict_sy.values()))
total=sum(typedict_sy.values())
ax=plt.gca()
ax.set_xlabel('Seyfert Type')
ax.set_ylabel('Frequency in Catalogue')
ax.set_title('Swift BAT 70-Month Catalogue - Seyfert Sources')
plt.xticks(ypos,typedict_sy.keys(), rotation='vertical')
plt.bar(ypos, typedict_sy.values(), align="center")
ax.text (4, 140, 'Total # of sources = ' + str(total), **paramhargs)
plt.subplots_adjust(bottom=0.2)
figseyfert.savefig("/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/seyfert_sources.pdf")

figbigseyfert = plt.figure (3)


typedict_sy=dict((k,v) for (k,v) in typedict_all.iteritems() if 'Sy' in k)
typedict_big_sy=dict((k,v) for (k,v) in typedict_sy.iteritems() if v>10)
total=sum(typedict_big_sy.values())
ypos=np.arange(len(typedict_big_sy.values()))
ax=plt.gca()
ax.set_xlabel('Seyfert Type')
ax.set_ylabel('Frequency in Catalogue')
ax.set_title('Swift BAT 70-Month Catalogue - Numerous Seyfert Sources')
plt.xticks(ypos,typedict_big_sy.keys(), rotation='vertical')
plt.bar(ypos, typedict_big_sy.values(), align="center")
ax.text (1, 140, 'Total # of sources = ' + str(total), **paramhargs)
plt.subplots_adjust(bottom=0.2)
figbigseyfert.savefig("/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/big_seyfert_sources.pdf")


