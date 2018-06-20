import numpy as np
from icecube.umdtools import misc, cache

filename =open('/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/4yr_Starburst/starburstcatalog_raw.txt','r')

rawdata = filename.read()
filename.close()
#get rid of the splitting for decimals and hexcode for negative sign, as well an extra \n in front of and behind ff's, and also exponentials, and finally a comma:
decdata = rawdata.replace('\n.\n','.').replace('\xe2\x88\x92\n','-').replace('\nff\n','ff').replace('\ne\n','e').replace('\n,\n',',')
rawtable = [map(str, row.split('||')) for row in decdata.strip().split("\n")][:]
#Okay, now we've got to reshape it the way we want...
table = np.reshape(rawtable,(128,10))

#Now we've got a table with the 127 starburst galaxies from Tessa's analysis. Let's make a dictionary so we can match this list up with another list of all the starburst galaxy information.

name,ra,dec,z,DL_Gpc,S12m,S25m,S60m,S100m,references = [],[],[],[],[],[],[],[],[],[]
for i in table[1:]:
  name.append(i[0])
  ra.append(float(i[1]))
  dec.append(float(i[2]))
  z.append(float(i[3]))
  DL_Gpc.append(float(i[4]))
  S12m.append(float(i[5]))
  S25m.append(float(i[6]))
  S60m.append(float(i[7]))
  S100m.append(float(i[8]))
  references.append(i[9])

#Finally let's save these values.

params = {'name':name,'ra':ra,'dec':dec,'z':z,'DL_Gpc':DL_Gpc,'S12m':S12m,'S25m':S25m,'S60m':S60m,'S100m':S100m,'references':references}

cache.save(params, '/data/user/brelethford/Data/starburst/pickle/params.pickle')



