import numpy as np
from icecube.umdtools import misc, cache

filename =open('/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/catalog_data/SNR_cloud/rawparams.txt','r')

rawdata = filename.read()
filename.close()
#get rid of the splitting for decimals and hexcode for negative sign, as well an extra \n in front of and behind ff's, and also exponentials, and finally a comma:
decdata = rawdata.replace('\n.\n','.').replace('\xe2\x88\x92\n','-').replace('\nff\n','ff').replace('\ne\n','e').replace('\n,\n',',')
rawtable = [map(str, row.split('||')) for row in decdata.strip().split("\n")][:]
#Okay, now we've got to reshape it the way we want...
table = [row[0].split(' ') for row in rawtable]
#Now we've got a table with the 127 starburst galaxies from Tessa's analysis. Let's make a dictionary so we can match this list up with another list of all the starburst galaxy information.

name,ra,dec,weight, extension = [],[],[],[],[]
for i in table[1:]:
  name.append(i[0])
  ra.append(float(i[1]))
  dec.append(float(i[2]))
  weight.append(float(i[3]))
  extension.append((float(i[4])))

#Finally let's save these values.

params = {'name':name,'ra':np.radians(ra),'dec':np.radians(dec),'weight':weight,'extension':np.radians(extension)}

cache.save(params, '/data/user/brelethford/Data/SNR_cloud/pickle/params.pickle')



