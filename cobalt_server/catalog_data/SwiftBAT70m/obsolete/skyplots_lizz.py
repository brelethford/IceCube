import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import itertools
import tables
import healpy as hp
import pylab as py
from scipy import optimize as opt
from scipy.stats import chi2, norm
from scipy import stats as st
from icecube.umdtools import cache
from icecube import icetray, dataclasses, neutrinoflux, histlite, NewNuFlux
#This defines big squares in the sky, and a set of offsets to move them by 
#so you have box a, and a copy of it at the values listed

region_A = dict (decs = np.array ([-0.298451302, 0.0523598776, 0.34906585]), 
                 ras = np.array ([[0.802851456, 1.32645023], [1.02101761, 1.31772359]]),
                 offsets = np.array ([-0.6, -2.2, -3.4, -1.6, 0.6, 1.8]))

region_B = dict (decs = np.array ([-0.0523598776, 0.261799388, 0.680678408, 0.959931089]), 
                 ras = np.array ([[2.07694181, 2.49582083], [1.9809487, 2.35619449], [2.0943951, 2.38237443]]),
                 offsets = np.array ([-0.6, -1.8, -2.6, -3.2, 0.6]))   

region_C = dict(decs = np.array ([0.4188790204786391, 0.5934119456780721, 0.5934119456780721]), 
                ras = np.array([[3.490658503988659, 3.7699111843077517],[0.,0.]]))

region_3 = dict(decs = np.array([0.2617993877991494, 0.5759586531581288, 0.7155849933176751, 0.9599310885968813]),
                ras = np.array([[4.310963252425994, 4.9218284906240095], [4.310963252425994, 4.590215932745087], [4.084070449666731, 4.4505895925855405]]))

k_2 = dict(decs = np.array([1.204277, 1.25663]), 
    ras = np.array([1.30899, 1.39626]) )
        
k_5 = dict(decs = np.array([1.23045, 1.2479104]),
     ras = np.array([1.34390, 1.37881]) )   


#These are kinda irrelevant, but may be useful for you at some point
#They determine if an event falls in the squares above

def is_in_region (region, data_dec, data_ra, offset=0):            
    decs =region['decs']                                                 
    ras = region['ras']
    out = np.zeros (len (data_dec), dtype=bool)
    for i in xrange (len (ras)):
        d1 = decs[i]
        d2 = decs[i+1]
        r1, r2 = (ras[i] + offset) % (2*np.pi)
        out += ((d1 <= data_dec) * (data_dec < d2) * (r1 <= data_ra) * (data_ra < r2))
    return out

def in_region (region, data_dec, data_ra, region_offset=0):            
    decs =region['decs']                                                 
    ras = region['ras']
    out = np.zeros (len (data_dec), dtype=bool)
    d1, d2 = decs
    r1, r2 = (ras + region_offset) % (2*np.pi)
    out += ((d1 <= data_dec) & (data_dec < d2) & (r1 <= data_ra) & (data_ra < r2))
    return out

#one was a dec band, so checked that

def is_in_dec(region, data_dec, data_ra):
    ras = region['ras']
    decs = region['decs']
    out = np.zeros( len(data_dec), dtype = bool)
    d1, d2 = decs
    out += ( (d1 <= data_dec) & (data_dec < d2) & ~in_region(region, data_dec, data_ra))
    return out

#This determines the resolution of your map, in terms of n (number of pixels)
n=hp.nside2npix(2**8)

#made an empty space the size of the resolution
binspace = np.zeros(n)

#and made one for coordinates
bin_dec = np.zeros(n)
bin_ra = np.zeros(n)

#so this makes a set of angles (phi, theta) for the pixels
#you literally have to calculate them...
#Then I appended the angles to the dec and ra arrays, where
#dec just needed to be converted from phi

for i in range(len(bin_dec)):
    [theta, phi] = hp.pix2ang(2**8, i)
    bin_dec[i] = np.pi/2. -theta
    bin_ra[i] = phi
    
#just checking in here:
print bin_dec.min(), bin_dec.max()
print bin_ra.min(), bin_ra.max()

#I am now coloring my squares, by asking it to see which pixels fall in which boxes
idx_A_off =np.sum([is_in_region(region_A, bin_dec, bin_ra, offset) for offset in region_A['offsets']], axis = 0) > 0
idx_B_off =np.sum([is_in_region(region_B, bin_dec, bin_ra, offset) for offset in region_B['offsets']], axis = 0) > 0

idx_A = is_in_region(region_A, bin_dec, bin_ra)
idx_B = is_in_region(region_B, bin_dec, bin_ra)
idx_C = is_in_region(region_C, bin_dec, bin_ra)
idx_3 = is_in_region(region_3, bin_dec, bin_ra)

idx_2 = in_region(k_2, bin_dec, bin_ra)
idx_5 = in_region(k_5, bin_dec, bin_ra)
idx_2_off = is_in_dec(k_2, bin_dec, bin_ra)
idx_5_off = is_in_dec(k_5, bin_dec, bin_ra)

#To get red, black and blue boxes I had to pick a convoluted color scheme
#And you then assign a value to the bin that gives you the color you want from the colormap
#yes this part fucking sucks
for i in range(len(binspace)):
    if idx_A[i] == 1:
        binspace[i] = -100
    if idx_B[i] == 1:
        binspace[i] = -100
    if idx_3[i] == 1:
        binspace[i] = 100
    if idx_C[i]  == 1:
        binspace[i] = 100
    if idx_5[i] == True:
        binspace[i] = -100
    #if idx_5[i] == 1:
    #    binspac[i] = 100
    #if idx_A_off[i] == 1:
    #    binspace[i] = -80
    #if idx_B_off[i] == 1:
    #    binspace[i] = -80
    if idx_5_off[i] == 1:
            binspace[i] = -80.5
    if binspace[i] == 0:
        binspace[i] = 80.5
        

#from matplotlib import cm

#This is said color map from above
from matplotlib import cm
flagmap = cm.flag
flagmap.set_under('w')

#This is the actual plotting part
#If you don't specify, it's in equatorial
#Rot = 180 shifts the axis from the center to the edge
#centered by default
hp.mollview(binspace, title = 'Equatorial Map of Regions', cbar = False,
            rot = 180, notext = False, cmap = flagmap)
hp.graticule()

#If you're smart, I'm sure you can write a loop for this
#you have to label the damn things by hand
hp.projtext(185, 2,'180', lonlat=True, fontweight = 'bold')
hp.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
hp.projtext(275, 2, '270', lonlat=True, fontweight = 'bold')
hp.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
hp.projtext(359, 2, '360', lonlat=True, fontweight = 'bold')
hp.projtext(193, -8, 'RA (deg)', lonlat=True, fontweight = 'bold')

hp.projtext(65, -8, 'A', lonlat = True, color = 'w', fontweight = 'bold')
hp.projtext(213, 28, 'C/4', lonlat = True, color= 'w', fontweight = 'bold')
hp.projtext(125, 25, 'B', lonlat = True, color = 'w', fontweight = 'bold')
hp.projtext(265, 25, '3', lonlat = True, color = 'w', fontweight = 'bold')
hp.projtext(95, 60, 'K_5', lonlat = True, color = 'k', fontweight = 'bold')
hp.projtext(180, 75, 'K_5 off', lonlat = True, color = 'k', fontweight = 'bold')


hp.projtext(350, 30, '30', lonlat=True, fontweight = 'bold')
hp.projtext(340, 60, '60', lonlat=True, fontweight = 'bold')
#hp.projtext(5, -5, 'Dec (deg)', lonlat=True)
hp.projtext(358, -33.5, '-30', lonlat=True, fontweight = 'bold')
hp.projtext(358, -63.5, '-60', lonlat=True, fontweight = 'bold')

plt.show()
exit()

#This is the galactic one, which required no rotation BUT
#it required a coordinate transform
#which is so much easier with healpy
#coord = ['C', 'G'] basically says take equatorial map (C) and make Galactic (G)
hp.mollview(binspace, title = 'Galactic Map of Regions', cbar = False,
            cmap = flagmap, notext = 'True', coord = ['C', 'G'])
hp.graticule()

hp.projtext(175, 2,'180', lonlat=True, fontweight = 'bold')
hp.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
hp.projtext(275, 2, '-90', lonlat=True, fontweight = 'bold')
hp.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
hp.projtext(-160, 2, '-180', lonlat=True, fontweight = 'bold')
hp.projtext(0, -8, 'Latitude (deg)', lonlat=True, fontweight = 'bold')

hp.projtext(55, 35, '3', lonlat = True, color = 'w', fontweight = 'bold')
hp.projtext(207, 30, 'B', lonlat = True, color= 'w', fontweight = 'bold')
hp.projtext(65, 60, 'C/4', lonlat = True, color = 'k', fontweight = 'bold')
hp.projtext(200, -40, 'A', lonlat = True, color = 'w', fontweight = 'bold')
hp.projtext(135, 15, 'K_5', lonlat = True, color = 'k', fontweight = 'bold')
hp.projtext(95, 25, 'K_5 off', lonlat = True, color = 'k', fontweight = 'bold')

hp.projtext(175, 30, '30', lonlat=True, fontweight = 'bold')
hp.projtext(165, 60, '60', lonlat=True, fontweight = 'bold')
#hp.projtext(5, -5, 'Longitude (deg)', lonlat=True)
hp.projtext(178, -33.5, '-30', lonlat=True, fontweight = 'bold')
hp.projtext(175, -63.5, '-60', lonlat=True, fontweight = 'bold')


plt.show()


