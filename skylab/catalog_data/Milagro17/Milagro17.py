import numpy as np
from icecube.umdtools import misc, cache

ra  = np.zeros(17)
dec = np.zeros(17)
w   = np.zeros(17)
ext = np.zeros(17)

# 17 milagro sources - ask mrichman where he got the extensions
ra[ 0] =  5.31766916; dec[ 0] =  0.64280476; w[ 0] = 1.00; ext[0] = np.radians(1.1)
ra[ 1] =  5.00350990; dec[ 1] =  0.10175270; w[ 1] = 1.00; ext[1] = np.radians(0)
ra[ 2] =  5.37439237; dec[ 2] =  0.70982541; w[ 2] = 1.00; ext[2] = np.radians(3)
ra[ 3] =  5.42762491; dec[ 3] =  0.63355452; w[ 3] = 1.00; ext[3] = np.radians(0)
ra[ 4] =  5.37125077; dec[ 4] =  0.63739424; w[ 4] = 1.00; ext[4] = np.radians(0)
ra[ 5] =  4.94137618; dec[ 5] =  0.00890118; w[ 5] = 1.00; ext[5] = np.radians(0)
ra[ 6] =  1.46049152; dec[ 6] =  0.38589230; w[ 6] = 1.00; ext[6] = np.radians(0)
ra[ 7] =  1.64689268; dec[ 7] =  0.39392081; w[ 7] = 1.00; ext[7] = np.radians(0)
ra[ 8] =  1.70955000; dec[ 8] =  0.18448130; w[ 8] = 1.00; ext[8] = np.radians(0)
ra[ 9] =  1.71251706; dec[ 9] =  0.30316369; w[ 9] = 1.00; ext[9] = np.radians(1.3)
ra[10] =  4.90507333; dec[10] = -0.06265732; w[10] = 1.00; ext[10] = np.radians(0)
ra[11] =  4.97436290; dec[11] =  0.06894051; w[11] = 1.00; ext[11] = np.radians(0)
ra[12] =  5.07489387; dec[12] =  0.24766222; w[12] = 1.00; ext[12] = np.radians(0)
ra[13] =  5.21172768; dec[13] =  0.50003683; w[13] = 1.00; ext[13] = np.radians(0)
ra[14] =  5.22778471; dec[14] =  0.50265482; w[14] = 1.00; ext[14] = np.radians(0)
ra[15] =  5.33023554; dec[15] =  0.70581115; w[15] = 1.00; ext[15] = np.radians(0)
ra[16] =  5.88490117; dec[16] =  1.06761790; w[16] = 1.00; ext[16] = np.radians(0)



#Finally let's save these values.

params = {'ra':list(ra),'dec':list(dec),'weight':list(w), 'ext':list(ext)}

cache.save(params, '/data/user/brelethford/Data/Milagro17/pickle/params.pickle')



