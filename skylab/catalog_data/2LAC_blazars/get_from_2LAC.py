from icecube.umdtools import misc, cache
import numpy as np

with open('/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/2LAC_blazars/raw2LAC.txt') as f:
    lines = map (str.strip, f.readlines())

import re
regex = r'\| ([-0-9.]+) \| ([-0-9.]+) \|.*\| [(<] ?([-0-9.]+)'
ra, dec, flux = np.array ([
    map (float, re.search (regex, line).groups())
    for line in lines[1:]
]).T

regex = r'\| [-0-9.]+ \| [-0-9.]+ \| ([\w ]+) \|.*\| [(<] ?[-0-9.]+'
src_class = np.array ([
    re.search (regex, line).groups()
    for line in lines[1:]
]).T[0]

mask = np.sum ([src_class == kind for kind in 'BL Lac/FSRQ/AGU'.split('/')], axis=0) > 0

params = {'src_dec':dec[mask],'src_ra':ra[mask],'flux':flux[mask],'class':src_class[mask]}

cache.save(params, '/data/user/brelethford/Data/2LAC/pickle/params.pickle')

