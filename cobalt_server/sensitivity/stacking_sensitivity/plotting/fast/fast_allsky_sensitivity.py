import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from skylab import utils
from scipy.stats import chi2, norm, exponweib
from icecube.umdtools import cache,misc
from icecube import icetray, dataclasses, histlite
from optparse import OptionParser
import argparse
fitfun = utils.FitDeltaChi2()

#couple fcns to get the plot to match csky
pad = 0.14
def icprelim (fig, x=pad + .02, y=1 - pad - .02, **kw):
    """Mark a figure as preliminary."""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text (x, y, 'IceCube Preliminary', **kw)

def getfig (fignum=None, aspect=None, width=None, figsize=None):
    aspect = aspect or 4/3.
    width = width or 7
    if figsize is None:
        figsize = (width, width / aspect)
    out = plt.figure (num=fignum, figsize=figsize)
    plt.clf ()
    return out

def pfig (*a, **kw):
    fig = getfig (*a, **kw)
    ax = fig.add_subplot (111)
    return fig, ax

##This script should be able to plot the bckg TS, along with the TS for sensdisc, given any catalog and number of years - we'll ask cat, n as args.

##arguments:

parser = OptionParser (usage = '%prog [options]')


parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years of data')

opts, args = parser.parse_args ()
years = opts.n

##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

alldata = '/data/user/brelethford/Output/fast/allsky_sensitivity/{}yr/'.format(str(years))

results = []
for dec in os.listdir(alldata):
  if dec.endswith('000'):
    try:
      results.append((float(dec), cache.load(alldata+dec+'/results.array')))
    except:
      print ('No results available for dec at {}'.format(dec))

results.sort()

decs, sens, disc, sens_err, disc_err = [], [], [], [], []
for result in results:
  decs.append(result[0])
  sens.append(result[1]['sensitivity_flux']*1e3)
  disc.append(result[1]['discovery_flux']*1e3)
  sens_err.append(result[1]['sensitivity_flux_error']*1e3)
  disc_err.append(result[1]['discovery_flux_error']*1e3)


my_sens = np.array(zip(list(np.sin(np.radians(decs))),sens))
my_disc = np.array(zip(list(np.sin(np.radians(decs))),disc))

official_folder = '/home/brelethford/stefan_data'
if years == 1:
   real_sens = np.genfromtxt('{0}/sens_86.csv'.format (official_folder), delimiter = ',')
   real_disc = np.genfromtxt('{0}/disc_86.csv'.format (official_folder), delimiter = ',')
elif years==4:
   real_sens = np.genfromtxt('{0}/sens4yr.csv'.format (official_folder), delimiter = ',')
   real_disc = np.genfromtxt('{0}/disc4yr.csv'.format (official_folder), delimiter = ',')
elif years==7:
   real_sens = np.load('{0}/track_sens.npy'.format (official_folder))
   real_disc = np.load('{0}/track_disc.npy'.format (official_folder))

colors = ['blue', 'red']
lw = 2
ls = '-'
#plot prep
#misc.tex_mpl_rc (False)

fig = plt.figure(100)
fig.clf()
gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
ax = plt.subplot (gs[0])
rax = plt.subplot (gs[1], sharex=ax)

my_label = 'skylab / $E^{{-{0}}}$'.format (2)

#my results
msd, ms2 = my_sens.T
mdd, md2 = my_disc.T

#jake's results
if years == 7:
  sd, s2 = np.sin(real_sens['dec']), real_sens['2']
  dd, d2 = np.sin(real_disc['dec']), real_disc['2']
else:
  sd, s2 = real_sens.T
  dd, d2 = real_disc.T

## Now to plot. ##
if years == 7:
  sys = 1.11
elif years == 4:
  sys = 1.21
else:
  sys = 1.
#mine
ms2 = ms2*sys
md2 = md2*sys
ax.semilogy (msd, ms2, label=my_label + ' - sens', lw=lw, ls=ls, color=colors[0])
ax.fill_between(msd,ms2-sens_err,ms2+sens_err, alpha = 0.2, color = colors[0])

ax.semilogy (mdd, md2, label=my_label + ' - disc', lw=lw, ls=ls, color=colors[1])
ax.fill_between(mdd,md2-disc_err,md2+disc_err, alpha = 0.2, color = colors[1])


#official

ax.semilogy (sd, s2,
           color=colors[0], lw=lw, ls=':', alpha = 0.5,
           label=r'official results / $E^{-2}$ - sens')

ax.semilogy (dd, d2,
           color=colors[1], lw=lw, ls=':', alpha = 0.5,
           label=r'official results / $E^{-2}$ - disc')

ax.set_ylabel (r'$E^2 '
        '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
        '\cdot dN/dE\,\,\,'
        '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
ax.set_ylim (10**-13., 10**-9.)
ymin, ymax = ax.get_ylim ()
ax.set_yticks (np.logspace (-13, -9, 5))


ax.text (.95, 10**-12.9, 'North', size='small', ha='right', va='bottom')
ax.text (-.95, 10**-12.9, 'South', size='small', ha='left', va='bottom')
leg = ax.legend (
    loc='upper right',
    prop=propsmall, ncol=2, handlelength=2.2,
)

#some Ratio plots
rat2 = np.interp (sd, msd, ms2) / s2
rat2disc = np.interp (dd, mdd, md2) / d2
rax.plot (sd, rat2, color=colors[0])
rax.plot (dd, rat2disc, color=colors[1])

rax.set_xlabel (r'$\sin(\delta)$')
rax.set_ylabel ('Drexel / Official')
rax.set_xlim (-1, 1)
rax.set_ylim (0.5,1.5)

ax.grid ()
rax.grid ()
fig.subplots_adjust (top=.93)
icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')

misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/fast/plots/allsky/{}yr/'.format(str(years)))
fig.savefig ('/data/user/brelethford/AGN_Core/Plots/fast/plots/allsky/{}yr/sensitivity.png'.format(str(years)))
fig.savefig ('/data/user/brelethford/AGN_Core/Plots/fast/plots/allsky/{}yr/sensitivity.pdf'.format(str(years)))
misc.ensure_dir('/home/brelethford/public_html/skylab/fast/plots/allsky/{}yr/'.format(str(years)))
fig.savefig('/home/brelethford/public_html/skylab/fast/plots/allsky/{}yr/sensitivity.png'.format(str(years)))
fig.savefig('/home/brelethford/public_html/skylab/fast/plots/allsky/{}yr/sensitivity.pdf'.format(str(years)))

