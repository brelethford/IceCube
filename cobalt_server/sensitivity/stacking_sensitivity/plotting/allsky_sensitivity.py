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

parser.add_option ('--datatype', dest = 'datatype', type = str,
                default = 'npz', metavar = 'DATA',
                help = 'Chooses my datatype or npz datatype')

opts, args = parser.parse_args ()
years = opts.n
datatype = opts.datatype

##Make plots prettier
misc.tex_mpl_rc()
w=4
propsmall = mpl.font_manager.FontProperties (size='small')
propxsmall = mpl.font_manager.FontProperties (size='x-small')

alldata = '/data/user/brelethford/Output/{0}/allsky_sensitivity/{1}yr/'.format(datatype,str(years))

sens_data = alldata + 'sens/'
disc_data = alldata + 'disc/'

decbins=35
decs = np.linspace(-85.,85.,decbins)
sens_decfolders = [sens_data+'dec_{0:4f}/'.format(dec) for dec in decs]
disc_decfolders = [disc_data+'dec_{0:4f}/'.format(dec) for dec in decs]


bins = 40
range = (0.0,40.0)

sens = []
disc = []

for sensfolder, discfolder in zip(sens_decfolders,disc_decfolders):
    sens.append(cache.load(sensfolder+'gamma2.0.array')[0]['flux'][0] * 1e3)
    disc.append(cache.load(discfolder+'gamma2.0.array')[0]['flux'][0] * 1e3)

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

ax.semilogy (mdd, md2, label=my_label + ' - disc', lw=lw, ls=ls, color=colors[1])

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

misc.ensure_dir ('/data/user/brelethford/AGN_Core/Plots/allsky/{}/sens/{}yr/'.format(datatype,str(years)))
fig.savefig ('/data/user/brelethford/AGN_Core/Plots/allsky/{}/sens/{}yr/sensitivity.png'.format(datatype,str(years)))
misc.ensure_dir('/home/brelethford/public_html/skylab/allsky/{}/sens/{}yr/'.format(datatype,str(years)))
fig.savefig('/home/brelethford/public_html/skylab/allsky/{}/sens/{}yr/sensitivity.png'.format(datatype,str(years)))

