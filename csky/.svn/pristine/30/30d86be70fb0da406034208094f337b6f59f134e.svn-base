#!/usr/bin/env python
# csky.py

from __future__ import print_function

from collections import defaultdict
import copy
from glob import glob
from itertools import izip
import os
import re
import socket
import sys
import time

import matplotlib as mpl
if socket.gethostname ()[:5] == 'cobol':
    mpl.use ('Agg')
elif socket.gethostname ()[:6] == 'condor':
    mpl.use ('Agg')
elif socket.gethostname () == 'squirrel':
    mpl.use ('Agg')
else:
    #print ('Using TkAgg')
    pass
import healpy
import matplotlib.pyplot as plt
import numpy as np
pi = np.pi
from optparse import OptionParser
from scipy import optimize, stats
import tables

from icecube import histlite
from icecube.umdtools.apps import Timed, Driver, command
from icecube.umdtools.conf import Configurator
from icecube.umdtools.submitter import Submitter
from icecube.umdtools.arrays import Arrays, combined_arrays
from icecube.umdtools.vars_class import Vars
from icecube.umdtools import cache, misc, quiet_healpy
ensure_dir = misc.ensure_dir

from icecube.csky import ana, pdf, hyp, llh, trial, bk, cat, coord, dists

no_emojify = lambda *a: a[0]

try:
    from emojify import emojify
except:
    emojify = no_emojify

if socket.gethostname () not in ('ordi', 'zgaskami', 'condor00'):
    emojify = no_emojify

job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

propsmall = mpl.font_manager.FontProperties (size=11)
propsmaller = mpl.font_manager.FontProperties (size=8)

def prush (*a, **kw):
    """Print and flush."""
    f = kw.get ('file', sys.stdout)
    print (*a, **kw)
    f.flush ()

def run (command):
    prush (command)
    os.system (command)

def getting (filename, *a, **kw):
    prush ('<-> {0} ...'.format (filename))
    return cache.get (filename, *a, **kw)

def loading (filename):
    prush ('<- {0} ...'.format (filename))
    return cache.load (filename)

def pattern_load (pattern):
    return [loading (filename)
            for filename in sorted (glob (pattern))]

def saving (obj, filename):
    prush ('-> {0} ...'.format (filename))
    return cache.save (obj, filename)

def savefig (fig, outdir, namebase, exts='png pdf', specialpdf=False):
    for ext in exts.split ():
        #if ext in ('eps', 'pdf'):
        #    plt.rc ('text', usetex=True)
        #else:
        #    plt.rc ('text', usetex=False)
        if specialpdf and ext == 'pdf':
            w, h = fig.get_figwidth (), fig.get_figheight ()
            if w > 10:
                factor = .70
            else:
                factor = .52
            fig.set_size_inches (factor * w, factor * h)
        fig.savefig ('{0}/{1}.{2}'.format (
            outdir, namebase, ext))
        if specialpdf and ext == 'pdf':
            fig.set_size_inches (w, h)

def savingfig (f, d, n, *a, **kw):
    prush ('-> {0}/{1} ...'.format (d, n))
    savefig (f, d, n, *a, **kw)
    #Need the same directory names as in condor:
    public_dir = '/home/brelethford/public_html/csky/' + d.split('brelethford/')[-1] + '/'
    print (public_dir)
    #eval scp '{0}/{1} ...'.format (d, n) brelethford@cobalt06:public_dir

def with_trailing_slash (dirname):
    dirname = dirname + '/'
    while dirname[-2:] == '//':
        dirname = dirname[:-1]
    return dirname

def without_trailing_slash (dirname):
    while dirname[-1] == '/':
        dirname = dirname[:-1]
    return dirname

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

def pjobdir (job_dir):
    prush ('\nJob dir is {0} .\n'.format (job_dir))

pad = 0.14
def icprelim (fig, x=pad + .02, y=1 - pad - .02, **kw):
    """Mark a figure as preliminary."""
    if 'color' not in kw:
        kw['color'] = 'red'
    if 'weight' not in kw:
        kw['weight'] = 'bold'
    fig.text (x, y, 'IceCube Preliminary', **kw)

soft_colors = ['#004466', '#c65353', '#5aeac0', '#dd9388', '#caca68']

class Csky6yr (Timed, Driver):

    def __init__ (self):
        Timed.__init__ (self)
        Driver.__init__ (self)

    def run (self, arglist=[]):

        usage = '%prog {[options]} [commands]\n' + self._command_help
        self.parser = parser = OptionParser (usage=usage)

        parser.add_option ('-c', '--conf-dirs', dest='conf_dirs',
                default=os.path.abspath ('conf'), metavar='DIRS',
                help='load configuration from comma-separated list of DIRS')

        parser.add_option ('--nside', dest='nside',
                default=16, type=int, metavar='NSIDE',
                help='use NSIDE with healpy')

        parser.add_option ('--sigma', dest='sigma',
                default=0, type=int, metavar='N',
                help='handle N-sigma calculations')

        parser.add_option ('--beta', dest='beta',
                default=None, type=float, metavar='BETA',
                help='must surpass threshold in BETA fraction of trials')

        parser.add_option ('--test-ext', dest='test_ext',
                default=0., type=float, metavar='EXT',
                help='test point has extension EXT degrees')

        parser.add_option ('--test-dec', dest='test_dec',
                default=None, type=float, metavar='dec',
                help='test point DEC in degrees')

        parser.add_option ('--i-source', dest='i_source',
                default=None, type=int, metavar='SOURCE',
                help='which catalog source to use')

        parser.add_option ('--src-ext', dest='src_ext',
                default=0., type=float, metavar='EXT',
                help='source point has extension EXT degrees')

        parser.add_option ('--src-gamma', dest='src_gamma',
                default=None, type=float, metavar='GAMMA',
                help='source has spectral index GAMMA')

        parser.add_option ('--src-thresh', dest='src_thresh',
                default=2, type=float, metavar='THRESH',
                help='source has spectral threshold THRESH')

        parser.add_option ('--src-cutoff', dest='src_cutoff',
                default=np.inf, type=float, metavar='CUTOFF',
                help='source has spectral cutoff CUTOFF')

        parser.add_option ('--n-inj', dest='n_inj',
                default=0, type=int, metavar='N',
                help='inject N signal events')

        parser.add_option ('--n-jobs', dest='n_jobs',
                default=10., type=float, metavar='N',
                help='perform N jobs (with --n-trials each)')

        parser.add_option ('--n-trials', dest='n_trials',
                default=100., type=float, metavar='N',
                help='perform N trials')

        parser.add_option ('--n-multi', dest='n_multi',
                default=2, type=int, metavar='NTEST',
                help='number of sources for multi-source test')

        parser.add_option ('--seed', dest='seed',
                default=0, type=int, metavar='SEED',
                help='initialize RNG with SEED')

        parser.add_option ('--just-diff', dest='just_diff',
                default=False, action='store_true',
                help='only do bg trials for diff sens dec(s)')

        parser.add_option ('--blacklist', dest='blacklist',
                default='cobol61,cobol65', help='nodes to avoid')

        parser.add_option ('--diff-ra', dest='diff_ra',
                default=False, action='store_true',
                help='toggles different ra for multi stacking test')

        parser.add_option ('--oneyear', dest='oneyear',
                default=False, action='store_true',
                help='toggles single year analysis treatment')

        parser.add_option ('--do-ext', dest='do_ext',
                default=False, action='store_true',
                help='do trials for extended sources')

        parser.add_option ('--sys-total', dest='sys_total',
                #default=1.21, type=float, metavar='FACTOR',
                default=1.11, type=float, metavar='FACTOR',
                help='plot fluxes including FACTOR systematic uncertainty')

        self.opts, self.commands = opts, commands = \
                parser.parse_args (arglist if arglist else sys.argv[1:])

        self.conf = Configurator (*self.opts.conf_dirs.split (','))
        self.mode = Vars ()
        self.mode.str = self.conf.root.mode
        self.mode.words = self.mode.str.split ('/')
        #self.mode.break_fit = 'break' in self.mode.words

        #self.go (commands, announcement=
        #         lambda s: emojify (
        #             ':penguin: :penguin: :penguin: '
        #             '{0} :penguin: :penguin: :penguin:'.format (s),
        #             False))

        self.go (commands, announcement=
                 lambda s: '<<< {} >>>'.format (s))

    # properties

    @property
    def root_dir (self):
        if socket.gethostname () in ('condor00'):
            return self.conf.root.root_dir
        else:
            return self.conf.root.cluster_root_dir

    @property
    def mode_dir (self):
        return ensure_dir ('{0}/{1}'.format (self.root_dir, self.mode.str))

    @property
    def mode_dir86 (self):
        return ensure_dir ('{0}/{1}/oneyear'.format (self.root_dir, self.mode.str))

    @property
    def cdata (self):
        try:
            return self._cdata
        except:
            self._cdata = loading ('{0}/cdata.vars'.format (self.mode_dir))
            return self._cdata

    @property
    def datasets (self):
        try:
            return self._datasets
        except:
            dataset40 = loading ('{0}/IC40.dataset'.format (self.mode_dir))
            dataset59 = loading ('{0}/IC59.dataset'.format (self.mode_dir))
            dataset79 = loading ('{0}/IC79.dataset'.format (self.mode_dir))
            dataset86 = loading ('{0}/IC86.dataset'.format (self.mode_dir))
            dataset3yr = loading ('{0}/IC3yr.dataset'.format (self.mode_dir))
            datasetMESE = loading ('{0}/MESE.dataset'.format (self.mode_dir))
            self._datasets = [dataset40,dataset59,dataset79,dataset86,dataset3yr,datasetMESE]
            return self._datasets

    @property
    def analysis (self):
        try:
            return self._analysis
        except:
            self._analysis = loading ('{0}/ps7yrs.analysis'.format (self.mode_dir))
            return self._analysis

    @property
    def dataset86 (self):
        try:
            return self._dataset86
        except:
            dataset86 = loading ('{0}/IC86_nopull.dataset'.format (self.mode_dir))
            self._dataset86 = [dataset86]
            return self._dataset86

    @property
    def analysis86 (self):
        try:
            return self._analysis86
        except:
            self._analysis86 = loading ('{0}/psIC86.analysis'.format (self.mode_dir86))
            return self._analysis86

    @property
    def ns_bounds (self):
        return (0, .999)


    @property
    def bg_tsds (self):
        if hasattr (self, '_bg_tsds'):
            return self._bg_tsds
        try:
            self._bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))
        except:
            self._bg_tsds = self.collect_bg_trials ()
        return self._bg_tsds

    @property
    def bg_tsds86 (self):
        if hasattr (self, '_bg_tsds86'):
            return self._bg_tsds86
        try:
            self._bg_tsds86 = loading ('{0}/bg_tsds.dict'.format (self.mode_dir86))
        except:
            self._bg_tsds86 = self.collect_bg_trials ()
        return self._bg_tsds86

    @property
    def sig_int (self):
        try:
            return loading ('{0}/sig_int.dict'.format (self.mode_dir))
        except:
            return self.collect_sig_int ()

    @property
    def sig_int86 (self):
        try:
            return loading ('{0}/sig_int.dict'.format (self.mode_dir86))
        except:
            return self.collect_sig_int ()

    def skymap_ra_dec (self, nside=None):
        if nside is None:
            nside = self.opts.nside
        zeniths, azimuths = np.array ([
            healpy.pix2ang (nside, i_pix)
            for i_pix in xrange (healpy.nside2npix(nside))]).T
        return azimuths, pi/2 - zeniths

    @property
    def blacklist (self):
        return [s+'.private.pa.umd.edu'
                for s in self.opts.blacklist.split (',')]

    # setting up raw data


    @command
    def build_cdata (self):
        """
        Build PS-only data.
        """
        cdata = Vars ()


        #ids_exp = ['{}'.format (s)
        #       for s in 'IC40 IC59 IC79 IC86I epinat_3yr'.split()]
        ids_exp = ['{}'.format (s)
               for s in 'IC40 IC59 IC79/noMESE IC86I/noMESE 3yr/noMESE MESE nopull/IC86I'.split()]
        filenames_exp = ['{}/psdata/{}/exp.pickle'.format (self.root_dir, i)
                     for i in ids_exp]

        #ids_mc = ['{}'.format (s)
        #       for s in 'IC40 IC59 IC79 IC86I epinat_3yr'.split()]
        ids_mc = ['{}'.format (s)
               for s in 'IC40 IC59 IC79/noMESE IC86I/noMESE 3yr/noMESE MESE nopull/IC86I'.split()]
        filenames_mc = ['{}/psdata/{}/mc.pickle'.format (self.root_dir, i)
                     for i in ids_mc]

        prush ('Loading MC and exp datasets...')
        contents_exp = map (cache.load, filenames_exp)
        contents_mc = map (cache.load, filenames_mc)

        nu40 = cdata.nu40 = Arrays ()
        nu59 = cdata.nu59 = Arrays ()
        nu79 = cdata.nu79 = Arrays ()
        nu86 = cdata.nu86 = Arrays ()
        nu86_nopull = cdata.nu86_nopull = Arrays ()
        nu3yr = cdata.nu3yr = Arrays ()
        nuMESE = cdata.nuMESE = Arrays ()
        nus = [nu40,nu59,nu79,nu86,nu3yr,nuMESE,nu86_nopull] 

        for nu,mc in zip(nus,contents_mc): 
            nu.true_energy =  mc['trueE']
            nu.true_zenith =  pi/2.+mc['trueDec']
            nu.true_azimuth = mc['trueRa']
            nu.energy =       10**mc['logE']
            nu.zenith =       np.pi/2.+np.arcsin(mc['sinDec'])
            nu.azimuth =      mc['ra']
            nu.oneweight =    mc['ow']
            nu.sigma =        mc['sigma']

        mu40 = cdata.mu40 = Arrays ()
        mu59 = cdata.mu59 = Arrays ()
        mu79 = cdata.mu79 = Arrays ()
        mu86 = cdata.mu86 = Arrays ()
        mu86_nopull = cdata.mu86_nopull = Arrays ()
        mu3yr = cdata.mu3yr = Arrays ()
        muMESE = cdata.muMESE = Arrays ()
        mus = [mu40,mu59,mu79,mu86,mu3yr,muMESE,mu86_nopull] 

        for mu,exp in zip(mus,contents_exp): 
            mu.energy =   10**exp['logE']
            mu.zenith =   np.pi/2.+np.arcsin(exp['sinDec'])
            mu.azimuth =  exp['ra']
            mu.sigma =    exp['sigma']

        #I think I need to give a bckg weight explicitly...
        #mu.weight = np.ones_like(mu.energy)
        for nu,mu in zip(nus,mus):
          nu.apply_cut ((1e1 < nu.energy) & (nu.energy < 1e10))
          mu.apply_cut ((1e1 < mu.energy) & (mu.energy < 1e10))
          for d in [nu, mu]:
            d.dec = d.zenith - pi/2
            d.ra = d.azimuth
            if 'true_zenith' in d:
                d.true_dec = d.true_zenith - pi/2
                d.true_ra = d.true_azimuth
        saving (cdata, '{0}/cdata.vars'.format (self.mode_dir))

    @command
    def setup_datasets (self):
        """Create Dataset."""
        nu40 = self.cdata.nu40
        nu59 = self.cdata.nu59
        nu79 = self.cdata.nu79
        nu86_nopull = self.cdata.nu86_nopull
        nu86 = self.cdata.nu86
        nu3yr = self.cdata.nu3yr
        nuMESE = self.cdata.nuMESE
        mu40 = self.cdata.mu40
        mu59 = self.cdata.mu59
        mu79 = self.cdata.mu79
        mu86_nopull = self.cdata.mu86_nopull
        mu86 = self.cdata.mu86
        mu3yr = self.cdata.mu3yr
        muMESE = self.cdata.muMESE
        livetime40  = 375.539 * 86400
        livetime59  = 348.138 * 86400
        livetime79  = 315.506 * 86400
        livetime86 = 332.61 *  86400
        livetime86II = 330.38 * 86400
        livetime86III = 359.95 * 86400
        livetime86IV = 367.21 * 86400
        livetime3yr = livetime86II+livetime86III+livetime86IV
        livetimeMESE = 1715 * 86400

        dataset40 = ana.Dataset ('IC40', livetime40, sig=nu40, data=mu40)
        dataset59 = ana.Dataset ('IC59', livetime59, sig=nu59, data=mu59)
        dataset79 = ana.Dataset ('IC79', livetime79, sig=nu79, data=mu79)
        dataset86 = ana.Dataset ('IC86', livetime86, sig=nu86, data=mu86)
        dataset86_nopull = ana.Dataset ('IC86_nopull', livetime86, sig=nu86_nopull, data=mu86_nopull)
        dataset3yr = ana.Dataset ('IC3yr', livetime3yr, sig=nu3yr, data=mu3yr)
        datasetMESE = ana.Dataset ('MESE', livetimeMESE, sig=nuMESE, data=muMESE)

        saving (dataset40, '{}/IC40.dataset'.format (self.mode_dir))
        saving (dataset59, '{}/IC59.dataset'.format (self.mode_dir))
        saving (dataset79, '{}/IC79.dataset'.format (self.mode_dir))
        saving (dataset86, '{}/IC86.dataset'.format (self.mode_dir))
        saving (dataset3yr, '{}/IC3yr.dataset'.format (self.mode_dir))
        saving (datasetMESE, '{}/MESE.dataset'.format (self.mode_dir))
        saving (dataset86_nopull, '{}/IC86_nopull.dataset'.format (self.mode_dir))

    @command
    def setup_ana (self):
        """
        Set up main analysis objects.
        """

        n_year = self.conf.root.n_year
        #nu = self.cdata.nu
        #atm = self.cdata.atm
        #diffuse = self.cdata.diffuse
        #mu = self.cdata.mu

        #fit = 'nospline' not in self.mode.words
        fit = False

        datasets = self.datasets
        #this actually references the IC86 without the factor of 1.1774
        dataset86 = self.dataset86

        analysis = ana.Analysis (
            datasets,
            dec_kw=dict (
                bins=41, range=(-1, 1),
                fit=fit
            ),
            energy_kw=dict (
                sindec_bins=41, logenergy_bins=36,
                logenergy_range=(1, 10),
                gammas=np.arange (1, 4.01, 0.125),
                fit=fit
            )#,
            #angres_kw=dict (
            #    apr_type=pdf.Kent1DAngResParameterization,
            #    sindec_bins=41, logenergy_bins=41,
            #    logenergy_range=(1, 10),
            #    fit=fit
            #)
        )
        analysis86 = ana.Analysis (
            dataset86,
            dec_kw=dict (
                bins=41, range=(-1, 1),
                fit=fit
            ),
            energy_kw=dict (
                sindec_bins=41, logenergy_bins=36,
                logenergy_range=(1, 10),
                gammas=np.arange (1, 4.01, 0.125),
                fit=fit
            )
        )

        saving (analysis, '{}/ps7yrs.analysis'.format (self.mode_dir))
        saving (analysis86, '{}/psIC86.analysis'.format (self.mode_dir86))

    @command
    def angres (self):
        """
        Make '2D' angular resolution plots.
        """
        if isinstance (self.analysis.apr, pdf.AngResParameterization):
            self.angres_2d ()
        else:
            self.angres_kent ()

    def angres_2d (self):
        """Make angular resolution plots."""
        misc.tex_mpl_rc ()
        ana = self.analysis
        plot_dir = misc.ensure_dir ('{0}/plots/pdfs'.format (self.mode_dir))

        hazimuth = ana.apr.hazimuth/pi*180
        hzenith = ana.apr.hzenith/pi*180

        fig = plt.figure (figsize=(12,5))
        axA = fig.add_subplot (1,2,1)
        fig.subplots_adjust (left=.06, right=.98)
        histlite.plot2d (axA, hazimuth,
                         vmin=0, vmax=180, cmap='Blues')
        histlite.label2d (axA, hazimuth, '.1f', size=8)
        axA.set_title (r'inferred azimuth uncertainty~($^\circ$)')
        axA.set_xlabel ('cos(reco zenith)')
        axA.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')

        axZ = fig.add_subplot (1,2,2)
        histlite.plot2d (axZ, hzenith,
                         vmin=0, vmax=180, cmap='Blues')
        histlite.label2d (axZ, hzenith, '.1f', size=8)
        axZ.set_title (r'inferred zenith uncertainty~($^\circ$)')
        axZ.set_xlabel ('cos(reco zenith)')
        axZ.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')
        savingfig (fig, plot_dir, 'angres_inferred')

        fig = plt.figure (figsize=(12,5))
        fig.subplots_adjust (left=.06, right=.98)
        axA = fig.add_subplot (1, 2, 1)
        axZ = fig.add_subplot (1, 2, 2)
        cz = np.linspace (-.99, .99, 100)
        le = np.linspace (3.01, 7.99, 100)
        CZ, LE = np.meshgrid (cz, le)
        axA.pcolormesh (CZ, LE, hazimuth.spline_fit() (CZ, LE),
                        cmap='Blues', vmin=0, vmax=180)
        axZ.pcolormesh (CZ, LE, hzenith.spline_fit() (CZ, LE),
                        cmap='Blues', vmin=0, vmax=180)

        axA.set_title (r'inferred azimuth uncertainty~($^\circ$)')
        axA.set_xlabel ('cos(reco zenith)')
        axA.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')
        axZ.set_title (r'inferred zenith uncertainty~($^\circ$)')
        axZ.set_xlabel ('cos(reco zenith)')
        axZ.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')
        plot_dir = misc.ensure_dir ('{0}/plots/pdfs'.format (self.mode_dir))
        savingfig (fig, plot_dir, 'angres_inferred_smooth')


    @command
    def angres_kent (self):
        """Make Kent angular resolution plots."""
        misc.tex_mpl_rc (True)
        plt.rc ('font', size=16)
        apr = self.analysis.apr['ps7yrs']
        plot_dir = misc.ensure_dir ('{0}/plots/pdfs'.format (self.mode_dir))

        h = apr.h/pi*180

        fig = plt.figure (figsize=(8,8))
        ax = fig.add_subplot (111)
        histlite.plot2d (ax, h,
                         vmin=0, vmax=180, cmap='Blues')
        histlite.label2d (ax, h, '.1f', size=4)
        ax.set_title (r'inferred uncertainty~($^\circ$)')
        ax.set_xlabel ('cos(reco zenith)')
        ax.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')

        savingfig (fig, plot_dir, 'angres_inferred')

        if not apr.fit:
            return

        fig = plt.figure (figsize=(8,8))
        ax = fig.add_subplot (111)
        cz = np.linspace (-.99, .99, 100)
        le = np.linspace (3.01, 7.99, 100)
        CZ, LE = np.meshgrid (cz, le)
        ax.pcolormesh (CZ, LE, h.spline_fit() (CZ, LE),
                       cmap='Blues', vmin=0, vmax=180)

        ax.set_title (r'inferred uncertainty~($^\circ$)')
        ax.set_xlabel ('cos(reco zenith)')
        ax.set_ylabel (r'$\log_{10}(\mathrm{reco}~E/\mathrm{GeV})$')
        savingfig (fig, plot_dir, 'angres_inferred_smooth')

    # background and signal trials

    @command
    def submit_do_bg_trials (self):
        """Submit bg TSDist jobs to cluster."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_bg_trials/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        if self.opts.just_diff:
            dec_degs = [-60]
        else:
            dec_degs = np.arange (-90, 90, 3)
        for test_ext in [0.]:
            for dec_deg in dec_degs:
                dec_deg = max (-89, min (89, dec_deg))
                for i_job in xrange (int (self.opts.n_jobs)):
                  if self.opts.oneyear:
                    command = '{} {} do_bg_trials --oneyear --conf-dirs={} ' \
                            ' --n-trials={}' \
                            ' --test-dec={:+07.2f}' \
                            ' --test-ext={:04.1f}' \
                            ' --seed={}'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.n_trials,
                                dec_deg,
                                test_ext,
                                i_job)
                  else:
                    command = '{} {} do_bg_trials --conf-dirs={} ' \
                            ' --n-trials={}' \
                            ' --test-dec={:+07.2f}' \
                            ' --test-ext={:04.1f}' \
                            ' --seed={}'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.n_trials,
                                dec_deg,
                                test_ext,
                                i_job)
                  label = 'do_bg_trials__dec_{:+07.2f}__ext__{:04.1f}__seed_{:08d}'.format (
                            dec_deg, test_ext, i_job)
                  commands.append (command)
                  labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def submit_do_bg_trials_multi (self):
        """Submit bg TSDist jobs for double (or more) test to cluster."""
        job_root = self.conf.root.cluster_job_dir
        n_multi = self.opts.n_multi
        print ("{} sources used in multi-test".format(n_multi))
        if self.opts.diff_ra:
          job_dir = '{0}/do_bg_trials_multi_diff/{1}'.format (job_root, job_id)
        else:
          job_dir = '{0}/do_bg_trials_multi_same/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        if self.opts.just_diff:
            dec_degs = [-60]
        else:
            dec_degs = np.arange (-90, 91, 15)
        for test_ext in [0.]:
            for dec_deg in dec_degs:
                dec_deg = max (-89, min (89, dec_deg))
                for i_job in xrange (int (self.opts.n_jobs)):
                    if self.opts.diff_ra:
                      command = '{} {} do_bg_trials_multi --diff-ra --n-multi={} ' \
                            ' --conf-dirs={} ' \
                            ' --n-trials={}' \
                            ' --test-dec={:+07.2f}' \
                            ' --test-ext={:04.1f}' \
                            ' --seed={}'.format (
                                env_shell,
                                this_script,
                                n_multi,
                                confs,
                                self.opts.n_trials,
                                dec_deg,
                                test_ext,
                                i_job)
                      label = 'do_bg_trials_multi_diff_dec_{:+07.2f}__ext__{:04.1f}__seed_{:08d}'.format (
                            dec_deg, test_ext, i_job)
                    else:
                      command = '{} {} do_bg_trials_multi --n-multi={} ' \
                            ' --conf-dirs={} ' \
                            ' --n-trials={}' \
                            ' --test-dec={:+07.2f}' \
                            ' --test-ext={:04.1f}' \
                            ' --seed={}'.format (
                                env_shell,
                                this_script,
                                n_multi,
                                confs,
                                self.opts.n_trials,
                                dec_deg,
                                test_ext,
                                i_job)
                      label = 'do_bg_trials_multi_same_dec_{:+07.2f}__ext__{:04.1f}__seed_{:08d}'.format (
                            dec_deg, test_ext, i_job)
                    commands.append (command)
                    labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def do_bg_trials (self):
        """
        Do background-only trials to get TSDists.
        """


        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        dec_deg = self.opts.test_dec
        dec = dec_deg / 180.*pi
        test_ext_deg = self.opts.test_ext
        test_ext = test_ext_deg / 180.*pi

        np.random.seed (seed)
        if self.opts.oneyear:
          data = hyp.DataHypothesis (self.analysis86)
          gammas = self.analysis86.pdfs_energy_sig['IC86_nopull'].gammas
          ps = hyp.PointSourceHypothesis (self.analysis86, dec, 0, 2)
        else:
          data = hyp.DataHypothesis (self.analysis)
          gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
          ps = hyp.PointSourceHypothesis (self.analysis, dec, 0, 2)

        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        c = tr.get_Chi2TSD (n_trials, 500)

        sm = bk.SavingModel (
                'test_ext_deg/dec_deg/tsd',
                '{:04.1f}/{:+07.2f}/{:08d}.chi2',
                )
        if self.opts.oneyear:
            filename = sm.save (
                c, '{0}/bg_tsds'.format (self.mode_dir86),
                test_ext_deg, dec_deg, seed)
	else:
            filename = sm.save (
		    c, '{0}/bg_tsds'.format (self.mode_dir),
		    test_ext_deg, dec_deg, seed)
        prush ('->', filename)

        #tsd_dir = ensure_dir (
        #    '{0}/bg_tsds/{1:04.1f}/{2:07.3f}'.format (
        #        self.mode_dir, test_ext_deg, zenith_deg
        #    ))

        #saving (tsd, '{0}/{1:08d}.tsdist'.format (tsd_dir, seed))

    @command
    def do_bg_trials_multi (self):
        """
        Do background-only trials to get TSDists for doublesource test.
        """
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        n_multi = int(self.opts.n_multi)
        dec_deg = self.opts.test_dec
        dec = dec_deg / 180.*pi
        multi_dec = [dec for i in range(n_multi)]
        test_ext_deg = self.opts.test_ext
        test_ext = test_ext_deg / 180.*pi
        np.random.seed (seed)

        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        if self.opts.diff_ra:
          multi_ra = np.linspace(0,2*pi,n_multi+1)[:-1] 
          test = 'multi_diff'
        else:
          multi_ra = [0 for i in range(n_multi)]
          test = 'multi_same'
        ps = hyp.PointSourceHypothesis (self.analysis, multi_dec, multi_ra, 2)
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        c = tr.get_Chi2TSD (n_trials, 500)

        sm = bk.SavingModel (
                'test_ext_deg/dec_deg/tsd',
                '{:04.1f}/{:+07.2f}/{:08d}.chi2',
                )
        filename = sm.save (
                c, '{0}/tests/{1}/{2}/bg_tsds'.format (self.mode_dir,test,n_multi),
                test_ext_deg, dec_deg, seed)
        prush ('->', filename)

    @command
    def do_bg_trials_northsouth (self):
        """
        Do background-only trials to get TSDists for northsouth test.
        """
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        dec_deg = [-30.,30.]
        dec = [d / 180.*pi for d in dec_deg]
        np.random.seed (seed)

        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, [0,0], 2)
        test = 'northsouth'
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        c = tr.get_Chi2TSD (n_trials, 500)

        #sm = bk.SavingModel (
        #        'test_ext_deg/dec_deg/tsd',
        #        '{:04.1f}/{:+07.2f}/{:08d}.chi2',
        #        )
        sm = bk.SavingModel (
                'tsd',
                '{:08d}.chi2',
                )
        filename = sm.save (
                c, '{0}/tests/{1}/bg_tsds'.format (self.mode_dir,test),
                seed)
        prush ('->', filename)

    @command
    def collect_bg_trials (self):
        """Collect bg_trials dict and cache it."""
        prush ('Collecting bg trials...')
        if self.opts.oneyear:
          bg_tsds = bk.get_all (
                '{0}/bg_tsds'.format (self.mode_dir86),
                '*.chi2')
          saving (bg_tsds, '{0}/bg_tsds.dict'.format (self.mode_dir86))
        else:
          bg_tsds = bk.get_all (
                '{0}/bg_tsds'.format (self.mode_dir),
                '*.chi2')
          saving (bg_tsds, '{0}/bg_tsds.dict'.format (self.mode_dir))
        return bg_tsds

    @command
    def collect_bg_trials_multi_same (self):
        """Collect bg_trials dict for multi test and cache it."""
        prush ('Collecting bg trials for multi test with same ra...')
        bg_tsds = bk.get_all (
                '{0}/tests/multi_same/{1}/bg_tsds'.format (self.mode_dir,self.opts.n_multi),
                '*.chi2')
        saving (bg_tsds, '{0}/tests/multi_same/{1}/bg_tsds.dict'.format (self.mode_dir,self.opts.n_multi))
        return bg_tsds

    @command
    def collect_bg_trials_multi_diff (self):
        """Collect bg_trials dict for multi test with different raand cache it."""
        prush ('Collecting bg trials for multi test with diff ra...')
        bg_tsds = bk.get_all (
                '{0}/tests/multi_diff/{1}/bg_tsds'.format (self.mode_dir,self.opts.n_multi),
                '*.chi2')
        saving (bg_tsds, '{0}/tests/multi_diff/{1}/bg_tsds.dict'.format (self.mode_dir,self.opts.n_multi))
        return bg_tsds

    @command
    def collect_bg_trials_northsouth (self):
        """Collect bg_trials dict for double test and cache it."""
        prush ('Collecting bg trials for double test...')
        bg_tsds = bk.get_all (
                '{0}/tests/northsouth/bg_tsds'.format (self.mode_dir),
                '*.chi2')
        saving (bg_tsds, '{0}/tests/northsouth/bg_tsds.dict'.format (self.mode_dir))
        return bg_tsds

    @command
    def submit_do_n_sig (self):
        """
        Submit n_sig jobs for some n-sigma and beta values.
        """
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_n_sig/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        spectra = [
                (2, np.inf),
                #(2.5, np.inf),
                #(2, 5),
                (3, np.inf)
        ]
        i_job = 0
        if self.opts.do_ext:
            src_exts = np.arange (0, 10.5)
        else:
            src_exts = [0]
        for src_ext in src_exts:
            test_exts = [0] if src_ext == 0 else [0, src_ext]
            #test_exts = [0]
            for test_ext in test_exts:
                for (src_gamma, src_cutoff) in spectra:
                    if self.opts.just_diff:
                        dec_degs = [150]
                    else:
                        dec_degs = np.arange (-90, 90, 3)
                    for dec_deg in dec_degs:
                        dec_deg = max (-89, min (89, dec_deg))
                        if self.opts.oneyear:
                          command = '{} {} do_n_sig --oneyear --conf-dirs={} ' \
                                ' --n-trials={}' \
                                ' --test-dec={:+07.2f}' \
                                ' --test-ext={:04.1f}' \
                                ' --src-ext={:04.1f} ' \
                                ' --src-gamma={:4.2f} ' \
                                ' --src-cutoff={} ' \
                                ' --sigma={:.0f}' \
                                ' --beta={:03.1f}' \
                                ' --seed={}'.format (
                                    env_shell,
                                    this_script,
                                    confs,
                                    self.opts.n_trials,
                                    dec_deg,
                                    test_ext,
                                    src_ext,
                                    src_gamma,
                                    src_cutoff,
                                    self.opts.sigma,
                                    self.opts.beta,
                                    i_job,
                                    )
                        else:
                          command = '{} {} do_n_sig --conf-dirs={} ' \
                                ' --n-trials={}' \
                                ' --test-dec={:+07.2f}' \
                                ' --test-ext={:04.1f}' \
                                ' --src-ext={:04.1f} ' \
                                ' --src-gamma={:4.2f} ' \
                                ' --src-cutoff={} ' \
                                ' --sigma={:.0f}' \
                                ' --beta={:03.1f}' \
                                ' --seed={}'.format (
                                    env_shell,
                                    this_script,
                                    confs,
                                    self.opts.n_trials,
                                    dec_deg,
                                    test_ext,
                                    src_ext,
                                    src_gamma,
                                    src_cutoff,
                                    self.opts.sigma,
                                    self.opts.beta,
                                    i_job,
                                    )
                        label = 'do_n_sig__' \
                                'dec_{0:+07.2f}__' \
                                'test_ext_{1:04.1f}__' \
                                'src_ext_{2:04.1f}__' \
                                'src_gamma_{3:4.2f}__' \
                                'src_cutoff_{4}__' \
                                'seed_{5:08d}'.format (
                                        dec_deg,
                                        test_ext,
                                        src_ext,
                                        src_gamma,
                                        src_cutoff,
                                        i_job)
                        commands.append (command)
                        labels.append (label)
                        i_job += 1

        submitter.submit_condor00 (
                commands, labels,
                blacklist=self.blacklist,
                max_per_interval=20)
        pjobdir (job_dir)

    @command
    def submit_do_n_sig_multi (self):
        """
        Submit n_sig jobs for some n-sigma and beta values.
        """
        job_root = self.conf.root.cluster_job_dir
        n_multi = self.opts.n_multi
        if self.opts.diff_ra:
          job_dir = '{0}/do_n_sig_multi_diff/{1}'.format (job_root, job_id)
        else:
          job_dir = '{0}/do_n_sig_multi_same/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        spectra = [
                (2, np.inf)#,
                #(2.5, np.inf),
                #(2, 5),
                #(3, np.inf),
        ]
        i_job = 0
        if self.opts.do_ext:
            src_exts = np.arange (0, 10.5)
        else:
            src_exts = [0]
        for src_ext in src_exts:
            test_exts = [0] if src_ext == 0 else [0, src_ext]
            #test_exts = [0]
            for test_ext in test_exts:
                for (src_gamma, src_cutoff) in spectra:
                    if self.opts.just_diff:
                        dec_degs = [150]
                    else:
                        dec_degs = np.arange (-90, 91, 15)
                    for dec_deg in dec_degs:
                        dec_deg = max (-89, min (89, dec_deg))
                        if self.opts.diff_ra:
                          command = '{} {} do_n_sig_multi --diff-ra --n-multi {} '\
                                '  --conf-dirs={} ' \
                                ' --n-trials={}' \
                                ' --test-dec={:+07.2f}' \
                                ' --test-ext={:04.1f}' \
                                ' --src-ext={:04.1f} ' \
                                ' --src-gamma={:4.2f} ' \
                                ' --src-cutoff={} ' \
                                ' --sigma={:.0f}' \
                                ' --beta={:03.1f}' \
                                ' --seed={}'.format (
                                    env_shell,
                                    this_script,
                                    n_multi,
                                    confs,
                                    self.opts.n_trials,
                                    dec_deg,
                                    test_ext,
                                    src_ext,
                                    src_gamma,
                                    src_cutoff,
                                    self.opts.sigma,
                                    self.opts.beta,
                                    i_job,
                                    )
                          label = 'do_n_sig_double_opp_' \
                                'dec_{0:+07.2f}__' \
                                'test_ext_{1:04.1f}__' \
                                'src_ext_{2:04.1f}__' \
                                'src_gamma_{3:4.2f}__' \
                                'src_cutoff_{4}__' \
                                'seed_{5:08d}'.format (
                                        dec_deg,
                                        test_ext,
                                        src_ext,
                                        src_gamma,
                                        src_cutoff,
                                        i_job)
                        else:
                          command = '{} {} do_n_sig_multi --n-multi {} '\
                                ' --conf-dirs={} ' \
                                ' --n-trials={}' \
                                ' --test-dec={:+07.2f}' \
                                ' --test-ext={:04.1f}' \
                                ' --src-ext={:04.1f} ' \
                                ' --src-gamma={:4.2f} ' \
                                ' --src-cutoff={} ' \
                                ' --sigma={:.0f}' \
                                ' --beta={:03.1f}' \
                                ' --seed={}'.format (
                                    env_shell,
                                    this_script,
                                    n_multi,
                                    confs,
                                    self.opts.n_trials,
                                    dec_deg,
                                    test_ext,
                                    src_ext,
                                    src_gamma,
                                    src_cutoff,
                                    self.opts.sigma,
                                    self.opts.beta,
                                    i_job,
                                    )
                          label = 'do_n_sig_double_' \
                                'dec_{0:+07.2f}__' \
                                'test_ext_{1:04.1f}__' \
                                'src_ext_{2:04.1f}__' \
                                'src_gamma_{3:4.2f}__' \
                                'src_cutoff_{4}__' \
                                'seed_{5:08d}'.format (
                                        dec_deg,
                                        test_ext,
                                        src_ext,
                                        src_gamma,
                                        src_cutoff,
                                        i_job)
                        commands.append (command)
                        labels.append (label)
                        i_job += 1

        submitter.submit_condor00 (
                commands, labels,
                blacklist=self.blacklist,
                max_per_interval=20)
        pjobdir (job_dir)

    @command
    def submit_do_n_sig_diff (self):
        """
        Submit n_sig jobs for some n-sigma and beta values.
        """
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_n_sig_diff/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        gamma = 2.
        i_job = 0
        #for zenith_deg in 120, 150:
        for zenith_deg in (150,):
            #for width in .25, .5, 1:
            for width in (.25,):
                #for logEmin in np.arange (3.+width, 7, width):
                for logEmin in np.arange (3.+width, 9, width):
                    logEmax = logEmin + width
                    command = '$IR4 {0} do_n_sig --conf-dirs={1} ' \
                            ' --n-trials={2}' \
                            ' --test-zenith={3:07.3f}' \
                            ' --test-ext=0.' \
                            ' --src-ext=0. ' \
                            ' --src-gamma=2. ' \
                            ' --src-thresh={4} ' \
                            ' --src-cutoff={5} ' \
                            ' --sigma={6:.0f}' \
                            ' --beta={7:03.1f}' \
                            ' --seed={8}'.format (
                                this_script,
                                confs,
                                self.opts.n_trials,
                                zenith_deg,
                                logEmin,
                                logEmax,
                                self.opts.sigma,
                                self.opts.beta,
                                i_job)
                    label = 'do_n_sig_diff__' \
                            'zen_{0:07.3f}__' \
                            'src_thresh_{1}__' \
                            'src_cutoff_{2}__' \
                            'seed_{3:08d}'.format (
                                    zenith_deg,
                                    logEmin,
                                    logEmax,
                                    i_job)
                    commands.append (command)
                    labels.append (label)
                    i_job += 1

        submitter.submit_condor00 (
                commands, labels,
                blacklist=self.blacklist,
                max_per_interval=20)
        pjobdir (job_dir)

    @command
    def do_n_sig (self):
        """Calculate n_sig for given n-sigma in beta fraction of trials."""
        # get parameters from command line
        prush ('Getting parameters from commandline...')
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        dec_deg = self.opts.test_dec
        dec = dec_deg / 180.*pi
        test_ext_deg = self.opts.test_ext
        test_ext = test_ext_deg / 180.*pi
        src_ext_deg = self.opts.src_ext
        src_ext = src_ext_deg / 180.*pi
        gamma = self.opts.src_gamma
        thresh = self.opts.src_thresh
        cutoff = self.opts.src_cutoff
        sigma = self.opts.sigma
        beta = self.opts.beta
        n_max = 7 if sigma < 3 else 12
        n_sig_guess = 15 if sigma < 3 else 45
        use_fit = sigma >= 3
        eps = 1e-5 if use_fit else 0

        # get the ts threshold
        if self.opts.oneyear:
          self.analysis86.ps_sindec_width = .2 if dec_deg < 60 else .05
          bg_tsd = bk.get_best (self.bg_tsds86, test_ext, dec_deg)
        else:
          self.analysis.ps_sindec_width = .2 if dec_deg < 60 else .03 #0.05
          bg_tsd = bk.get_best (self.bg_tsds, test_ext, dec_deg)

        ts = bg_tsd.isf (stats.norm.sf (sigma) + eps, fit=use_fit)

        prush ('Setting up injection...')
        if self.opts.oneyear:
          data = hyp.DataHypothesis (self.analysis86)
          gammas = self.analysis86.pdfs_energy_sig['IC86_nopull'].gammas
          ps = hyp.PointSourceHypothesis (self.analysis86, dec, 0, gamma,
                                        energy_range=(10**thresh, 10**cutoff))
        else:
          data = hyp.DataHypothesis (self.analysis)
          gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
          ps = hyp.PointSourceHypothesis (self.analysis, dec, 0, gamma,
                                        energy_range=(10**thresh, 10**cutoff))
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        prush ('Finding n_inj...')
        prush ('- ts > {}'.format (ts))
        prush ('- beta = {0:.3f}'.format (beta))
        prush ('- dec = {0:.3f}'.format (dec_deg))
        prush ('- test_ext = {0:.3f}'.format (test_ext_deg))
        prush ('- src_ext = {0:.3f}'.format (src_ext_deg))
        prush ('- gamma = {0:.3f}'.format (gamma))
        prush ('- thresh = {0:.3f}'.format (thresh))
        prush ('- cutoff = {0:.3f}'.format (cutoff))
        prush ('- ns_bounds = {0:.3f}, {1:.3f}'.format (*self.ns_bounds))
        prush ('- n_sig_guess = {}'.format (n_sig_guess))

        sens = tr.get_sens (ts, beta, n_batch=200, tol=.02)
        flux100TeV = ps.to_flux (sens['n_sig_best'], 100, 1e3) # E0^2 Phi(E0) in TeV/cm^2/s
        fluxGeV = ps.to_flux (sens['n_sig_best']) # E0^2 Phi(E0) in TeV/cm^2/s
        prush ('Sensitivity flux is E0^2 Phi(E0) '
               '= {:.3e} TeV/cm^2/s'.format (flux100TeV))
        prush ('Sensitivity flux is Phi(1GeV) '
               '= {:.3e}'.format (fluxGeV))

        outname = 'sens'

        sm = bk.SavingModel (
            'sigma/beta/test_ext/src_ext/gamma/thresh/cutoff/dec/sens',
            '{:1.0f}/{:3.1f}/{:04.1f}/{:04.1f}/{:4.2f}/{:.2f}/{:.2f}/'
            '{:+07.2f}/{}.pickle'
        )
        if self.opts.oneyear:
            filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/sig_int/'.format (self.mode_dir86),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname)
        else:
            filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/sig_int/'.format (self.mode_dir),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname)
        prush ('->', filename)

    @command
    def do_n_sig_multi (self):
        """Calculate n_sig for doublesource test for given n-sigma in beta fraction of trials."""
        # get parameters from command line
        prush ('Getting parameters from commandline...')
        seed = self.opts.seed
        n_multi = int(self.opts.n_multi)
        n_trials = int (self.opts.n_trials)
        dec_deg = self.opts.test_dec
        dec = dec_deg / 180.*pi
        multi_dec = [dec for i in range(n_multi)]
        test_ext_deg = self.opts.test_ext
        test_ext = test_ext_deg / 180.*pi
        src_ext_deg = self.opts.src_ext
        src_ext = src_ext_deg / 180.*pi
        gamma = self.opts.src_gamma
        thresh = self.opts.src_thresh
        cutoff = self.opts.src_cutoff
        sigma = self.opts.sigma
        beta = self.opts.beta
        n_max = 7 if sigma < 3 else 12
        n_sig_guess = 15 if sigma < 3 else 45

        self.analysis.ps_sindec_width = .2 if dec_deg < 60 else .05

        # get the ts threshold - have to reference the right directory
        #bg_tsd = bk.get_best (self.bg_tsds, test_ext, dec_deg)
        if self.opts.diff_ra:
          bg_tsd = bk.get_best (loading ('{0}/tests/multi_diff/{1}/bg_tsds.dict'.format (self.mode_dir,n_multi)), test_ext, dec_deg)
          multi_ra = np.linspace(0,2*pi,n_multi+1)[:-1] 
          test = 'multi_diff'
        else:
          bg_tsd = bk.get_best (loading ('{0}/tests/multi_same/{1}/bg_tsds.dict'.format (self.mode_dir,n_multi)), test_ext, dec_deg)
          multi_ra = np.array([0 for i in range(n_multi)])
          test = 'multi_same'
        use_fit = sigma >= 3
        eps = 1e-5 if use_fit else 0
        ts = bg_tsd.isf (stats.norm.sf (sigma) + eps, fit=use_fit)

        prush ('Setting up injection...')
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, multi_dec, multi_ra, gamma,
                                        energy_range=(10**thresh, 10**cutoff))
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        prush ('Finding n_inj for multi test...')
        prush ('- ts > {}'.format (ts))
        prush ('- beta = {0:.3f}'.format (beta))
        prush ('- dec = {0:.3f}'.format (dec_deg))
        prush ('- test_ext = {0:.3f}'.format (test_ext_deg))
        prush ('- src_ext = {0:.3f}'.format (src_ext_deg))
        prush ('- gamma = {0:.3f}'.format (gamma))
        prush ('- thresh = {0:.3f}'.format (thresh))
        prush ('- cutoff = {0:.3f}'.format (cutoff))
        prush ('- ns_bounds = {0:.3f}, {1:.3f}'.format (*self.ns_bounds))
        prush ('- n_sig_guess = {}'.format (n_sig_guess))

        sens = tr.get_sens (ts, beta, n_sig_guess=n_sig_guess, n_batch=200, tol=.02)
        flux100TeV = ps.to_flux (sens['n_sig_best'], 100, 1e3) # E0^2 Phi(E0) in TeV/cm^2/s
        fluxGeV = ps.to_flux (sens['n_sig_best']) # E0^2 Phi(E0) in TeV/cm^2/s
        prush ('Sensitivity flux is E0^2 Phi(E0) '
               '= {:.3e} TeV/cm^2/s'.format (flux100TeV))
        prush ('Sensitivity flux is Phi(1GeV) '
               '= {:.3e}'.format (fluxGeV))

        outname = 'sens'

        sm = bk.SavingModel (
            'sigma/beta/test_ext/src_ext/gamma/thresh/cutoff/dec/sens',
            '{:1.0f}/{:3.1f}/{:04.1f}/{:04.1f}/{:4.2f}/{:.2f}/{:.2f}/'
            '{:+07.2f}/{}.pickle'
        )

        filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/tests/{1}/{2}/sig_int/'.format (self.mode_dir,test,n_multi),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname)
        prush ('->', filename)

    @command
    def do_n_sig_northsouth (self):
        """Calculate n_sig for doublesource test for given n-sigma in beta fraction of trials."""
        # get parameters from command line
        prush ('Getting parameters from commandline...')
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        dec_deg = [-30.,30.]
        dec = [d / 180.*pi for d in dec_deg]
        test_ext_deg = self.opts.test_ext
        test_ext = test_ext_deg / 180.*pi
        src_ext_deg = self.opts.src_ext
        src_ext = src_ext_deg / 180.*pi
        gamma = self.opts.src_gamma
        thresh = self.opts.src_thresh
        cutoff = self.opts.src_cutoff
        sigma = self.opts.sigma
        beta = self.opts.beta
        n_max = 7 if sigma < 3 else 12
        n_sig_guess = 15 if sigma < 3 else 45

        self.analysis.ps_sindec_width = .2 if dec_deg < 60 else .05

        # get the ts threshold - have to reference the right directory
        #bg_tsd = bk.get_best (self.bg_tsds, test_ext, dec_deg)
        bg_tsd = loading ('{0}/tests/northsouth/bg_tsds.dict'.format (self.mode_dir))
        use_fit = sigma >= 3
        eps = 1e-5 if use_fit else 0
        ts = bg_tsd.isf (stats.norm.sf (sigma) + eps, fit=use_fit)

        prush ('Setting up injection...')
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, [0,0], gamma,
                                        energy_range=(10**thresh, 10**cutoff))
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        prush ('Finding n_inj for doubletest...')
        prush ('- ts > {}'.format (ts))
        prush ('- beta = {0:.3f}'.format (beta))
        prush ('- dec = {0}'.format (dec_deg))
        prush ('- test_ext = {0}'.format (test_ext_deg))
        prush ('- src_ext = {0}'.format (src_ext_deg))
        prush ('- gamma = {0:.3f}'.format (gamma))
        prush ('- thresh = {0:.3f}'.format (thresh))
        prush ('- cutoff = {0:.3f}'.format (cutoff))
        prush ('- ns_bounds = {0:.3f}, {1:.3f}'.format (*self.ns_bounds))
        prush ('- n_sig_guess = {}'.format (n_sig_guess))

        sens = tr.get_sens (ts, beta, n_sig_guess=n_sig_guess, n_batch=200, tol=.02)
        flux100TeV = ps.to_flux (sens['n_sig_best'], 100, 1e3) # E0^2 Phi(E0) in TeV/cm^2/s
        fluxGeV = ps.to_flux (sens['n_sig_best']) # E0^2 Phi(E0) in TeV/cm^2/s
        prush ('Sensitivity flux is E0^2 Phi(E0) '
               '= {:.3e} TeV/cm^2/s'.format (flux100TeV))
        prush ('Sensitivity flux is Phi(1GeV) '
               '= {:.3e}'.format (fluxGeV))

        outname = 'sens'

        sm = bk.SavingModel (
            'sigma/beta/test_ext/src_ext/gamma/thresh/cutoff/dec/sens',
            '{:1.0f}/{:3.1f}/{}.pickle'
        )

        filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/tests/northsouth/sig_int/'.format (self.mode_dir),
                sigma, beta, outname)
        prush ('->', filename)

    @command
    def collect_sig_int (self):
        """Collect signal injection results and cache."""
        prush ('Collecting signal results...')
        if self.opts.oneyear:
          d1 = bk.get_all (
                '{0}/sig_int'.format (self.mode_dir86),
                'sens.pickle',
                lambda x: x[0])
          saving (d1, '{0}/sig_int.dict'.format (self.mode_dir86))

        else:
          d1 = bk.get_all (
                '{0}/sig_int'.format (self.mode_dir),
                'sens.pickle',
                lambda x: x[0])
          saving (d1, '{0}/sig_int.dict'.format (self.mode_dir))

        # prush ('Collecting sys_spread results...')
        # d = bk.get_all (
        #         '{0}/sig_int'.format (self.mode_dir),
        #         '*sens_sys_spread.pickle',
        #         lambda x: x[0])
        # saving (d, '{0}/sig_int_sys_spread.dict'.format (self.mode_dir))

        # prush ('Collecting sys_spread_late results...')
        # d = bk.get_all (
        #         '{0}/sig_int'.format (self.mode_dir),
        #         '*sens_sys_spread_late.pickle',
        #         lambda x: x[0])
        # saving (d, '{0}/sig_int_sys_spread_late.dict'.format (self.mode_dir))

        # prush ('Collecting sys_domeff results...')
        # d = bk.get_all (
        #         '{0}/sig_int'.format (self.mode_dir),
        #         '*sens_sys_domeff.pickle',
        #         lambda x: x[0])
        # saving (d, '{0}/sig_int_sys_domeff.dict'.format (self.mode_dir))

        # prush ('Collecting sys_domeff_noq results...')
        # d = bk.get_all (
        #         '{0}/sig_int'.format (self.mode_dir),
        #         '*sens_sys_domeff_noq.pickle',
        #         lambda x: x[0])
        # saving (d, '{0}/sig_int_sys_domeff_noq.dict'.format (self.mode_dir))
        return d1

    @command
    def collect_sig_int_multi_same (self):
        """Collect signal injection results for double test and cache."""
        prush ('Collecting signal results...')
        d1 = bk.get_all (
                '{0}/tests/multi_same/{1}/sig_int'.format (self.mode_dir,self.opts.n_multi),
                'sens.pickle',
                lambda x: x[0])
        saving (d1, '{0}/tests/multi_same/{1}/sig_int.dict'.format (self.mode_dir,self.opts.n_multi))
        return d1

    @command
    def collect_sig_int_multi_diff (self):
        """Collect signal injection results for double_opp test and cache."""
        prush ('Collecting signal results...')
        d1 = bk.get_all (
                '{0}/tests/multi_diff/{1}/sig_int'.format (self.mode_dir,self.opts.n_multi),
                'sens.pickle',
                lambda x: x[0])
        saving (d1, '{0}/tests/multi_diff/{1}/sig_int.dict'.format (self.mode_dir,self.opts.n_multi))
        return d1

    @command
    def submit_do_bg_skies (self):
        """Submit bg sky jobs to cluster."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_bg_skies/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        for i_job in xrange (int (self.opts.n_jobs)):
            command = '$IR4 {0} do_bg_skies --conf-dirs={1} ' \
                    ' --n-trials={2}' \
                    ' --nside={3}' \
                    ' --seed={4}'.format (
                        this_script,
                        confs,
                        self.opts.n_trials,
                        self.opts.nside,
                        i_job)
            label = 'do_bg_skies__seed_{0:08d}'.format (i_job)
            commands.append (command)
            labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def do_bg_skies (self):
        """Generate background only full-sky TSDist."""
        ana = self.ana
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        np.random.seed (seed)

        bg_tsds = self.bg_tsds[0]

        nside = self.opts.nside
        zeniths, azimuths = np.array ([
            healpy.pix2ang (nside, i_pix)
            for i_pix in xrange (healpy.nside2npix(nside))]).T

        gtss = []
        prush ('Doing trials...')
        for i in xrange (n_trials):
            e = ens.Ensemble (ana.pdfs, ana.bg)
            #map_ps = e.pretrial_p (bg_tsds, zeniths, azimuths)
            #print (map_ps)
            #gtss.append (-np.log10 (np.min (map_ps)))
            gtss.append (
                    -np.log10 (e.pretrial_p_best (bg_tsds, zeniths, azimuths)))
            if i % 10 == 0:
                prush ('... {0}%'.format (100. * i / n_trials))
        prush ('Done.')
        gtss = np.array (gtss)

        sm = bk.SavingModel (
                'nside/gts_array',
                '{0:04d}/gts_{1:08d}.array',
                )

        filename = sm.save (
                gtss, '{0}/bg_skies/'.format (self.mode_dir), nside, seed)
        prush ('->', filename)

    @command
    def collect_bg_gts (self):
        """
        Collect bg-only -log10(pretrial_p_best) distribution.
        """
        gtsd = ens.TSDist (np.concatenate (map (
            cache.load,
            sorted (glob ('{0}/bg_skies/{1:04d}/gts_*array'.format (
                self.mode_dir, self.opts.nside
            )))
        )))
        saving (gtsd, '{0}/bg_gts.dist'.format (self.mode_dir))
        return gtsd


    @command
    def submit_do_bg_cats (self):
        """Submit bg sky jobs to cluster."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_bg_cats/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        for i_job in xrange (int (self.opts.n_jobs)):
            command = '$IR4 {0} do_bg_cats --conf-dirs={1} ' \
                    ' --n-trials={2}' \
                    ' --nside={3}' \
                    ' --seed={4}'.format (
                        this_script,
                        confs,
                        self.opts.n_trials,
                        self.opts.nside,
                        i_job)
            label = 'do_bg_cats__seed_{0:08d}'.format (i_job)
            commands.append (command)
            labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def do_bg_cats (self):
        """Generate background only full-sky TSDist."""
        ana = self.ana
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        np.random.seed (seed)

        bg_tsds = self.bg_tsds[0]

        nside = self.opts.nside

        catalog = cat.catalog_from_file (
            '{0}/source_list_full.txt'.format (self.root_dir))

        zeniths = catalog.zeniths
        azimuths = catalog.azimuths

        gtss = []
        prush ('Doing trials...')
        for i in xrange (n_trials):
            e = ens.Ensemble (ana.pdfs, ana.bg)
            gtss.append (
                    -np.log10 (e.pretrial_p_best (bg_tsds, zeniths, azimuths)))
            if i % 10 == 0:
                prush ('... {0}%'.format (100. * i / n_trials))
        prush ('Done.')
        gtss = np.array (gtss)

        sm = bk.SavingModel (
                'nside/gts_array',
                'gts_{0:08d}.array',
                )

        filename = sm.save (
                gtss, '{0}/bg_cats/'.format (self.mode_dir), seed)
        prush ('->', filename)

    @command
    def collect_bg_cat_gts (self):
        """
        Collect bg-only -log10(pretrial_p_best) distribution.
        """
        gtsd = ens.TSDist (np.concatenate (map (
            cache.load,
            sorted (glob ('{0}/bg_cats/gts_*array'.format (
                self.mode_dir
            )))
        )))
        saving (gtsd, '{0}/bg_cat_gts.dist'.format (self.mode_dir))
        return gtsd

    # plotting
    
    # PAPER
    @command
    def sens_vs_track (self):
        """Plot integrated sensitivity."""

        sig_int = self.sig_int

        #spectra = [
        #        (2, 5), (2, 6), (2, np.inf),
        #        (2.23, 6), (2.23, np.inf),
        #        (2.46, np.inf),
        #        (2.69, np.inf),
        #        #(2.92, np.inf)
        #        ]

        misc.tex_mpl_rc (False)

        #fig = getfig (aspect=16/10., width=6)
        #fig = getfig (aspect=1., width=5)
        fig = plt.figure(100)
        fig.clf()
        gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
        #ax = fig.add_subplot (111)
        ax = plt.subplot (gs[0])
        rax = plt.subplot (gs[1], sharex=ax)
        colors = soft_colors
        lw = 2
        cutoff = np.inf
        if cutoff == 5:
            ls = '--'
        else:
            ls = '-'
        #sens - gamma = 2
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label2 = 'PS / csky / $E^{{-{0}}}${1}'.format (
                #int (self.conf.root.n_year), gamma, cutoff_str)
                2, cutoff_str)
        sensdir = self.mode_dir+'/sig_int/'
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        x = sig_int[0][.9][0][0][2][2][np.inf]
        sindec2 = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        sys = self.opts.sys_total
        sens2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        mask2 = np.r_[False,sens2[1:] / sens2[:-1] > .5]
        ax.semilogy (sindec2[mask2], sens2[mask2],
                label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #sens - gamma = 3
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label3 = 'PS / csky / $E^{{-{0}}}${1}'.format (
                3, cutoff_str)
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        x = sig_int[0][.9][0][0][3][2][np.inf]
        sindec3 = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        sys = self.opts.sys_total
        sens3 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        mask3 = np.r_[False,sens3[1:] / sens3[:-1] > .5]
        ax.semilogy (sindec3[mask3], sens3[mask3],
                label=label3 + ' - sens', lw=lw, ls=ls, color=colors[2])
        #disc - gamma=2
        x = sig_int[5][.5][0][0][2][2][np.inf]
        sindec2_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec2_disc, disc2,
                label=label2 + ' - disc', lw=lw, ls=ls, color=colors[1])
        #disc - gamma=
        x = sig_int[5][.5][0][0][3][2][np.inf]
        sindec3_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc3 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec3_disc, disc3,
                label=label2 + ' - disc', lw=lw, ls=ls, color=colors[3])
        #ax.semilogy (sindec_w, sens_w,
        #        lw=1, alpha=.5, color=color)

        #Stefan's results
        lw = 2
        thing_sens = np.load ('{0}/data/track_sens.npy'.format (self.root_dir))
        thing_disc = np.load ('{0}/data/track_disc.npy'.format (self.root_dir))
        sd, s2, s3 = thing_sens['dec'], thing_sens['2'], thing_sens['3']
        ax.semilogy (np.sin (sd), s2,
                   color=colors[0], lw=lw, ls=':', alpha = 0.5,
                   label=r'PS + MESE / skylab / $E^{-2}$ - sens')
        ax.semilogy (np.sin (sd), s3 / 100,
                   color=colors[2], lw=lw, ls=':', alpha = 0.5,
                   label=r'PS + MESE / skylab / $E^{-3}$ - sens')

        dd, d2, d3 = thing_disc['dec'], thing_disc['2'], thing_disc['3']
        ax.semilogy (np.sin (dd), d2,
                   color=colors[1], lw=lw, ls=':', alpha = 0.5,
                   label=r'PS + MESE / skylab / $E^{-2}$ - disc')
        ax.semilogy (np.sin (dd), d3 / 100,
                   color=colors[3], lw=lw, ls=':', alpha = 0.5,
                   label=r'PS + MESE / skylab / $E^{-3}$ - disc')
        #ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13.5, 10**-7.)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-13, -7, 7))
        #plt.yticks (np.logspace (-13, -7, 7))

        ax.text (.95, 10**-13.4, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-13.4, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)

        #some Ratio plots
        rat2 = np.interp (sd, np.arcsin(sindec2), sens2) / s2 
        rat3 = np.interp (sd, np.arcsin(sindec3), sens3) / (s3 / 100)
        rat2disc = np.interp (dd, np.arcsin(sindec2_disc), disc2) / d2
        rat3disc = np.interp (dd, np.arcsin(sindec3_disc), disc3) / (d3 / 100)
        rax.plot (np.sin(sd), rat2, color=colors[0])
        rax.plot (np.sin(sd), rat3, color=colors[2])
        rax.plot (np.sin(dd), rat2disc, color=colors[1])
        rax.plot (np.sin(dd), rat3disc, color=colors[3])
        
        rax.set_xlabel (r'$\sin(\delta)$')
        rax.set_ylabel ('Drexel / Stefan')
        rax.set_xlim (-1, 1)
        rax.set_ylim (0.5,1.5)

        ax.grid ()
        rax.grid ()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'my7yr_vs_track_sensitivity_sq')

    @command
    def sens_tests (self):
        """Plot allsky integrated sensitivity and multitests."""

        sig_int = self.sig_int

        misc.tex_mpl_rc (False)
        fig = plt.figure(100)
        #fig = getfig (aspect=1., width=5)

        gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
        #ax = fig.add_subplot (111)
        ax = plt.subplot (gs[0])
        rax = plt.subplot (gs[1], sharex=ax)
        colors = soft_colors
        lw = 2
        cutoff = np.inf
        ls = '-'
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        #7yr sens
        label = 'PS allsky sens - $E^{{-{0}}}${1}'.format (
                #int (self.conf.root.n_year), gamma, cutoff_str)
                2, cutoff_str)
        sensdir = self.mode_dir+'/sig_int/'
        sys = self.opts.sys_total
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        x = sig_int[0][.9][0][0][2][2][np.inf]
        sindec = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        sens = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec, sens,
                label=label, lw=lw, ls=ls, color=colors[0])
        #same / diff tests
        multis = [2,3,4]
        for n in multis:
          label_same = '{} sources, same ra'.format (n)
          label_diff = '{} sources, diff ra'.format (n)
          sensdir_same = self.mode_dir+'/tests/multi_same/{}/sig_int/'.format(n)
          sensdir_diff = self.mode_dir+'/tests/multi_diff/{}/sig_int/'.format(n)
          sig_int_same = bk.get_all(
                sensdir_same,'sens.pickle', lambda x:x[0])
          sig_int_diff = bk.get_all(
                sensdir_diff,'sens.pickle', lambda x:x[0])
          x_same = sig_int_same[0][.9][0][0][2][2][np.inf]
          x_diff = sig_int_diff[0][.9][0][0][2][2][np.inf]
          sindec_same = np.array ([
            np.sin (np.radians(k)) for k in sorted (x_same)])
          sindec_diff = np.array ([
            np.sin (np.radians(k)) for k in sorted (x_diff)])
          sens_same = sys*np.array ([
            x_same[k][-1] for k in sorted (x_same)])
          sens_diff = sys*np.array ([
            x_diff[k][-1] for k in sorted (x_diff)])
          ax.scatter (sindec_same, sens_same,
                label=label_same, color=colors[n])
          ax.scatter (sindec_diff, sens_diff,
                label=label_diff, marker = '+',color=colors[n])
          #rat = sens_diff / sens_same 
          rat_same = sens_same /  np.interp (sindec_same, sindec, sens)
          rat_diff = sens_diff / np.interp (sindec_diff, sindec, sens)
          rax.plot (sindec_same, rat_same, color=colors[n])
          rax.plot (sindec_diff, rat_diff, linestyle = '--', color=colors[n])
        #ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13, 10**-10.)
        ax.set_xlim (-1., 1.)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-13, -10, 4))
        #plt.yticks (np.logspace (-13, -10, 4))

        ax.text (.95, 10**-13.4, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-13.4, 'South', size='small', ha='left', va='bottom')
        #ax.text (.95, 2*10**-13.0, 'North', size='small', ha='right', va='bottom')
        #ax.text (-.95, 2*10**-13.0, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)
        #some Ratio plots
        
        rax.set_xlabel (r'$\sin(\delta)$')
        rax.set_ylabel ('Csky / Skylab')
        rax.set_xlim (-1, 1)
        rax.set_ylim (0.5,1.5)


        ax.grid ()
        rax.grid()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'my7yr_stacking_tests')


    @command
    def sens_oneyear (self):
        """Plot IC86 integrated sensitivity."""

        sig_int = self.sig_int86

        misc.tex_mpl_rc (False)

        fig = plt.figure(100)
        fig.clf()
        gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
        ax = plt.subplot (gs[0])
        rax = plt.subplot (gs[1], sharex=ax)
        colors = soft_colors
        lw = 2
        cutoff = np.inf
        if cutoff == 5:
            ls = '--'
        else:
            ls = '-'
        #sens - gamma = 2
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label2 = 'PS / csky / $E^{{-{0}}}${1}'.format (
                #int (self.conf.root.n_year), gamma, cutoff_str)
                2, cutoff_str)
        sensdir = self.mode_dir+'/oneyear/sig_int/'
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        x = sig_int[0][.9][0][0][2][2][np.inf]
        sindec2 = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        sys = self.opts.sys_total
        sens2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec2, sens2,
                label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #disc - gamma=2
        x = sig_int[5][.5][0][0][2][2][np.inf]
        sindec2_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec2_disc, disc2,
                label=label2 + ' - disc', lw=lw, ls=ls, color=colors[1])

        #Stefan's results
        thing_sens = np.genfromtxt('{0}/data/sens_86.csv'.format (self.root_dir), delimiter = ',')
        thing_disc = np.genfromtxt('{0}/data/disc_86.csv'.format (self.root_dir), delimiter = ',')
        lw = 2
        sd, s2 = thing_sens.T
        dd, d2 = thing_disc.T
        ax.semilogy (sd, s2,
                   color=colors[0], lw=lw, ls=':', alpha = 0.5,
                   label=r'IC86  / skylab / $E^{-2}$ - sens')

        ax.semilogy (dd, d2,
                   color=colors[1], lw=lw, ls=':', alpha = 0.5,
                   label=r'IC86 / skylab / $E^{-2}$ - disc')
        #ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13.5, 10**-7.)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-13, -7, 7))
        #plt.yticks (np.logspace (-13, -7, 7))

        ax.text (.95, 10**-13.4, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-13.4, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)

        #some Ratio plots
        rat2 = np.interp (sd, sindec2, sens2) / s2 
        rat2disc = np.interp (dd, sindec2_disc, disc2) / d2
        rax.plot (sd, rat2, color=colors[0])
        rax.plot (dd, rat2disc, color=colors[1])
        
        rax.set_xlabel (r'$\sin(\delta)$')
        rax.set_ylabel ('Drexel / Stefan')
        rax.set_xlim (-1, 1)
        rax.set_ylim (0.5,1.5)

        ax.grid ()
        rax.grid ()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'oneyear_vs_track_sensitivity')

    @command
    def disc_vs_track (self):
        """Plot integrated sensitivity."""

        sig_int = self.sig_int

        #spectra = [
        #        (2, 5), (2, 6), (2, np.inf),
        #        (2.23, 6), (2.23, np.inf),
        #        (2.46, np.inf),
        #        (2.69, np.inf),
        #        #(2.92, np.inf)
        #        ]
        spectra = [(2, np.inf), (2, 5),
                   #(2.5, np.inf),
                   (3, np.inf)
                  ]

        misc.tex_mpl_rc (False)

        #fig = getfig (aspect=16/10., width=6)
        fig = getfig (aspect=1., width=5)
        ax = fig.add_subplot (111)
        colors = soft_colors
        lw = 2
        for (gamma, cutoff) in spectra:
            #if cutoff == np.inf:
            #    lw, ls = 2, '-'
            #elif cutoff == 6:
            #    lw, ls = 1, '--'
            #else:
            #    lw, ls = 1, ':'
            #colors = {2.: 'b', 2.23: 'g', 2.46: 'orange', 2.69: 'r'}
            if cutoff == 5:
                ls = '--'
            else:
                ls = '-'
            if gamma == 2:
                color = colors[0]
            elif gamma == 2.5:
                color = colors[2]
            elif gamma == 3:
                color = colors[1]
            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at $10^{{{0:.0f}}}$ GeV'.format (cutoff)
            label = '{0} Year Cascades, $E^{{-{1}}}${2}'.format (
                    int (self.conf.root.n_year), gamma, cutoff_str)
            x = bk.get_best (
                    sig_int, 5, 0.5,
                    self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            sindec = np.array ([
                np.sin (k/180.*pi) for k in sorted (x)])
            sys = self.opts.sys_total
            sens = sys * np.array ([
                x[k][-1] for k in sorted (x)])
            if cutoff == 5:
                mask = sindec != np.sin (np.radians (27))
                ax.semilogy (sindec[mask], sens[mask],
                        label=label, lw=lw, ls=ls, color=color)
            else:
                ax.semilogy (sindec, sens,
                        label=label, lw=lw, ls=ls, color=color)
            #ax.semilogy (sindec_w, sens_w,
            #        lw=1, alpha=.5, color=color)

        lw = 2
        thing = np.load ('{0}/data/track_disc.npy'.format (self.root_dir))
        d, s2, s3 = thing['dec'], thing['2'], thing['3']
        ax.semilogy (np.sin (d), s2,
                     color=colors[0], lw=lw, ls=':',
                     label=r'7 Year Tracks, $E^{-2}$')
        ax.semilogy (np.sin (d), s3 / 100,
                     color=colors[1], lw=lw, ls=':',
                     label=r'7 Year Tracks, $E^{-3}$')
        ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13.5, 10**-8)
        ymin, ymax = ax.get_ylim ()
        plt.yticks (np.logspace (-13, -9, 5))
        ax.text (.95, 10**-13.4, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-13.4, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper center',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)
        ax.grid ()
        plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'cascade_vs_track_disc_sq')

    @command
    def sens_vs_old (self):
        """Plot integrated sensitivity."""

        sig_int = self.sig_int
        sig_int2 = cache.load ('/home/mike/work/i3/scripts/meseps/csky2yr/anaroot/paper/yrs2/sig_int.dict')
        sig_int26 = cache.load ('/home/mike/work/i3/scripts/meseps/csky2yr/anaroot/paper/yrs6/sig_int.dict')

        #spectra = [
        #        (2, 5), (2, 6), (2, np.inf),
        #        (2.23, 6), (2.23, np.inf),
        #        (2.46, np.inf),
        #        (2.69, np.inf),
        #        #(2.92, np.inf)
        #        ]
        spectra = [(2, np.inf),
                   #(2.5, np.inf),
                   (3, np.inf),
                   #(2, 5)
                  ]
        sys = self.opts.sys_total

        lss = ['-', ':', '--', ':']
        lss = ['-', '--']
        lw = 2
        colors = soft_colors
        misc.tex_mpl_rc (False)

        fig = getfig (aspect=1., width=5)
        ax = fig.add_subplot (111)

        spectra = [(2, np.inf),
                   (3, np.inf), ]
        lss = ['-', '--']
        for (gamma, cutoff) in spectra:
            #if cutoff == np.inf:
            #    lw, ls = 2, '-'
            #elif cutoff == 6:
            #    lw, ls = 1, '--'
            #else:
            #    lw, ls = 1, ':'
            #colors = {2.: 'b', 2.23: 'g', 2.46: 'orange', 2.69: 'r'}
            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
            label = '{0} Year Cascades, $E^{{-{1}}}${2}'.format (
                    2, gamma, cutoff_str)
            x2 = bk.get_best (
                    sig_int2, 0, 0, .9,
                    self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            sindec2 = np.array ([
                np.cos (-k/180.*pi) for k in sorted (x2)])
            sys = self.opts.sys_total
            sens2 = sys * np.array ([
                x2[k][-1] for k in sorted (x2)])
            ax.semilogy (sindec2, sens2,
                         label=label, lw=lw,
                         ls='--', color=colors[gamma-2],
                         zorder=5)

        for ((gamma, cutoff), ls) in izip (spectra, lss):
            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
            label = '{0} Year Cascades, $E^{{-{1}}}${2}'.format (
                    int (self.conf.root.n_year), gamma, cutoff_str)
            x = bk.get_best (
                    sig_int, 0, .9,
                    self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            sindec = np.array ([
                np.sin (k/180.*pi) for k in sorted (x)])
            sens = sys * np.array ([
                x[k][-1] for k in sorted (x)])
            ax.semilogy (sindec, sens,
                         label=label, lw=lw,
                         ls='-', color=colors[gamma-2],
                         zorder=10)

            #x = bk.get_best (
            #        sig_int26, 0, 0, .9,
            #        self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            #sindec = np.array ([
            #    np.cos (-k/180.*pi) for k in sorted (x)])
            #sens = sys * np.array ([
            #    x[k][-1] for k in sorted (x)])
            #ax.semilogy (sindec, sens,
            #             lw=1, ls=ls, color=colors[0],
            #             zorder=8)

        ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-12, 10**-9.8)
        ymin, ymax = ax.get_ylim ()
        plt.yticks (np.logspace (-12, -10, 3))
        ax.text (.95, 10**-11.9, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-11.9, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper center',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)
        ax.grid ()
        plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.1 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'cascade_vs_old_sensitivity_sq')

    @command
    def sens_vs_mainz (self):
        """Plot integrated sensitivity."""

        sig_int = self.sig_int
        sig_intm = cache.load (
            '/home/mike/work/i3/scripts/meseps/csky6yr/anaroot/old/'
            'mainz_tev/respe/sig_int.dict')

        #spectra = [
        #        (2, 5), (2, 6), (2, np.inf),
        #        (2.23, 6), (2.23, np.inf),
        #        (2.46, np.inf),
        #        (2.69, np.inf),
        #        #(2.92, np.inf)
        #        ]
        spectra = [(2, np.inf),
                   #(2.5, np.inf),
                   (3, np.inf),
                   #(2, 5)
                  ]
        sys = self.opts.sys_total

        lss = ['-', ':', '--', ':']
        lss = ['-', '--']
        lw = 2
        colors = soft_colors
        rcolors = ['.5', 'k']
        misc.tex_mpl_rc (False)

        fig = getfig (aspect=1., width=5)

        gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
        ax = plt.subplot (gs[0])
        rax = plt.subplot (gs[1], sharex=ax)

        spectra = [(2, np.inf),
                   (3, np.inf), ]
        lss = ['-', '--']

        for (gamma, cutoff) in spectra:
            #if cutoff == np.inf:
            #    lw, ls = 2, '-'
            #elif cutoff == 6:
            #    lw, ls = 1, '--'
            #else:
            #    lw, ls = 1, ':'
            #colors = {2.: 'b', 2.23: 'g', 2.46: 'orange', 2.69: 'r'}
            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
            label = '{0} Year Cascades, $E^{{-{1}}}${2} (Mainz)'.format (
                    6, gamma, cutoff_str)
            xm = bk.get_best (
                    sig_intm, 0, 0, .9,
                    self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            sindecm = np.array ([
                np.cos (-k/180.*pi) for k in sorted (xm, reverse=True)])
            sys = self.opts.sys_total
            sensm = sys * np.array ([
                xm[k][-1] for k in sorted (xm, reverse=True)])
            ax.semilogy (sindecm, sensm,
                         label=label, lw=lw,
                         ls='--', color=colors[gamma-2],
                         zorder=5)

            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
            label = '{0} Year Cascades, $E^{{-{1}}}${2}'.format (
                    int (self.conf.root.n_year), gamma, cutoff_str)
            x = bk.get_best (
                    sig_int, 0, .9,
                    self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            sindec = np.array ([
                np.sin (k/180.*pi) for k in sorted (x)])
            sens = sys * np.array ([
                x[k][-1] for k in sorted (x)])
            ax.semilogy (sindec, sens,
                         label=label, lw=lw,
                         ls='-', color=colors[gamma-2],
                         zorder=10)

            sd = np.linspace (-.99, .99, 100)
            s = np.interp (sd, sindec, sens)
            sm = np.interp (sd, sindecm, sensm)
            #print (s[-4:])
            #print (sm[-4:])
            rax.plot (sd, s / sm, color=colors[gamma-2], lw=2)

            #x = bk.get_best (
            #        sig_int26, 0, 0, .9,
            #        self.opts.test_ext, self.opts.src_ext, gamma, 3, cutoff)
            #sindec = np.array ([
            #    np.cos (-k/180.*pi) for k in sorted (x)])
            #sens = sys * np.array ([
            #    x[k][-1] for k in sorted (x)])
            #ax.semilogy (sindec, sens,
            #             lw=1, ls=ls, color=colors[0],
            #             zorder=8)

        rax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-12.1, 10**-9.8)
        ymin, ymax = ax.get_ylim ()
        plt.sca(ax)
        plt.yticks (np.logspace (-12, -10, 3))
        ax.text (.95, 10**-12, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-12, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper center',
            prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        ax.grid ()
        #leg.get_frame ().set_linewidth (0)

        rax.axhline (1, color='.5', ls='--')
        rax.set_ylim (.6, 1.4)
        rax.set_ylabel ('Current/Mainz')
        rax.grid ()
        plt.sca (rax)
        plt.yticks ([.8, 1, 1.2])

        #plt.tight_layout ()
        fig.subplots_adjust (top=.93, left=.15, right=.95)
        icprelim (ax, x=-1, y=1.1 * ymax, ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        savingfig (fig, plot_dir, 'cascade_vs_mainz_sensitivity_sq')


    @command
    def sens (self):
        """Plot integrated sensitivity."""

        sig_int = self.sig_int
        #spectra = [
        #        (2, 5), (2, 6), (2, np.inf),
        #        (2.23, 6), (2.23, np.inf),
        #        (2.46, np.inf),
        #        (2.69, np.inf),
        #        #(2.92, np.inf)
        #        ]
        max_gamma = 4.
        spectra = [(2, max_gamma)]
                   #(2.5, max_gamma),
                   #(3, max_gamma)]

        misc.tex_mpl_rc ()

        fig = getfig (aspect=16/9., width=6)
        ax = fig.add_subplot (111)
        for (gamma, cutoff) in spectra:
            #if cutoff == np.inf:
            #    lw, ls = 2, '-'
            #elif cutoff == 6:
            #    lw, ls = 1, '--'
            #else:
            #    lw, ls = 1, ':'
            #colors = {2.: 'b', 2.23: 'g', 2.46: 'orange', 2.69: 'r'}
            if gamma == 2:
                lw = 2
                ls = '-'
            elif gamma == 2.5:
                lw = 1
                ls = '-'
            elif gamma == 2.7:
                lw = 1
                ls = '--'
            color = 'k'
            #colors = {2.: 'b', 2.23: 'g', 2.46: 'orange', 2.69: 'r'}
            cutoff_str = '' if cutoff == np.inf \
                    else ' cutoff at $10^{{{0:.0f}}}$ GeV'.format (cutoff)
            #label = 'sensitivity, $E^{{-{0:.2f}}}${1}'.format (
            #        gamma, cutoff_str)
            label = '{0} Year Sensitivity ($E^{{-{1}}}${2})'.format (
                int (self.conf.root.n_year),
                '2' if gamma == 2 else format (np.round (gamma, 1), '.1f'),
                cutoff_str)
            sensdir = self.mode_dir+'/sig_int/'
            sig_int = bk.get_all(
                    sensdir,'sens.pickle', lambda x:x[0])
            x = sig_int[0][.9][0][0][2][2][np.inf]
            sindec = np.array ([
                np.cos (-k/180.*pi) for k in sorted (x)])
            sens = np.array ([
                x[k][-1] for k in sorted (x)])
            ax.semilogy (sindec, sens,
                    label=label, lw=lw, ls=ls, color=color)
        ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13., 10**-10.)
        ax.legend (loc='upper right', prop=propsmaller, handlelength=4)
        ax.grid ()
        fig.subplots_adjust (bottom=.14, top=.91, right=.97)

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        savingfig (fig, plot_dir, 'ps_7yr_sensitivity')


    def get_diffsens (self, zenith_deg, width):
        """Get the differential sensitivity for a zenith and bin width."""

        sig_int = self.sig_int

        #logEmins = np.arange (3. + width, 7, width)
        #logEmins = np.arange (3. + width, 9, width)
        logEmins = np.arange (3. + width, 8.2, width)
        logEmaxs = logEmins + width
        xs = [bk.get_best (
                sig_int, 0, 0, .9,
                0, 0, 2, logEmin, logEmax, zenith_deg)
             for (logEmin, logEmax) in izip (logEmins, logEmaxs)]

        bins = np.r_[logEmins, logEmaxs[-1]]
        values = [x[-1] for x in xs]
        return histlite.Hist (10**bins, values)

    @command
    def diff_cmp (self):
        """Compare differential sensitivities."""

        misc.tex_mpl_rc ()

        # load up old calculations
        sens_dict = loading (
                '{0}/data/other_sens_and_disco_dec_minus60_yrs3.dict'.format (
                    self.root_dir))
        sens_hists = {}
        sens_dict['MEST'] = sens_dict['MESE']
        del sens_dict['MESE']
        thru = 'Throughgoing tracks'
        comb = 'Combined IceCube+ANTARES tracks'
        sens_dict[thru] = sens_dict['IC86']
        del sens_dict['IC86']
        for k in sens_dict:
            bins = np.r_[np.array (sens_dict[k][0]) - .125,
                    sens_dict[k][0][-1] + .25]
            values = sens_dict[k][1]
            sens_hists[k] = histlite.Hist (10**bins, values)

        ## fetch IceCube+ANTARES sens
        ##comb_bins = 10**np.arange (2.25, 8.1, .25)
        ##comb_values = np.loadtxt (
        ##        '{0}/data/comb_sens.txt'.format (self.root_dir)).T[1]
        ##sens_hists[comb] = histlite.Hist (comb_bins, comb_values)
        ant_bins = 10**np.arange (2.25, 8.1, .25)
        ant_values = np.loadtxt (
                '{0}/data/ant_sens.txt'.format (self.root_dir)).T[1]
        sens_hists['ANTARES'] = histlite.Hist (ant_bins, ant_values)

        # fetch MESC values
        sens_hists['MESC'] = self.opts.sys_total * self.get_diffsens (150, .25)

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))

        # plot time
        for (first, suff) in [(2, '_lessbusy'), (0, '')]:
            fig = getfig ()
            ax = fig.add_subplot (111)
            ax.loglog ()
            colors = ('b', 'orange', 'g', 'purple', 'magenta', 'r')
            curves = ('LESE', 'STeVE', 'MEST', thru, 'ANTARES', 'MESC')
            #colors = ('b', 'orange', 'g', 'purple', 'r')
            #curves = ('LESE', 'STeVE', 'MEST', thru, 'MESC')
            for i, k in enumerate (curves[first:]):
                histlite.plot1d (ax, sens_hists[k],
                        label=k,
                        color='r' if k == 'MESC' else colors[i],
                        lw=2 if k == 'MESC' else 1)

            if first:
                ty = 10**-6.5
            else:
                ty = 10**-5.5
            ax.text (10**4.53, ty, r'$\delta=-60^\circ$',
                    size=14, ha='center', va='center')
            ax.set_xlim (1e2, 1e9)
            ax.set_ylim (1e-11, 1e-5)
            ax.grid ()
            ax.set_xlabel (r'$E$\,\,\,[GeV]')
            ax.set_ylabel (r'$E^2 \cdot dN/dE\,\,\,'
                    '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
            ax.legend (loc='upper right', prop=propsmall, handlelength=3)
            savingfig (fig, plot_dir, 'diff_sens' + suff)

        # one more, pub quality
        fig = getfig (aspect=16/9., width=6)
        ax = fig.add_subplot (111)
        ax.loglog ()
        curves = ('LESE', 'STeVE', 'MEST', thru, 'ANTARES', 'MESC')
        histlite.plot1d (ax, sens_hists['ANTARES'],
                         label='ANTARES', color='.5', ls='--')
        histlite.plot1d (ax, sens_hists['MEST'],
                         label='Starting Tracks', color='k' )
        histlite.plot1d (ax, sens_hists[thru],
                         label='Throughgoing Tracks', color='.5')
        histlite.plot1d (ax, sens_hists['MESC'],
                         label='Cascades (New)', color='k', lw=2)
        ax.text (10**4.51, 10**-6.5, r'$\delta=-60^\circ$',
                size=14, ha='center', va='center')
        ax.set_xlim (1e3, 1e9)
        ax.set_ylim (3e-11, 1e-6)
        ax.grid ()
        ax.set_xlabel (r'$E$\,\,\,[GeV]')
        ax.set_ylabel (r'$E^2 \cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.legend (loc='upper right', prop=propsmall, handlelength=3)
        fig.subplots_adjust (bottom=.14, top=.95, right=.97)
        prefix = 'MESC_{0}yr'.format (int (self.conf.root.n_year))
        savingfig (fig, plot_dir, '{0}_diff_sens'.format (prefix))


if __name__ == '__main__':
    self = Csky6yr ()
    self.run ()
    try:
        import pandas as pd
    except:
        pass


hl = histlite


