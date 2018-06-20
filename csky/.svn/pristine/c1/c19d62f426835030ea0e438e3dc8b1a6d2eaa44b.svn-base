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

class teststack (Timed, Driver):

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
                default=2.0, type=float, metavar='GAMMA',
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

        parser.add_option ('--seed', dest='seed',
                default=0, type=int, metavar='SEED',
                help='initialize RNG with SEED')

        parser.add_option ('--just-diff', dest='just_diff',
                default=False, action='store_true',
                help='only do bg trials for diff sens dec(s)')

        parser.add_option ('--blacklist', dest='blacklist',
                default='cobol61,cobol65', help='nodes to avoid')

        parser.add_option ('--cat', dest='cat',
                default=None,
                help='specifies a catalog for stacking')

        parser.add_option ('--weights', dest='weights',
                default=None,
                help='specifies a weighting for stacking')

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
    def psdata (self):
        try:
            return self._psdata
        except:
            self._psdata = loading ('{0}/psdata.vars'.format (self.mode_dir))
            return self._psdata

    @property
    def datasets (self):
        try:
            return self._datasets
        except:
            datasetfolder = self.root_dir+'/data/datasets/'
            dataset40   = loading ('{}IC40.dataset'.format (datasetfolder))
            dataset59   = loading ('{}IC59.dataset'.format (datasetfolder))
            dataset79   = loading ('{}IC79.dataset'.format (datasetfolder))
            dataset86   = loading ('{}IC86.dataset'.format (datasetfolder))
            dataset2012 = loading ('{}IC86II.dataset'.format (datasetfolder))
            dataset2013 = loading ('{}IC86III.dataset'.format (datasetfolder))
            dataset2014 = loading ('{}IC86IV.dataset'.format (datasetfolder))
            datasetMESE = loading ('{}MESE.dataset'.format (datasetfolder))
            if 'yrs1' in self.mode.words:
              self._datasets = [dataset86]
            elif 'yrs4' in self.mode.words:
              self._datasets = [dataset40,dataset59,dataset79sirin,dataset86]
            else:
              self._datasets = [dataset40,dataset59,dataset79,dataset86,dataset2012,dataset2013,dataset2014,datasetMESE]
            return self._datasets

    @property
    def analysis (self):
        try:
            return self._analysis
        except:
            #if self.conf.
            self._analysis = loading ('{0}/ps.analysis'.format (self.mode_dir))
            return self._analysis

    @property
    def ns_bounds (self):
        return (0, .999)


    @property
    def bg_tsds (self):
        if hasattr (self, '_bg_tsds'):
            return self._bg_tsds
        try:
          if self.opts.cat:
            self._bg_tsds = loading ('{0}/cats/{1}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat))
          else:
            self._bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))
        except:
            self._bg_tsds = self.collect_bg_trials ()
        return self._bg_tsds

    @property
    def sig_int (self):
        try:
          if self.opts.cat:
            return loading ('{0}/cats/{1}/sig_int.dict'.format (self.mode_dir,cat))
          else:
            return loading ('{0}/sig_int.dict'.format (self.mode_dir))
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

    # setting up stefan's npz data

    @command
    def build_psdata (self):
        """
        Build stefan's npz data.
        """
        psdata = Vars ()

        #ids_exp = ['{}'.format (s)
        #       for s in 'IC40 IC59 IC79 IC86I epinat_3yr'.split()]
        ids_exp = ['{}'.format (s)
               for s in 'IC40 IC59 IC79b IC86 IC86-2012 IC86-2013 IC86-2014 MESE MESE_followup'.split()]
        filenames_exp = ['{}/stefandata/{}_exp.npy'.format (self.root_dir, i)
                     for i in ids_exp]

        #ids_mc = ['{}'.format (s)
        #       for s in 'IC40 IC59 IC79 IC86I epinat_3yr'.split()]
        ids_mc = ['{}'.format (s)
               for s in 'IC40 IC59 IC79b IC86 IC86-2012 MESE'.split()]
        filenames_mc = ['{}/stefandata/{}_MC.npy'.format (self.root_dir, i)
                     for i in ids_mc]

        prush ('Loading MC and exp datasets...')
        contents_exp = map (np.load, filenames_exp)
        contents_mc = map (np.load, filenames_mc)

        nu40 = psdata.nu40 = Arrays ()
        nu59 = psdata.nu59 = Arrays ()
        nu79 = psdata.nu79 = Arrays ()
        nu79sirin = psdata.nu79sirin = Arrays ()
        nu86 = psdata.nu86 = Arrays ()
        nu2012 = psdata.nu2012 = Arrays ()
        nuMESE = psdata.nuMESE = Arrays ()
        nus = [nu40,nu59,nu79,nu86,nu2012,nuMESE] 

        for nu,mc in zip(nus,contents_mc): 
            nu.true_energy =  mc['trueE']
            nu.true_zenith =  pi/2.+mc['trueDec']
            nu.true_azimuth = mc['trueRa']
            nu.energy =       10**mc['logE']
            nu.zenith =       np.pi/2.+np.arcsin(mc['sinDec'])
            nu.azimuth =      mc['ra']
            nu.oneweight =    mc['ow']
            nu.sigma =        mc['sigma']
            if 'dist' in mc.dtype.names:
              nu.dist =       mc['dist']
            else:
              nu.dist =       np.zeros_like(nu.energy)

        mu40 = psdata.mu40 = Arrays ()
        mu59 = psdata.mu59 = Arrays ()
        mu79 = psdata.mu79 = Arrays ()
        mu86 = psdata.mu86 = Arrays ()
        mu2012 = psdata.mu2012 = Arrays ()
        mu2013 = psdata.mu2013 = Arrays ()
        mu2014 = psdata.mu2014 = Arrays ()
        muMESE_first = Arrays ()
        muMESEfollowup = Arrays ()

        mus = [mu40,mu59,mu79,mu86,mu2012,mu2013,mu2014,muMESE_first,muMESEfollowup] 

        for mu,exp in zip(mus,contents_exp): 
            mu.energy =   10**exp['logE']
            mu.zenith =   np.pi/2.+np.arcsin(exp['sinDec'])
            mu.azimuth =  exp['ra']
            mu.sigma =  exp['sigma']
            if 'dist' in exp.dtype.names:
              mu.dist =   exp['dist']
            else:
              mu.dist =   np.zeros_like(mu.energy)

        #also need to add the followup to the exp 
        muMESE = psdata.muMESE = combined_arrays((muMESE_first,muMESEfollowup))

        mus = [mu40,mu59,mu79,mu86,mu2012,mu2013,mu2014,muMESE] 

        for d in nus+mus:
            d.dec = d.zenith - pi/2
            d.ra = d.azimuth
            if 'true_zenith' in d:
                d.true_dec = d.true_zenith - pi/2
                d.true_ra = d.true_azimuth
        saving (psdata, '{0}/psdata.vars'.format (self.mode_dir))


    @command
    def setup_datasets (self):
        """Create Dataset."""
        nu40 = self.psdata.nu40
        nu59 = self.psdata.nu59
        nu79 = self.psdata.nu79
        nu86 = self.psdata.nu86
        nu2012 = self.psdata.nu2012
        nuMESE = self.psdata.nuMESE
        mu40 = self.psdata.mu40
        mu59 = self.psdata.mu59
        mu79 = self.psdata.mu79
        mu86 = self.psdata.mu86
        mu2012 = self.psdata.mu2012
        mu2013 = self.psdata.mu2013
        mu2014 = self.psdata.mu2014
        muMESE = self.psdata.muMESE
        livetime40  = 375.539 * 86400
        livetime59  = 348.138 * 86400
        livetime79  = 315.506 * 86400
        livetime86 = 332.61 *  86400
        livetime86II = 330.38 * 86400
        livetime86III = 359.95 * 86400
        livetime86IV = 367.21 * 86400
        livetimeMESE = 1715 * 86400

        dataset40 = ana.Dataset ('IC40', livetime40, sig=nu40, data=mu40)
        dataset59 = ana.Dataset ('IC59', livetime59, sig=nu59, data=mu59)
        dataset79 = ana.Dataset ('IC79', livetime79, sig=nu79, data=mu79)
        dataset86 = ana.Dataset ('IC86', livetime86, sig=nu86, data=mu86)
        dataset86II = ana.Dataset ('IC86II', livetime86II, sig=nu2012, data=mu2012)
        dataset86III = ana.Dataset ('IC86III', livetime86III, sig=nu2012, data=mu2013)
        dataset86IV = ana.Dataset ('IC86IV', livetime86IV, sig=nu2012, data=mu2014)
        datasetMESE = ana.Dataset ('MESE', livetimeMESE, sig=nuMESE, data=muMESE)
        datasetfolder = self.root_dir+'/data/datasets/'
        saving (dataset40, '{}IC40.dataset'.format (datasetfolder))
        saving (dataset59, '{}IC59.dataset'.format (datasetfolder))
        saving (dataset79, '{}IC79.dataset'.format (datasetfolder))
        saving (dataset86, '{}IC86.dataset'.format (datasetfolder))
        saving (dataset86II, '{}IC86II.dataset'.format (datasetfolder))
        saving (dataset86III, '{}IC86III.dataset'.format (datasetfolder))
        saving (dataset86IV, '{}IC86IV.dataset'.format (datasetfolder))
        saving (datasetMESE, '{}MESE.dataset'.format (datasetfolder))

    @command
    def setup_ana (self):
        """
        Set up main analysis objects.
        """

        n_year = self.conf.root.n_year
        print ("Making analysis object for " + str(n_year) + " years of data")
        #nu = self.cdata.nu
        #atm = self.cdata.atm
        #diffuse = self.cdata.diffuse
        #mu = self.cdata.mu


        datasets = self.datasets
        if 'yrs4' in self.mode.words:
          print("4yr")
          analysis = ana.Analysis (
            datasets, 
            dec_kw=dict (
                bins=101, range=(-1, 1),
                keys = ('IC40','IC59','IC79','IC86'),
            ),
            energy_kw=dict (
                sindec_bins=101, logenergy_bins=36,
                logenergy_range=(1, 10),
                keys = ('IC40','IC59','IC79','IC86')
            )
          )
        elif 'yrs1' in self.mode.words:
          print("one year")
          analysis = ana.Analysis (
            datasets, 
            dec_kw=dict (
                bins=101, range=(-1, 1),
                keys = (['IC86']),
            ),
            energy_kw=dict (
                sindec_bins=101, logenergy_bins=36,
                logenergy_range=(1, 10),
                #smooth=4, smooth_bins=[300, 300, 1],
                keys = (['IC86'])
            )
          )
        else:
          analysis = ana.Analysis (
            datasets, 
            dec_kw=dict (
                bins=101, range=(-1, 1),
                keys = ('IC40','IC59','IC79','IC86','IC86II','IC86III','IC86IV'),
            ),
            energy_kw=dict (
                #sindec_bins=41, logenergy_bins=36,
                sindec_bins=101, logenergy_bins=36,
                logenergy_range=(1, 10),
                #smooth=4, smooth_bins=[300, 300, 1],
                keys = ('IC40','IC59','IC79','IC86','IC86II','IC86III','IC86IV')
            )
          )
          analysis.set_dec_pdf(bins=101, range=(-1, 0),
                keys = ['MESE'],
                ) 
          analysis.set_energy_pdfs (
            sindec_bins=1, logenergy_bins=72, dist_bins=50,
            sindec_range=(-1, -0.087), logenergy_range=(1, 10),
            pre_smooth=.5, smooth=4, smooth_bins=[1, 300, 300], keys=['MESE']
          )
        saving (analysis, '{}/ps.analysis'.format (self.mode_dir))

    @command
    def sindecpdf (self):
        """
        Make '1D' sindec pdf plots.
        """
        analysis = self.analysis
        for year in analysis.datasets.keys():
            sindecs = np.linspace(-1.,1.,101)
            plot_dir  = misc.ensure_dir ('{0}/plots/sindecpdfs/{1}'.format (self.mode_dir,year))
            fig, ax = pfig()#plt.subplots (num=10)
            ax.cla()
            histlite.plot1d (ax, analysis.pdf_space_bg[year].h)
            ax.set_title("bg space pdf - {}".format(year))
            savingfig (fig, plot_dir, 'bg space pdf - {}'.format(year))  
            plt.close(fig)

    @command
    def energypdf (self):
        """
        Make '2D' energy pdf plots.
        """
        analysis = self.analysis
        for g in np.linspace(1,4,7):#analysis.gammas:
          for year in analysis.datasets.keys():
            plot_dir  = misc.ensure_dir ('{0}/plots/energypdfs/{1}'.format (self.mode_dir,year))
            fig, ax = pfig()#plt.subplots (num=10)
            fig_bg, ax_bg = pfig()
            fig_sig, ax_sig = pfig()
            ax.cla()
            ax_bg.cla()
            ax_sig.cla()
            if year != 'MESE':
              histlite.plot2d (ax, analysis.pdfs_energy_sig[year].epdfs[g].h[:,:,0] / analysis.pdf_energy_bg[year].h[:,:,0],
              cmap='seismic', log=True, vmin=1e-4, vmax=1e4, cbar=True)
            else:
              histlite.plot2d (ax, analysis.pdfs_energy_sig[year].epdfs[g].h[-0.5].T / analysis.pdf_energy_bg[year].h[-0.5].T,
              cmap='RdBu_r', log=True, vmin=1e-2, vmax=1e2,cbar=True)
              ax.set_ylim(2,8.5)
              histlite.plot2d (ax_bg, analysis.pdf_energy_bg[year].h[-0.5].T,
              cmap='OrRd_r', log=True, vmin=10**(-4.5), vmax=1e-2, cbar=True)
              histlite.plot2d (ax_sig, analysis.pdfs_energy_sig[year].epdfs[g].h[-0.5].T,
              cmap='OrRd_r', log=True, vmin=10**(-4.5), vmax=1e-2, cbar=True)
              ax_bg.set_ylim(2,8.5) ; ax_sig.set_ylim(2,8.5)
              xlab = r'starting distance $[\mathrm{m}]$'
              ylab = r'$log(E[\mathrm{GeV}])$'
              ax_bg.set_xlabel(xlab) ; ax_sig.set_xlabel(xlab) ; ax.set_xlabel(xlab)
              ax_bg.set_ylabel(ylab) ; ax_sig.set_ylabel(ylab) ; ax.set_ylabel(ylab)
              ax_bg.set_title('Experimental Data')
              ax_sig.set_title(r'Monte Carlo simulation, $\gamma=${}'.format(g))
              ax.set_title('Smoothed Ratio MC/EXP')
              savingfig (fig_bg, plot_dir, 'E-{}_bg'.format(g))  
              savingfig (fig_sig, plot_dir, 'E-{}_sig'.format(g))  
            savingfig (fig, plot_dir, 'E-{}'.format(g))  
            plt.close(fig)



    # background and signal trials

    @command
    def do_one_TS (self):
        """
        Do a single unscrambled test for TS, fitted ns.
        """

        seed = self.opts.seed
        weight = self.opts.weights
        n_trials = int (self.opts.n_trials)
        n_year = self.conf.root.n_year
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
        if not weight:
          weight = 'equal'
          weights = np.ones_like(params['dec'])
        else:
          weights = params[str(self.opts.weights)]
        if 'extension' in params.keys():
            test_ext = params['extension']
        else:
            test_ext = np.zeros_like(params['dec'])
        dec_deg = np.degrees(params['dec']) 
        ra_deg = np.degrees(params['ra']) 
        test_ext_deg = [ext for ext in test_ext]
        dec = dec_deg / 180.*pi
        ra = ra_deg / 180.*pi

        np.random.seed (seed)
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, 2, extensions = test_ext_deg, weights = weights, sigsub = False)#self.conf.root.sigsub)

        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        c = tr.get_one_ts_ns_gamma (TRUTH=True)

        sm = bk.SavingModel (
                'ts',
                'one_ts.dict',
                )
        filename = sm.save (
                 c, '{0}/cats/{1}/{2}/one_ts_truth'.format (self.mode_dir,self.opts.cat,weight))

        prush ('->', filename)

    @command
    def submit_do_bg_trials (self):
        """Submit bg TSDist stacking jobs to cluster."""
        cat = self.opts.cat
        print ("Submitting trials for catalog".format(cat))
        job_root = self.conf.root.cluster_job_dir
        job_dir = misc.ensure_dir ('{0}/do_bg_trials_stack/{1}'.format (job_root, job_id))
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))

        #Weighting
        weights = self.opts.weights
        if not weights:
           weights = 'equal'

        for i_job in xrange (int (self.opts.n_jobs)):
          command = '{} {} do_bg_trials --conf-dirs={} ' \
                    ' --cat={}' \
                    ' --weights={}' \
                    ' --n-trials={}' \
                    ' --seed={}'.format (
                        env_shell,
                        this_script,
                        confs,
                        cat,
                        weights,
                        self.opts.n_trials,
                        i_job)
          label = 'do_bg_trials_cat_{}_seed_{:08d}'.format (
                    cat, i_job)
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
        n_year = self.conf.root.n_year
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
        if self.opts.weights == 'equal':
          weights = np.ones_like(params['dec'])
        else:
          weights = params[str(self.opts.weights)]
        if 'extension' in params.keys():
            test_ext = params['extension']
        else:
            test_ext = np.zeros_like(params['dec'])
        dec_deg = np.degrees(params['dec']) 
        ra_deg = np.degrees(params['ra']) 
        test_ext_deg = [ext for ext in test_ext]
        dec = dec_deg / 180.*pi
        ra = ra_deg / 180.*pi

        np.random.seed (seed)
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, 2, extensions = test_ext_deg, weights = weights, sigsub = False)#self.conf.root.sigsub)

        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        c = tr.get_Chi2TSD (n_trials, 500)

        if self.opts.cat:
          sm = bk.SavingModel (
                'tsd',
                '{:08d}.chi2',
                )
          filename = sm.save (
                 c, '{0}/cats/{1}/{2}/bg_tsds'.format (self.mode_dir,self.opts.cat,self.opts.weights),
	         seed)
        else:
          sm = bk.SavingModel (
                'test_ext_deg/dec_deg/tsd',
                '{:04.1f}/{:+07.2f}/{:08d}.chi2',
                )
          filename = sm.save (
		    c, '{0}/bg_tsds'.format (self.mode_dir),
		    test_ext_deg, dec_deg, seed)
        prush ('->', filename)

    @command
    def collect_bg_trials (self):
        """Collect bg_trials dict and cache it."""
        prush ('Collecting bg trials...')
        if self.opts.cat:
            if not self.opts.weights:
                self.opts.weights = 'equal'
            bg_tsds = bk.get_all (
                '{0}/cats/{1}/{2}/bg_tsds'.format (self.mode_dir,self.opts.cat,self.opts.weights),
                '*.chi2')
            saving (bg_tsds, '{0}/cats/{1}/{2}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        else:
            bg_tsds = bk.get_all (
                '{0}/bg_tsds'.format (self.mode_dir),
                '*.chi2')
            saving (bg_tsds, '{0}/bg_tsds.dict'.format (self.mode_dir))
        return bg_tsds

    @command
    def bg_plots (self):
        """
        Make plots for TS distribution. If allsky, also plot eta and dof across sky.
        """
        if self.opts.cat:
            bg_tsds = loading ('{0}/cats/{1}/{2}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}/plots/bg_tsds'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            fig, ax = pfig()
            tsd = bg_tsds[self.opts.cat][self.conf.root.n_year]
            x = np.linspace(1e-3,20,100)
            h = tsd.get_hist(bins=100).normalize()
            histlite.plot1d (ax, h, label = '{} background trials'.format(tsd.n_total)) ; plt.plot (x,tsd.eta*tsd.chi2.pdf(x), color= 'k', ls='--',
                            label = r'$\tilde{\chi}^2 (n_{dof} = $' + str(np.round(tsd.chi2.args[0],2)) + r', $\eta = $' + str(np.round(tsd.eta,2)) + ')')
            ax.semilogy()
            ax.set_title("bg TS distribution - {0}yr {1}".format(self.conf.root.n_year,self.opts.cat))
            ax.set_ylim (0, max(tsd.values))
            ax.set_xlim (0, h.bins[0][-1])
            plt.legend(loc='upper right', prop=propsmall)
            ax.set_xlabel('TS')
            savingfig (fig, plot_dir, 'tsdist')  
            plt.close(fig)
        else:
            bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))[0]
            plot_dir  = misc.ensure_dir ('{0}/plots/bg_tsds'.format (self.mode_dir))
            sindec = np.sin (np.radians (sorted (bg_tsds)))
            etas, ndofs = np.array ([(bg_tsds[d].eta, bg_tsds[d].ndof) for d in sorted (bg_tsds)]).T
            figall = plt.figure(1) 
            plt.clf() ; plt.plot (sindec, etas, label = 'eta') ; plt.plot (sindec, ndofs, label = 'ndof')
            plt.ylim (0, 2)
            plt.legend(loc='upper right', prop=propsmall)
            plt.xlabel(r'$\sin(\delta)$')
            savingfig (figall, plot_dir, 'bg_tsds_eta_ndof')  
            figall, ax = pfig()
            plt.close(figall)
            #now tsdist for each dec
            for tsd, dec in zip(bg_tsds.values(), bg_tsds.keys()):
              fig, ax = pfig()
              plot_dir  = misc.ensure_dir ('{0}/plots/bg_tsds/{1}'.format (self.mode_dir,dec))
              x = np.linspace(1e-3,20,100)
              h = tsd.get_hist(bins=100).normalize()
              histlite.plot1d (ax, h, label = '{} background trials'.format(tsd.n_total)) ; plt.plot (x,tsd.eta*tsd.chi2.pdf(x), color= 'k', ls='--',
                            label = r'$\tilde{\chi}^2 (n_{dof} = $' + str(np.round(tsd.chi2.args[0],2)) + r', $\eta = $' + str(np.round(tsd.eta,2)) + ')')
              ax.semilogy()
              ax.set_title("bg TS distribution - {0}yr - dec {1}".format(self.conf.root.n_year, dec))
              ax.set_ylim (1e-6, 10)
              ax.set_xlim (0, 20)
              plt.legend(loc='upper right', prop=propsmall)
              ax.set_xlabel('TS')
              savingfig (fig, plot_dir, 'tsdist')  
              plt.close(fig)

    @command
    def submit_do_n_sig (self):
        """
        Submit n_sig stacking jobs for some n-sigma and beta values.
        """
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_n_sig_stacking/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        spectra = [
                (2, np.inf)
                #(2.5, np.inf),
                #(2, 5),
                #(3, np.inf)
        ]
        i_job = 0
        cat = self.opts.cat
        #Weighting
        weights = self.opts.weights
        if not weights:
            weights = 'equal'

        for (src_gamma, src_cutoff) in spectra:
                command = '{} {} do_n_sig --cat {} '\
                        ' --weights={} ' \
                        ' --conf-dirs={} ' \
                        ' --n-trials={}' \
                        ' --src-gamma={:4.2f} ' \
                        ' --src-cutoff={} ' \
                        ' --sigma={:.0f}' \
                        ' --beta={:03.1f}' \
                        ' --seed={}'.format (
                            env_shell,
                            this_script,
                            cat,
                            weights,
                            confs,
                            self.opts.n_trials,
                            src_gamma,
                            src_cutoff,
                            self.opts.sigma,
                            self.opts.beta,
                            i_job,
                            )
                label = 'do_n_sig__' \
                        'cat_{0}__' \
                        'src_gamma_{1:4.2f}__' \
                        'src_cutoff_{2}__' \
                        'seed_{3:08d}'.format (
                                cat,
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
    def do_n_sig (self):
        """Calculate n_sig for given n-sigma in beta fraction of trials."""
        # get parameters from command line
        prush ('Getting parameters from commandline...')
        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
            if 'extension' in params.keys():
                test_ext = params['extension']
                src_ext = params['extension']
            else:
                test_ext = np.zeros_like(params['dec'])
                src_ext = np.zeros_like(params['dec'])
            if self.opts.weights == 'equal':
                weights = np.ones_like(params['dec'])
            else:
                weights = params[str(self.opts.weights)]
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            test_ext_deg = [ext for ext in test_ext]
        else:
            dec_deg = self.opts.test_dec
            test_ext_deg = self.opts.test_ext
            test_ext = test_ext_deg / 180.*pi
            ra_deg = 0
            weights = 1. #put in to eliminate redundancy for ps injection 
            src_ext_deg = self.opts.src_ext
            src_ext = src_ext_deg / 180.*pi
        dec = dec_deg / 180.*pi
        ra = ra_deg / 180.*pi
        gamma = self.opts.src_gamma
        thresh = self.opts.src_thresh
        cutoff = self.opts.src_cutoff
        sigma = self.opts.sigma
        beta = self.opts.beta
        n_max = 7 if sigma < 3 else 12
        n_sig_guess = 15 if sigma < 3 else 45
        use_fit = sigma >= 3
        eps = 1e-5 if use_fit else 0
        #Use random seed
        np.random.seed(seed)

        # get the ts threshold
        if self.opts.cat:
          self.analysis.ps_sindec_width = .05
          bg_tsd = self.bg_tsds
        else:
          self.analysis.ps_sindec_width = .01# = .2 if dec_deg < 60 else 0.05
          bg_tsd = bk.get_best (self.bg_tsds, test_ext, dec_deg)

        ts = bg_tsd.isf (stats.norm.sf (sigma) + eps, fit=use_fit)

        prush ('Setting up injection...')
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdfs_energy_sig['IC86'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, gamma, weights = weights,
                                        energy_range=(10**thresh, 10**cutoff),
                                        sigsub = False)#self.conf.root.sigsub)
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas))

        prush ('Finding n_inj...')
        prush ('- ts > {}'.format (ts))
        prush ('- beta = {0:.3f}'.format (beta))
        if self.opts.cat:
            prush ('- cat = {}'.format (self.opts.cat))
        else:
            prush ('- dec = {0:.3f}'.format (dec_deg))
            prush ('- test_ext = {0:.3f}'.format (test_ext_deg))
            prush ('- src_ext = {0:.3f}'.format (src_ext_deg))
        prush ('- gamma = {0:.3f}'.format (gamma))
        prush ('- thresh = {0:.3f}'.format (thresh))
        prush ('- cutoff = {0:.3f}'.format (cutoff))
        prush ('- ns_bounds = {0:.3f}, {1:.3f}'.format (*self.ns_bounds))
        prush ('- n_sig_guess = {}'.format (n_sig_guess))

        sens = tr.get_sens (ts, beta, n_batch=200, tol=.04,full_output=True)
        flux100TeV = ps.to_flux (sens['n_sig_best'], 100, 1e3) # E0^2 Phi(E0) in TeV/cm^2/s
        fluxGeV = ps.to_flux (sens['n_sig_best']) # E0^2 Phi(E0) in TeV/cm^2/s
        prush ('Sensitivity flux is E0^2 Phi(E0) '
               '= {:.3e} TeV/cm^2/s'.format (flux100TeV))
        prush ('Sensitivity flux is Phi(1GeV) '
               '= {:.3e}'.format (fluxGeV))

        outname = 'sens'
        outname_ts = 'tsdist'

        if self.opts.cat:
          sm = bk.SavingModel (
                'sigma/beta/gamma/thresh/cutoff/sens',
                '{:1.0f}/{:3.1f}/{:4.2f}/{:.2f}/'
                '{:+07.2f}/{}.pickle'
                )
          filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/cats/{1}/{2}/sig_int/'.format (self.mode_dir,self.opts.cat,self.opts.weights),
                sigma, beta,
                gamma, thresh, cutoff,
                outname)
          filename_ts = sm.save (
                sens['tss'],
                '{0}/cats/{1}/{2}/sig_int/'.format (self.mode_dir,self.opts.cat,self.opts.weights),
                sigma, beta,
                gamma, thresh, cutoff,
                outname_ts)
        else:
          sm = bk.SavingModel (
            'sigma/beta/test_ext/src_ext/gamma/thresh/cutoff/dec/sens',
            '{:1.0f}/{:3.1f}/{:04.1f}/{:04.1f}/{:4.2f}/{:.2f}/{:.2f}/'
            '{:+07.2f}/{}.pickle'
          )
          filename = sm.save (
                (sens['n_sig_best'], fluxGeV, flux100TeV),
                '{0}/sig_int/'.format (self.mode_dir),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname)
          filename_ts = sm.save (
                sens['tss'],
                '{0}/sig_int/'.format (self.mode_dir),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname_ts)
        prush ('->', filename)
        prush ('->', filename_ts)

    @command
    def collect_sig_int (self):
        """Collect signal injection results and cache."""
        prush ('Collecting signal results...')
        if self.opts.cat:
            if not self.opts.weights:
              self.opts.weights = 'equal'
            d1 = bk.get_all (
                '{0}/cats/{1}/{2}/sig_int'.format (self.mode_dir,self.opts.cat,self.opts.weights),
                'sens.pickle',
                lambda x: x[0])
            saving (d1, '{0}/cats/{1}/{2}/sig_int.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        else:
            d1 = bk.get_all (
                '{0}/sig_int'.format (self.mode_dir),
                'sens.pickle',
                lambda x: x[0])
            saving (d1, '{0}/sig_int.dict'.format (self.mode_dir))
        bg_tsds = cache.load (
            '{0}/cats/{1}/{2}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        ts,ns,gamma = cache.load ('{0}/cats/{1}/{2}/one_ts_truth/one_ts.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights)) 
        tsmed, eta, ndof = float(bg_tsds.median()), bg_tsds.eta, bg_tsds.ndof
        result_dir = misc.ensure_dir('/data/i3home/brelethford/csky/stacktest/{0}/results/{1}/'.format (self.opts.cat,self.opts.weights))
        sensresults = d1[0][0.9][2.0][2.0]['+000inf']
        discresults = d1[5][0.5][2.0][2.0]['+000inf']
        outfile = open(result_dir+"sensdisc.txt","w")
        outfile.write('single unscrambled test: TS = {0:.2f}, fitted ns = {1:.2f}, fitted gamma = {2:.2f}'.format(ts,ns,gamma))
        outfile.write('\nbackground TS median = {0:.2f}, eta = {1:.2f}, ndof = {2:.2f}'.format(tsmed,eta,ndof))
        outfile.write('\nsens mu = {0:.2f}, sens flux (TeV) = {1:5.3}'.format(sensresults[0],sensresults[-1]))
        outfile.write('\ndisc mu = {0:.2f}, disc flux (TeV) = {1:5.3}'.format(discresults[0],discresults[-1]))
        outfile.close()

        return d1

    @command
    def sig_plots (self):
        """
        Make plots for TS distribution for injected trials. If allsky, also plot eta and dof across sky.
        """
        if self.opts.cat:
            if not self.opts.weights:
              self.opts.weights = 'equal'
            bg_tsd = loading ('{0}/cats/{1}/{2}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            sens_ts = loading ('{0}/cats/{1}/{2}/sig_int/0/0.9/2.00/2.00/+000inf/tsdist.pickle'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            disc_ts = loading ('{0}/cats/{1}/{2}/sig_int/5/0.5/2.00/2.00/+000inf/tsdist.pickle'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}/plots/sig_int'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            result_dir = misc.ensure_dir('/data/i3home/brelethford/csky/stacktest/{0}/results/{1}/'.format (self.opts.cat,self.opts.weights))
            fig, ax = pfig()
            x = np.linspace(1e-3,20,100)
            h_bg = bg_tsd.get_hist(bins=100).normalize()
            h_sens = histlite.hist(sens_ts, bins=100).normalize()
            h_disc = histlite.hist(disc_ts, bins=100).normalize()
            histlite.plot1d (ax, h_bg, label = '{} background trials'.format(bg_tsd.n_total)) ; plt.plot (x,bg_tsd.eta*bg_tsd.chi2.pdf(x), color= 'k', ls='--',
                            label = r'$\tilde{\chi}^2 (n_{dof} = $' + str(np.round(bg_tsd.chi2.args[0],2)) + r', $\eta = $' + str(np.round(bg_tsd.eta,2)) + ')')
            histlite.plot1d (ax, h_sens, color = 'red', label = '{} sens injection trials'.format(len(sens_ts)))
            histlite.plot1d (ax, h_disc, color = 'green', label = '{} disc injection trials'.format(len(disc_ts)))
            ax.semilogy()
            ax.set_title("TS distributions - {0}yr {1}".format(self.conf.root.n_year,self.opts.cat))
            ax.set_ylim (1e-4, 10)
            ax.set_xlim (0, 40)
            plt.legend(loc='upper right', prop=propsmall)
            ax.set_xlabel('TS')
            savingfig (fig, plot_dir, 'sig_tsdist')  
            savingfig (fig, result_dir, 'sig_tsdist')  
            plt.close(fig)
        else:
            bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))[0]
            plot_dir  = misc.ensure_dir ('{0}/plots/bg_tsds'.format (self.mode_dir))
            sindec = np.sin (np.radians (sorted (bg_tsds)))
            etas, ndofs = np.array ([(bg_tsds[d].eta, bg_tsds[d].ndof) for d in sorted (bg_tsds)]).T
            figall = plt.figure(1) 
            plt.clf() ; plt.plot (sindec, etas, label = 'eta') ; plt.plot (sindec, ndofs, label = 'ndof')
            plt.ylim (0, 2)
            plt.legend(loc='upper right', prop=propsmall)
            plt.xlabel(r'$\sin(\delta)$')
            savingfig (figall, plot_dir, 'bg_tsds_eta_ndof')  
            figall, ax = pfig()
            plt.close(figall)
            #now tsdist for each dec
            for tsd, dec in zip(bg_tsds.values(), bg_tsds.keys()):
              fig, ax = pfig()
              plot_dir  = misc.ensure_dir ('{0}/plots/bg_tsds/{1}'.format (self.mode_dir,dec))
              x = np.linspace(1e-3,20,100)
              h = tsd.get_hist(bins=100).normalize()
              histlite.plot1d (ax, h, label = '{} background trials'.format(tsd.n_total)) ; plt.plot (x,tsd.eta*tsd.chi2.pdf(x), color= 'k', ls='--',
                            label = r'$\tilde{\chi}^2 (n_{dof} = $' + str(np.round(tsd.chi2.args[0],2)) + r', $\eta = $' + str(np.round(tsd.eta,2)) + ')')
              ax.semilogy()
              ax.set_title("bg TS distribution - {0}yr - dec {1}".format(self.conf.root.n_year, dec))
              ax.set_ylim (1e-6, 10)
              ax.set_xlim (0, 20)
              plt.legend(loc='upper right', prop=propsmall)
              ax.set_xlabel('TS')
              savingfig (fig, plot_dir, 'tsdist')  
              plt.close(fig)

if __name__ == '__main__':
    self = teststack ()
    self.run ()
    try:
        import pandas as pd
    except:
        pass


hl = histlite


