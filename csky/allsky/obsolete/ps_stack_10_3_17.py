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
import subprocess

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
import pylab as py

from icecube import histlite
from icecube.umdtools.apps import Timed, Driver, command
from icecube.umdtools.conf import Configurator
from icecube.umdtools.submitter import Submitter
from icecube.umdtools.arrays import Arrays, combined_arrays
from icecube.umdtools.vars_class import Vars
from icecube.umdtools import cache, misc#, quiet_healpy
from skylab.ps_llh import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector import PointSourceInjector
from skylab.llh_models import ClassicLLH, EnergyLLH

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

hem = np.sin(np.radians(-5.))

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

def savefig (fig, outdir, namebase, exts='png', specialpdf=False):
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
    public_dir = '/home/brelethford/public_html/csky/' + d.split('brelethford/')[-1]
    print (public_dir+'/'+n)
    subprocess.call(['ssh', 'cobalt06', 'mkdir', '-p', '{0}/'.format(public_dir)])
    subprocess.call(["scp", '-r','{0}/{1}.png'.format (d,n), 'brelethford@cobalt06:{}'.format(public_dir)])
    #subprocess.call(["scp", '-r','{0}/{1}.pdf'.format (d,n), 'brelethford@cobalt06:{}'.format(public_dir)])

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

soft_colors = ['#004466', '#d06050', '#2aca80', '#dd9388', '#caca68']

class Csky6yr (Timed, Driver):

    def __init__ (self):
        Timed.__init__ (self)
        Driver.__init__ (self)

    def run (self, arglist=[]):

        usage = '%prog {[options]} [commands]\n' + self._command_help
        self.parser = parser = OptionParser (usage=usage)

        parser.add_option ('-c', '--conf-dirs', dest='conf_dirs',
                default=os.path.abspath ('conf_7yr'), metavar='DIRS',
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

        parser.add_option ('--test-ra', dest='test_ra',
                default=None, type=float, metavar='ra',
                help='test point RA in degrees')

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
                default=50., type=float, metavar='N',
                help='perform N jobs (with --n-trials each)')

        parser.add_option ('--n-trials', dest='n_trials',
                default=200., type=float, metavar='N',
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

        parser.add_option ('--dec-order', dest='dec_order',
                default=False, action='store_true',
                help='orders unblinding test to add sources according to dec')

        parser.add_option ('--weight-order', dest='weight_order',
                default=False, action='store_true',
                help='orders unblinding test to add sources with decreasing source weight')

        parser.add_option ('--ns', dest='ns',
                default=False, action='store_true',
                help='order stacking sources from north to south, rather than south to north')

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
          if 'corrected' in self.mode.words[0]:
            print ("For 'corrected' npz files")
            if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/datasets_corrected'
                print ("Stefan's data, included pullfactor")
            else:
                datasetfolder = self.root_dir+'/data/datasets_nopullfactor_corrected'
                print ("Stefan's data, no pullfactor")

          elif self.conf.root.npzdata: #whose data is being used? mine or stefan's?

            if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/datasets'
                print ("Stefan's data, included pullfactor")
            else:
                datasetfolder = self.root_dir+'/data/datasets_nopullfactor'
                print ("Stefan's data, no pullfactor")

          else:
            if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/my_datasets'
                print ("My data, included pullfactor")
            else:
                datasetfolder = self.root_dir+'/data/my_datasets_nopullfactor'
                print ("My data, no pullfactor")

          self._psdata = loading ('{0}/psdata.vars'.format (datasetfolder))

        return self._psdata

    @property
    def datasets (self):
        try:
            return self._datasets
        except:
            samples = self.conf.root.samplelist.split()
            if 'corrected' in self.mode.words[0]:
              if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/datasets_corrected'
                print ("For 'corrected' npz files")
              else:
                datasetfolder = self.root_dir+'/data/datasets_nopullfactor_corrected'
                print ("For 'corrected' npz files - nopullfactor")
            elif self.conf.root.npzdata: #whose data is being used? mine or stefan's?
              if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/datasets'
                print ("Stefan's data, included pullfactor")
              else:
                datasetfolder = self.root_dir+'/data/datasets_nopullfactor'
                print ("Stefan's data, no pullfactor")
            else:
              if self.conf.root.pullfactor:
                datasetfolder = self.root_dir+'/data/my_datasets'
                print ("My data, included pullfactor")
              else:
                datasetfolder = self.root_dir+'/data/my_datasets_nopullfactor'
                print ("My data, no pullfactor")

            samplenames = ['{}/{}.dataset'.format (datasetfolder, i) for i in samples]
            self._datasets = map (loading, samplenames)

            return self._datasets

    @property
    def analysis (self):
        try:
            return self._analysis
        except:
            #if self.conf.
            self._analysis = loading ('{0}/ps.analysis'.format (self.mode_dir))
            return self._analysis

    #like analysis, but for skylab (and only llhmodel - need to create psllh object
    @property
    def llhs (self):
        try:
            return self._llhs
        except:
            self._llhs = loading ('{0}/ps.llhs'.format (self.mode_dir))
            return self._llhs

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
            return loading ('{0}/cats/{1}/sig_int.dict'.format (self.mode_dir,self.opts.cat))
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
        #Is pullcorrected data used?
        pullfactor = self.conf.root.pullfactor

        ids_exp = ['{}'.format (s)
               for s in 'IC40 IC59 IC79_noMESE IC79b IC86_noMESE IC86-2012_noMESE IC86-2013_noMESE IC86-2014_noMESE MESE MESE_followup'.split()]
        filenames_exp = ['{}/data/stefandata/{}_exp.npy'.format (self.root_dir, i)
                     for i in ids_exp]

        ids_mc = ['{}'.format (s)
               for s in 'IC40 IC59 IC79_noMESE IC79b_noMESE IC86_noMESE IC86-2012_noMESE MESE'.split()]

        if 'corrected' in self.mode.words[0]:
          filenames_mc = ['{}/data/stefandata/{}_corrected_MC.npy'.format (self.root_dir, i)
                     for i in ids_mc]
        else:
          filenames_mc = ['{}/data/stefandata/{}_MC.npy'.format (self.root_dir, i)
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
        nus = [nu40,nu59,nu79sirin,nu79,nu86,nu2012,nuMESE] 

        for nu,mc in zip(nus,contents_mc): 
            nu.true_energy =  mc['trueE']
            nu.true_zenith =  pi/2.+mc['trueDec']
            nu.true_azimuth = mc['trueRa']
            nu.energy =       10**mc['logE']
            nu.zenith =       np.pi/2.+np.arcsin(mc['sinDec'])
            nu.azimuth =      mc['ra']
            nu.oneweight =    mc['ow']
            if pullfactor:
              nu.sigma =        mc['sigma']
            else:
              nu.sigma =        mc['sigma']/1.1774
              print("Reintroducing pullcorrect factor for earlier years...")
            if 'dist' in mc.dtype.names:
              nu.dist =       mc['dist']
            else:
              nu.dist =       np.zeros_like(nu.energy)

        mu40 = psdata.mu40 = Arrays ()
        mu59 = psdata.mu59 = Arrays ()
        mu79 = psdata.mu79 = Arrays ()
        mu79sirin = psdata.mu79sirin = Arrays ()
        mu86 = psdata.mu86 = Arrays ()
        mu2012 = psdata.mu2012 = Arrays ()
        mu2013 = psdata.mu2013 = Arrays ()
        mu2014 = psdata.mu2014 = Arrays ()
        muMESE_first = Arrays ()
        muMESEfollowup = Arrays ()

        mus = [mu40,mu59,mu79sirin,mu79,mu86,mu2012,mu2013,mu2014,muMESE_first,muMESEfollowup] 

        for mu,exp in zip(mus,contents_exp): 
            mu.energy =   10**exp['logE']
            mu.zenith =   np.pi/2.+np.arcsin(exp['sinDec'])
            mu.azimuth =  exp['ra']
            if pullfactor:
              mu.sigma =  exp['sigma']
            else:
              mu.sigma =  exp['sigma']/1.1774
            if 'dist' in exp.dtype.names:
              mu.dist =   exp['dist']
            else:
              mu.dist =   np.zeros_like(mu.energy)

        #also need to add the followup to the exp 
        muMESE = psdata.muMESE = combined_arrays((muMESE_first,muMESEfollowup))

        mus = [mu40,mu59,mu79sirin,mu79,mu86,mu2012,mu2013,mu2014,muMESE] 

        for d in nus+mus:
            d.apply_cut ((1e1 < d.energy) & (d.energy < 1e10))
            d.dec = d.zenith - pi/2
            d.ra = d.azimuth
            if 'true_zenith' in d:
                d.true_dec = d.true_zenith - pi/2
                d.true_ra = d.true_azimuth

        #Where does this data belong?
        if 'corrected' in self.mode.words[0]:
          print ("For 'corrected' npz files")
          if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/datasets_corrected'
            print ("Stefan's data, included pullfactor")
          else:
            datasetfolder = self.root_dir+'/data/datasets_nopullfactor_corrected'
            print ("Stefan's data, no pullfactor")

        else: #whose data is being used? mine or stefan's?
          if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/datasets'
            print ("Stefan's data, included pullfactor")
          else:
            datasetfolder = self.root_dir+'/data/datasets_nopullfactor'
            print ("Stefan's data, no pullfactor")

        saving (psdata, '{0}/psdata.vars'.format (misc.ensure_dir(datasetfolder)))

    @command
    def build_mydata (self):
        """
        Build my ps data.
        My own pullcorrections, no extra cuts, decorr.
        """
        psdata = Vars ()
        pullfactor = self.conf.root.pullfactor
        #Currently leaving MESE out of these datasets.
        ids = ['{}'.format (s)
               for s in 'IC40 IC59 IC79/noMESE sirin_IC79/noMESE IC86I/noMESE 3yr/noMESE MESE'.split()]
        filenames_exp = ['{}/data/mydata/{}/exp.pickle'.format (self.root_dir, i)
                     for i in ids]
        filenames_mc = ['{}/data/mydata/{}/mc.pickle'.format (self.root_dir, i)
                     for i in ids]

        prush ('Loading MC and exp datasets...')
        contents_exp = map (cache.load, filenames_exp)
        contents_mc = map (cache.load, filenames_mc)

        nu40 = psdata.nu40 = Arrays ()
        nu59 = psdata.nu59 = Arrays ()
        nu79 = psdata.nu79 = Arrays ()
        nu79sirin = psdata.nu79sirin = Arrays ()
        nu86 = psdata.nu86 = Arrays ()
        nu3yr = psdata.nu3yr = Arrays ()
        nuMESE = psdata.nuMESE = Arrays ()
        nus = [nu40,nu59,nu79,nu79sirin,nu86,nu3yr,nuMESE]

        for nu,mc in zip(nus,contents_mc):
            nu.true_energy =  mc['trueE']
            nu.true_zenith =  pi/2.+mc['trueDec']
            nu.true_azimuth = mc['trueRa']
            nu.energy =       10**mc['logE']
            nu.zenith =       np.pi/2.+np.arcsin(mc['sinDec'])
            nu.azimuth =      mc['ra']
            nu.oneweight =    mc['ow']
            if pullfactor:
              nu.sigma =  mc['sigma']*1.1774
            else:
              nu.sigma =  mc['sigma']

        mu40 = psdata.mu40 = Arrays ()
        mu59 = psdata.mu59 = Arrays ()
        mu79 = psdata.mu79 = Arrays ()
        mu79sirin = psdata.mu79sirin = Arrays ()
        mu86 = psdata.mu86 = Arrays ()
        mu3yr = psdata.mu3yr = Arrays ()
        muMESE = psdata.muMESE = Arrays ()
        mus = [mu40,mu59,mu79,mu79sirin,mu86]
        mus = [mu40,mu59,mu79,mu79sirin,mu86,mu3yr,muMESE]
        for mu,exp in zip(mus,contents_exp):
            mu.energy =   10**exp['logE']
            mu.zenith =   np.pi/2.+np.arcsin(exp['sinDec'])
            mu.azimuth =  exp['ra']
            if pullfactor:
              mu.sigma =  exp['sigma']*1.1774
            else:
              mu.sigma =  exp['sigma']

        for nu,mu in zip(nus,mus):
          nu.apply_cut ((1e1 < nu.energy) & (nu.energy < 1e10))
          mu.apply_cut ((1e1 < mu.energy) & (mu.energy < 1e10))
          for d in [nu, mu]:
            d.dec = d.zenith - pi/2
            d.ra = d.azimuth
            if 'true_zenith' in d:
                d.true_dec = d.true_zenith - pi/2
                d.true_ra = d.true_azimuth

        #Put this data in appropriate folder, according to pullfactor.
        if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/my_datasets'
            print ("My data, included pullfactor")
        else:
            datasetfolder = self.root_dir+'/data/my_datasets_nopullfactor'
            print ("My data, no pullfactor")

        saving (psdata, '{0}/psdata.vars'.format (misc.ensure_dir(datasetfolder)))

    @command
    def setup_datasets (self):
        """Create Dataset."""
        #define datafolder - my or stefan's data
        if 'corrected' in self.mode.words[0]:
          print ("For 'corrected' npz files")
          if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/datasets_corrected/'
          else:
            datasetfolder = self.root_dir+'/data/datasets_nopullfactor_corrected/'
            print ("Stefan's data, no pullfactor")
        elif self.conf.root.npzdata: #whose data is being used? mine or stefan's?
          if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/datasets/'
            print ("Stefan's data, included pullfactor")
          else:
            datasetfolder = self.root_dir+'/data/datasets_nopullfactor/'
            print ("Stefan's data, no pullfactor")
        else:
          if self.conf.root.pullfactor:
            datasetfolder = self.root_dir+'/data/my_datasets/'
            print ("My data, included pullfactor")
          else:
            datasetfolder = self.root_dir+'/data/my_datasets_nopullfactor/'
            print ("My data, no pullfactor")

        nu40 = self.psdata.nu40
        nu59 = self.psdata.nu59
        nu79 = self.psdata.nu79
        nu79sirin = self.psdata.nu79sirin
        nu86 = self.psdata.nu86
        mu40 = self.psdata.mu40
        mu59 = self.psdata.mu59
        mu79 = self.psdata.mu79
        mu79sirin = self.psdata.mu79sirin
        mu86 = self.psdata.mu86
        livetime40  = 375.539 * 86400
        livetime59  = 348.138 * 86400
        livetime79  = 315.506 * 86400
        livetime86 = 332.61 *  86400
        livetime86II = 330.38 * 86400
        livetime86III = 359.95 * 86400
        livetime86IV = 367.21 * 86400
        livetime3yr = livetime86II + livetime86III + livetime86IV
        livetimeMESE = 1715 * 86400
        dataset40 = ana.Dataset ('IC40', livetime40, sig=nu40, data=mu40)
        dataset59 = ana.Dataset ('IC59', livetime59, sig=nu59, data=mu59)
        dataset79 = ana.Dataset ('IC79', livetime79, sig=nu79, data=mu79)
        dataset79sirin = ana.Dataset ('IC79sirin', livetime79, sig=nu79sirin, data=mu79sirin)
        dataset86 = ana.Dataset ('IC86', livetime86, sig=nu86, data=mu86)

        saving (dataset40, '{}IC40.dataset'.format (datasetfolder))
        saving (dataset59, '{}IC59.dataset'.format (datasetfolder))
        saving (dataset79, '{}IC79.dataset'.format (datasetfolder))
        saving (dataset79sirin, '{}IC79sirin.dataset'.format (datasetfolder))
        saving (dataset86, '{}IC86.dataset'.format (datasetfolder))
        #add remaining years if using stefan's data
        if self.conf.root.npzdata:
            nu2012 = self.psdata.nu2012
            mu2012 = self.psdata.mu2012
            mu2013 = self.psdata.mu2013
            mu2014 = self.psdata.mu2014
            dataset86II = ana.Dataset ('IC86II', livetime86II, sig=nu2012, data=mu2012)
            dataset86III = ana.Dataset ('IC86III', livetime86III, sig=nu2012, data=mu2013)
            dataset86IV = ana.Dataset ('IC86IV', livetime86IV, sig=nu2012, data=mu2014)
            saving (dataset86II, '{}IC86II.dataset'.format (datasetfolder))
            saving (dataset86III, '{}IC86III.dataset'.format (datasetfolder))
            saving (dataset86IV, '{}IC86IV.dataset'.format (datasetfolder))
        #if not npz, use the 3yr from mydata.
        else:
            nu3yr = self.psdata.nu3yr
            mu3yr = self.psdata.mu3yr
            dataset3yr = ana.Dataset ('IC86_3yr', livetime3yr, sig=nu3yr, data=mu3yr)
            saving (dataset3yr, '{}IC86_3yr.dataset'.format (datasetfolder))
        nuMESE = self.psdata.nuMESE
        muMESE = self.psdata.muMESE
        datasetMESE = ana.Dataset ('MESE', livetimeMESE, sig=nuMESE, data=muMESE)
        saving (datasetMESE, '{}MESE.dataset'.format (datasetfolder))

    @command
    def setup_ana (self):
        """
        Set up main analysis objects.
        """
        #Let's determine whose binning and which years we're using.
        n_year = self.conf.root.n_year
        print ("Making analysis object for " + str(n_year) + " years of data")

        datasets = self.datasets
        analysis = ana.Analysis(datasets, build_pdfs=False)
        hem = np.sin(np.radians(-5.))
        samples = self.conf.root.samplelist.split()
        if self.conf.root.sledge:
          cap_threshold = 0.5
        else:
          cap_threshold = 0

        binning = self.conf.root.binning
        #binning can be 'npz_restrict' (meaning restricted), 'npz_full' (unrestricted npz), or 'uniform'.
        llhmodels = []

        #Now we'll populate both the analysis object and llhmodels, depending on the binning and samples.
        
        #Stefan's Binning: note that only logE_range and Ebins change depending on the binning used (restricted or full) 
        if binning != 'uniform':
          #IC40
          dec_bins_40 = np.unique(np.concatenate([
                       np.linspace(-1., -0.25, 10 + 1),
                       np.linspace(-0.25, 0.0, 5 + 1),
                       np.linspace(0.0, 1., 10 + 1),
                       ]))
          dec_bins_logE_40 = np.unique(np.concatenate([
                       np.linspace(-1., -0.25, 10 + 1),
                       np.linspace(-0.25, 0.0, 10 + 1),
                       np.linspace(0.0, 1., 10 + 1),
                       ]))
          Ebins_40_r = np.linspace(2.6, 8.4, 35 + 1)

          #IC59
          dec_bins_59 = np.unique(np.concatenate([
                               np.linspace(-1., -0.95, 2 + 1),
                               np.linspace(-0.95, -0.25, 25 + 1),
                               np.linspace(-0.25, 0.05, 15 + 1),
                               np.linspace(0.05, 1., 10 + 1),
                               ]))
          
          dec_bins_logE_59 = np.unique(np.concatenate([
                               np.linspace(-1., 0.05, 20 + 1),
                               np.linspace(0.05, 1., 10 + 1),
                               ]))

          Ebins_59_r = np.linspace(2.5, 8.5, 50 + 1)

          #IC79sirin
          dec_bins_79sirin = np.unique(np.concatenate([
                       np.linspace(-1., -0.75, 10 + 1),
                       np.linspace(-0.75, 0.0, 15 + 1),
                       np.linspace(0.0, 1., 20 + 1),
                       ]))
          dec_bins_logE_79sirin = dec_bins_79sirin
          Ebins_79sirin_r = np.linspace(2.5, 8.5, 40 + 1)

          #IC79
          dec_bins_79 = np.unique(np.concatenate([
                               np.linspace(-1., -0.75, 10 + 1),
                               np.linspace(-0.75, 0.0, 15 + 1),
                               np.linspace(0.0, 1., 20 + 1),
                               ]))
          dec_bins_logE_79 = dec_bins_79
          Ebins_79_r = np.linspace(2.55, 8.55, 40 + 1)

          #IC86
          dec_bins_86 = np.unique(np.concatenate([
                             np.linspace(-1., -0.2, 10 + 1),
                             np.linspace(-0.2, hem, 4 + 1),
                             np.linspace(hem, 0.2, 5 + 1),
                             np.linspace(0.2, 1., 10),
                             ]))
          dec_bins_logE_86 = dec_bins_86 
          Ebins_86_r = np.linspace(1.5, 9.5, 40 + 1)

          #skylab always has restricted, but csky's is dependent on the binning choice.
          if binning == 'npz_restrict':
            print ('binning = restricted')
            Ebins_40 = Ebins_40_r 
            Ebins_59 = Ebins_59_r 
            Ebins_79sirin = Ebins_79sirin_r 
            Ebins_79 = Ebins_79_r 
            Ebins_86 = Ebins_86_r 
            logE_range_40 = (2.6, 8.4)
            logE_range_59 = (2.5, 8.5)
            logE_range_79sirin = (2.5, 8.5)
            logE_range_79 = (2.55, 8.55)
            logE_range_86 = (1.5, 9.5)

          elif binning == 'npz_full':
            print ('binning = full')
            Ebins_40 = np.linspace(1., 10., 35 + 1)
            Ebins_59 = np.linspace(2., 9.5, 50 + 1)
            Ebins_79sirin = np.linspace(2., 9., 67 + 1)
            Ebins_79 = np.linspace(2., 9., 67 + 1)
            Ebins_86 = np.linspace(1., 10., 40 + 1)
            logE_range_40 = (1., 10.)
            logE_range_59 = (2., 9.5)
            logE_range_79sirin = (2., 9.)
            logE_range_79 = (2., 9.)
            logE_range_86 = (1., 10.)

          #Fill in analysis object and llhmodel.
          #IC40
          if 'IC40' in samples:
              #skylab
              energy_bins = [Ebins_40_r, dec_bins_logE_40]
              llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins_40)
              llhmodels.append(llhmodel)
              #csky
              analysis.set_dec_pdf(bins=dec_bins_40, fit=True, range=(-1, 1), keys = ['IC40'])
              analysis.set_energy_pdfs (
                sindec_bins=dec_bins_logE_40,
                logenergy_bins=Ebins_40,
                sindec_range=(-1, 1), logenergy_range=logE_range_40,
                cap_threshold=cap_threshold,
                keys=['IC40']
              )

              #IC59
          if 'IC59' in samples:
              #skylab
              energy_bins = [Ebins_59_r, dec_bins_logE_59]
              llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins_59)
              llhmodels.append(llhmodel)
              #csky
              analysis.set_dec_pdf(bins=dec_bins_59, fit=True, range=(-1, 1), keys = ['IC59'])
              analysis.set_energy_pdfs (
                sindec_bins=dec_bins_logE_59,
                logenergy_bins=Ebins_59,
                sindec_range=(-1, 1), logenergy_range=logE_range_59,
                cap_threshold=cap_threshold,
                keys=['IC59']
              )

              #IC79sirin
          if 'IC79sirin' in samples:
              energy_bins = [Ebins_79sirin_r, dec_bins_logE_79sirin]
              llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins_79sirin)
              llhmodels.append(llhmodel)
              #csky
              analysis.set_dec_pdf(bins=dec_bins_79sirin, fit=True, range=(-1, 1), keys = ['IC79sirin'])
              analysis.set_energy_pdfs (
                sindec_bins=dec_bins_logE_79sirin,
                logenergy_bins=Ebins_79sirin,
                sindec_range=(-1, 1), logenergy_range=logE_range_79sirin,
                cap_threshold=cap_threshold,
                keys=['IC79sirin']
              )
              #IC79
          if 'IC79' in samples:
              #skylab
              energy_bins = [Ebins_79_r, dec_bins_logE_79]
              llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins_79)
              llhmodels.append(llhmodel)
              #csky
              analysis.set_dec_pdf(bins=dec_bins_79, fit=True, range=(-1, 1), keys = ['IC79'])
              analysis.set_energy_pdfs (
                sindec_bins=dec_bins_logE_79,
                logenergy_bins=Ebins_79,
                sindec_range=(-1, 1), logenergy_range=logE_range_79,
                cap_threshold=cap_threshold,
                keys=['IC79']
              )

              #IC86
          if 'IC86' in samples:
              #skylab
              energy_bins = [Ebins_86_r, dec_bins_logE_86]
              llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins_86)
              llhmodels.append(llhmodel)
              #csky
              analysis.set_dec_pdf(bins=dec_bins_86, fit=True, range=(-1, 1), keys = ['IC86'])
              analysis.set_energy_pdfs (
                sindec_bins=dec_bins_logE_86,
                logenergy_bins=Ebins_86,
                sindec_range=(-1, 1), logenergy_range=logE_range_86,
                cap_threshold=cap_threshold,
                keys=['IC86']
              )

          '''
          #Stefan's Binning - Full range
          elif binning == 'npz_full':
            print("Using stefan's binning for csky - skylab still restricted!")
            #IC40
            if 'IC40' in samples:
                analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.25, 10 + 1),
                                 np.linspace(-0.25, 0.0, 5 + 1),
                                 np.linspace(0.0, 1., 10 + 1)
                                 ])), fit=True, range=(-1, 1), keys = ['IC40']
                )
                analysis.set_energy_pdfs (
                  sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.25, 10 + 1),
                                 np.linspace(-0.25, 0.0, 10 + 1),
                                 np.linspace(0.0, 1., 10 + 1)
                                 ])),
                  logenergy_bins=np.linspace(1., 10., 35 + 1),
                  sindec_range=(-1, 1), logenergy_range=(1., 10.),
                  cap_threshold=cap_threshold,
                  keys=['IC40']
                )

                #IC59
            if 'IC59' in samples:
                analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.95, 2 + 1),
                                 np.linspace(-0.95, -0.25, 25 + 1),
                                 np.linspace(-0.25, 0.05, 15 + 1),
                                 np.linspace(0.05, 1., 10 + 1)
                                 ])), fit=True, range=(-1, 1), keys = ['IC59']
                )
                analysis.set_energy_pdfs (
                  sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., 0.05, 20 + 1),
                                 np.linspace(0.05, 1., 10 + 1)
                                 ])),
                  logenergy_bins=np.linspace(2., 9.5, 50 + 1),
                  sindec_range=(-1, 1), logenergy_range=(2, 9.5),
                  cap_threshold=cap_threshold,
                  keys=['IC59']
                )

                #IC79
            if 'IC79sirin' in samples:
                analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.75, 10 + 1),
                                 np.linspace(-0.75, 0.0, 15 + 1),
                                 np.linspace(0.0, 1., 20 + 1)
                                 ])), fit=True, range=(-1, 1), keys = ['IC79sirin']
                )
                analysis.set_energy_pdfs (
                  sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.75, 10 + 1),
                                 np.linspace(-0.75, 0.0, 15 + 1),
                                 np.linspace(0.0, 1., 20 + 1)
                                 ])),
                  logenergy_bins=np.linspace(2., 9., 67 + 1),
                  sindec_range=(-1, 1), logenergy_range=(2, 9.),
                  cap_threshold=cap_threshold,
                  keys=['IC79sirin']
                )

            if 'IC79' in samples:
                analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.75, 10 + 1),
                                 np.linspace(-0.75, 0.0, 15 + 1),
                                 np.linspace(0.0, 1., 20 + 1)
                                 ])), fit=True, range=(-1, 1), keys = ['IC79']
                )
                analysis.set_energy_pdfs (
                  sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.75, 10 + 1),
                                 np.linspace(-0.75, 0.0, 15 + 1),
                                 np.linspace(0.0, 1., 20 + 1)
                                 ])),
                  logenergy_bins=np.linspace(2., 9., 67 + 1),
                  sindec_range=(-1, 1), logenergy_range=(2, 9.),
                  cap_threshold=cap_threshold,
                  keys=['IC79']
                )

                #IC86
            if 'IC86' in samples:
                analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.2, 10 + 1),
                                 np.linspace(-0.2, hem, 4 + 1),
                                 np.linspace(hem, 0.2, 5 + 1),
                                 np.linspace(0.2, 1., 10)
                                 ])), fit=True, range=(-1, 1), keys = ['IC86']
                )
                analysis.set_energy_pdfs (
                  sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.2, 10 + 1),
                                 np.linspace(-0.2, hem, 4 + 1),
                                 np.linspace(hem, 0.2, 5 + 1),
                                 np.linspace(0.2, 1., 10)
                                 ])),
                  logenergy_bins=np.linspace(1., 10., 40 + 1),
                  sindec_range=(-1, 1), logenergy_range=(1, 10),
                  cap_threshold=cap_threshold,
                  keys=['IC86']
                )

                  #IC86II-IV
            if 'IC86II' in samples:
                  analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])), range=(-1, 1), keys = ['IC86II']
                  )
                  analysis.set_energy_pdfs (
                    sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])),
                    logenergy_bins=np.linspace(1., 9.5, 50 + 1),
                    fit_dims=(0,1),
                    sindec_range=(-1, 1), logenergy_range=(1, 9.5),
                    keys=['IC86II']
                  )

            if 'IC86III' in samples:
                  analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])), range=(-1, 1), keys = ['IC86III']
                  )
                  analysis.set_energy_pdfs (
                    sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])),
                    logenergy_bins=np.linspace(1., 9.5, 50 + 1),
                    sindec_range=(-1, 1), logenergy_range=(1, 9.5),
                    fit_dims=(0,1),
                    keys=['IC86III']
                  )

            if 'IC86IV' in samples:
                  analysis.set_dec_pdf(bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])), range=(-1, 1), keys = ['IC86IV']
                  )
                  analysis.set_energy_pdfs (
                    sindec_bins=np.unique(np.concatenate([
                                 np.linspace(-1., -0.93, 4 + 1),
                                 np.linspace(-0.93, -0.3, 10 + 1),
                                 np.linspace(-0.3, 0.05, 9 + 1),
                                 np.linspace(0.05, 1., 18 + 1)
                                 ])),
                    logenergy_bins=np.linspace(1., 9.5, 50 + 1),
                    sindec_range=(-1, 1), logenergy_range=(1, 9.5),
                    fit_dims=(0,1),
                    keys=['IC86IV']
                  )

            if 'MESE' in samples:
                  #MESE
                  analysis.set_dec_pdf(bins= np.unique(np.concatenate([np.linspace(-1., -0.93, 4 + 1),
                                       np.linspace(-0.93, hem, 12 + 1)])),
                                       range=(-1, hem), keys = ['MESE']
                  )
                  analysis.set_energy_pdfs (
                    sindec_bins=np.linspace(-1., hem, 4 + 1),
                    logenergy_bins=np.linspace(2., 8.5, 67 + 1), dist_bins=50,
                    sindec_range=(-1, hem), logenergy_range=(2, 8.5),
                    pre_smooth=.5, smooth=4, smooth_bins=[1, 300, 300], keys=['MESE']
                  )
          '''
        #Uniform Binning
        else:
          print("Uniform Binning")
          #skylab
          dec_bins = np.unique(np.linspace(-1., 1, 40 + 1))
          energy_bins = [np.linspace( 2.5,8.5,24+1),dec_bins]
          llhmodel = EnergyLLH(energy_bins, sinDec_bins=dec_bins)
          for sample in samples:
            llhmodels.append(llhmodel)

          #csky
          print ('Currently, no fit applied to uniform binning')
          analysis = ana.Analysis (
              datasets, 
              dec_kw=dict (
                  bins=40+1, range=(-1, 1),
                  keys = tuple(self.conf.root.samplelist.split()),
              ),
              energy_kw=dict (
                  sindec_bins=40+1, logenergy_bins=24+1,
                  logenergy_range=(1, 10),
                  cap_threshold=cap_threshold,
                  keys = tuple(self.conf.root.samplelist.split()),
              )
            )

        #save the csky analysis object
        prush('Saving analysis object...')
        saving (analysis, '{}/ps.analysis'.format (self.mode_dir))

        #configure and save the skylab object with appropriate data
        npz_dir = self.root_dir + '/data/stefandata'

        exp40 = np.load(npz_dir + '/IC40_exp.npy')
        exp59 = np.load(npz_dir + '/IC59_exp.npy')
        exp79 = np.load(npz_dir + '/IC79b_exp.npy')
        exp79sirin = np.load(npz_dir + '/IC79_noMESE_exp.npy')
        exp86 = np.load(npz_dir + '/IC86_noMESE_exp.npy')
          
        mc40 = np.load(npz_dir + '/IC40_corrected_MC.npy')
        mc59 = np.load(npz_dir + '/IC59_corrected_MC.npy')
        mc79 = np.load(npz_dir + '/IC79b_corrected_MC.npy')
        mc79sirin = np.load(npz_dir + '/IC79_noMESE_corrected_MC.npy')
        mc86 = np.load(npz_dir + '/IC86_noMESE_corrected_MC.npy')

        exps = {'IC40':exp40, 'IC59':exp59, 'IC79sirin':exp79sirin, 'IC79':exp79, 'IC86':exp86}
        mcs = {'IC40':mc40, 'IC59':mc59, 'IC79sirin':mc79sirin, 'IC79':mc79, 'IC86':mc86}
        livetimes = {'IC40':375.532,'IC59':348.138, 'IC79sirin':315.506, 'IC79':315.506,'IC86':332.61}

        if not self.conf.root.pullfactor:
          for sample in samples:
            exps[sample]['sigma'] /= 1.1774
            mcs[sample]['sigma'] /= 1.1774
          print("Reintroducing pullcorrect factor for earlier years...")

        #Skylab: We have to constrain gamma to be 2?

        for llhmodel in llhmodels:
          llhmodel.params = {'gamma': (2.0, (1.,4.))}

        llhmodels = dict(zip(samples,llhmodels)) 

        llhs =[]
        #now define ps object (unweighted) for each year.
        for sample in samples:
          exp = exps['{}'.format(sample)] 
          mc = mcs['{}'.format(sample)] 
          ps_llh = PointSourceLLH(exp, mc, livetime=livetimes['{}'.format(sample)], llh_model=llhmodels[sample], mode = 'all', seed=0, scramble=False)
          llhs.append(ps_llh)

        llhs = dict(zip(samples,llhs)) 

        prush('Saving llh object...')
        saving (llhs, '{}/ps.llhs'.format (self.mode_dir))

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
        '''
        for g in np.linspace(1,4,7):#analysis.gammas:
          for year in analysis.datasets.keys():
            plot_dir  = misc.ensure_dir ('{0}/plots/energypdfs/{1}'.format (self.mode_dir,year))
            fig, ax = pfig()#plt.subplots (num=10)
            ax.cla()
            if year != 'MESE':
              histlite.plot2d (ax, analysis.pdf_ratios_energy[year].epdfratios[g].h[:,:,0],
              cmap='seismic', log=True, vmin=1e-4, vmax=1e4, cbar=True)
            else:
              histlite.plot2d (ax, analysis.pdf_ratios_energy[year].epdfratios[g].h[-0.5].T,
              #histlite.plot2d (ax, analysis.pdfs_energy_sig[year].epdfs[g].h[-0.5].T / analysis.pdf_energy_bg[year].h[-0.5].T,
              cmap='RdBu_r', log=True, vmin=1e-2, vmax=1e2,cbar=True)
              ax.set_ylim(2,8.5)
              xlab = r'starting distance $[\mathrm{m}]$'
              ylab = r'$log(E[\mathrm{GeV}])$'
              ax.set_title('Energy Ratio (MC/EXP)')
            savingfig (fig, plot_dir, 'E-{}'.format(g))  
            plt.close(fig)
        '''

        #Let's plot out the skylab pdf.
        #breakpoint
        llhs = self.llhs
        for g in np.linspace(1,4,7):#analysis.gammas:
          for sample in analysis.datasets.keys():
            ps_llh = llhs[sample]
            plot_dir  = misc.ensure_dir ('{0}/plots/energypdfs/{1}'.format (self.mode_dir,sample))
            exp = ps_llh.exp
            spline = ps_llh.llh_model._w_spline_dict[(('gamma', g),)]
            saving (spline, '{0}/splinetests/{1}.spline'.format(self.mode_dir,sample))

            #SoB_skylab_energy = np.exp(ps_llh.llh_model._spline_eval(spline,exp))

            fig, ax = pfig()#plt.subplots (num=10)
            ax.cla()
            #spline_xbins=np.linspace(-1.,1.,40+1)
            #spline_ybins=np.logspace(1,10,28+1)
            xbins = np.linspace(-1.,1., np.shape(spline.values)[0])
            ybins = np.logspace(-1.,10., np.shape(spline.values)[1])
            #dummy_hist=hl.Hist((spline_xbins,spline_ybins),np.zeros((200,200)))
            #dense_hist2d=hl.hist_from_function((spline_xbins,spline_ybins),spline)
            #xmesh,ymesh=np.meshgrid(*dummy_hist.centers)

            #splinecounts=spline(xmesh,ymesh)
            #splinehist=hl.Hist(dummy_hist.bins,splinecounts.T)
            #exmaple:
            h = histlite.Hist((xbins,ybins), spline.values)
            #fig = plt.figure (figsize=(8,8))
            #ax = fig.add_subplot (111)
            #cz = np.linspace (-.99, .99, 100)
            #le = np.linspace (3.01, 7.99, 100)
            #CZ, LE = np.meshgrid (cz, le)
            #ax.pcolormesh (CZ, LE, h.spline_fit() (CZ, LE),
            #           cmap='Blues', vmin=0, vmax=180)

            histlite.plot2d (ax, h,
              cmap='RdBu_r', log=True, vmin=1e-2, vmax=1e2, cbar=True)
            ax.set_title('Energy Ratio (MC/EXP)')
            savingfig (fig, plot_dir, 'E-{}_skylab'.format(g))  
            plt.close(fig)

    @command
    def energypdf_compare (self):
        """
        Make '2D' energy pdf plots for csky and skylab - 4yrs in a grid.
        """
        analysis = self.analysis
        llhs = self.llhs
        #comparing the first four in a grid.
        analysis = self.analysis
        samples = self.analysis.datasets.keys()[0:4]
        for g in np.linspace(1,4,7):#analysis.gammas:
          fig, ax = pfig()
          ax.cla()
          gs = mpl.gridspec.GridSpec (2, 2, hspace=.3, wspace = .2)
          plot_dir  = misc.ensure_dir ('{0}/plots/energypdfs/compare4yr'.format (self.mode_dir))
          #csky - ratio
          for n,sample in zip(range(len(samples)),samples):
            ax = plt.subplot (gs[n])
            hsb = histlite.hist_from_function((np.linspace(-1.,1.,400,), np.linspace(1.,10.,400), [0,1]),
            self.analysis.pdf_ratios_energy[sample].epdfratios[g])[:,:,0]
            histlite.plot2d (ax, hsb, cmap='RdBu_r', log=True, vmin=1e-1, vmax=1e1, cbar=True)
            ax.set_xlim(-1.,1.)
            ax.set_ylim(1,10)
            xlab = r'$\sin{\delta}]$'
            ylab = r'$log(E[\mathrm{GeV}])$'
            ax.set_title('{}'.format(sample))
          plt.suptitle(r'Energy Ratio (MC/EXP) - $\gamma = ${}'.format(g))
          savingfig (fig, plot_dir, 'E-{}_csky'.format(g))  
          plt.close(fig)
          #csky - bg
          gs1 = mpl.gridspec.GridSpec (2, 2, hspace=.3, wspace = .2)
          fig, ax = pfig()
          ax.cla()
          for n,sample in zip(range(len(samples)),samples):
            ax = plt.subplot (gs1[n])
            hb = self.analysis.pdf_ratios_energy[sample].epdfratios[g].h_bg_0[:,:,0]
            #need a line for 68% containment?
            hb_c = hb.contain(1,frac=0.5) #1 for energy axis
            histlite.plot1d(ax, hb_c)
            histlite.plot2d (ax, hb, cmap='seismic', log=True, vmin=1e-1, vmax=1e1, cbar=True)
            ax.set_xlim(-1.,1.)
            ax.set_ylim(1,10)
            xlab = r'$\sin{\delta}]$'
            ylab = r'$log(E[\mathrm{GeV}])$'
            ax.set_title('{}'.format(sample))
          plt.suptitle(r'Energy Ratio (EXP) - $\gamma = ${}'.format(g))
          savingfig (fig, plot_dir, 'E-{}_csky_bg'.format(g))  
          plt.close(fig)
          
          #csky - sig
          gs2 = mpl.gridspec.GridSpec (2, 2, hspace=.3, wspace = .2)
          fig, ax = pfig()
          ax.cla()
          for n,sample in zip(range(len(samples)),samples):
            ax = plt.subplot (gs2[n])
            hb = self.analysis.pdf_ratios_energy[sample].epdfratios[g].h_sig_0[:,:,0]
            histlite.plot2d (ax, hb, cmap='OrRd_r', log=True, vmin=1e-1, vmax=1e1, cbar=True)
            ax.set_xlim(-1.,1.)
            ax.set_ylim(1,10)
            xlab = r'$\sin{\delta}]$'
            ylab = r'$log(E[\mathrm{GeV}])$'
            ax.set_title('{}'.format(sample))
          plt.suptitle(r'Energy Ratio (MC) - $\gamma = ${}'.format(g))
          savingfig (fig, plot_dir, 'E-{}_csky_sig'.format(g))  
          plt.close(fig)
          #skylab
          fig, ax = pfig()
          ax.cla()
          gs = mpl.gridspec.GridSpec (2, 2, hspace=.3, wspace = .2)
          for n,sample in zip(range(len(samples)),samples):
            ax = plt.subplot (gs[n])
            SoB = llhs[sample].llh_model._w_spline_dict[(('gamma', g),)]
            def hlSoB (sample, g):
              def f (sd, le, d):
                f = SoB
                out = f(np.transpose ([le.ravel(), sd.ravel()]))
                return 10**out.reshape (sd.shape)
              return f
            shsb = histlite.hist_from_function ((np.linspace (-1,1,400), np.linspace (1,10,400), [0, 1]),
            hlSoB(sample,g))[:,:,0]
            ax.set_xlim(-1.,1.)
            ax.set_ylim(1,10)
            histlite.plot2d (ax, shsb, cmap='RdBu_r', log=True, vmin=1e-1, vmax=1e1, cbar=True)
            xlab = r'$\sin{\delta}]$'
            ylab = r'$log(E[\mathrm{GeV}])$'
            ax.set_title('{}'.format(sample))
           
          plt.suptitle(r'Energy Ratio (MC/EXP) - $\gamma = ${}'.format(g))
          savingfig (fig, plot_dir, 'E-{}_skylab'.format(g))  
          plt.close(fig)

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


    # background and signal trials

    @command
    def submit_do_bg_trials (self):
        """Submit bg TSDist jobs to cluster."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_bg_trials/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))


        if self.opts.just_diff:
            dec_degs = [-60]
        else:
            #dec_degs = [16]
            dec_degs = np.arange (-90, 90, 3) #but, let's try narrowing in on the horizon.
            dec_degs = np.concatenate([np.linspace(-90.,-6.,29), np.linspace(-5.,5.,11), np.linspace(6,90,29)])


        #for test_ext in np.arange(11):
        for test_ext in [0.]:
            for dec_deg in dec_degs:
                dec_deg = max (-89, min (89, dec_deg))
                for i_job in xrange (int (self.opts.n_jobs)):
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
    def submit_do_bg_trials_stacking (self):
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
            raise ValueError('Must provide a weighting scheme for a stacking analysis')

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
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
            if self.opts.weights == 'equal':
                print ("No weighting given - use equal weights")
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
        else:
            dec_deg = self.opts.test_dec
            test_ext_deg = self.opts.test_ext
            ra_deg = 0 
            weights = 1. #to avoid redundancy in establishing a ps hypothesis

        test_ext = np.radians(test_ext_deg) #to allow float or list
        dec = np.radians(dec_deg)
        ra = np.radians(ra_deg)

        analysis=self.analysis

        #analysis.set_energy_pdfs (disable=True) #Include to test space-only distributions.
        np.random.seed (seed)
        data = hyp.DataHypothesis (analysis)
        gammas = np.linspace(1.,4,25)

        ps = hyp.PointSourceHypothesis (analysis, dec, ra, 2, extensions = test_ext, weights = weights, sigsub = self.conf.root.sigsub)

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
                raise ValueError('No weighting scheme chosen for specified catalog')
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
            tsd = bg_tsds#[self.opts.cat][self.conf.root.n_year]
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
        Submit n_sig jobs for some n-sigma and beta values.
        """
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_n_sig/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        spectra = [(2, np.inf)]
        i_job = 0
        if self.opts.do_ext:
            src_exts = np.arange (0, 10.5)
        else:
            src_exts = [0]
            #src_exts = np.arange(6)
        betas = [0.9,0.5]
        sigmas = [0,5]
        for sigma, beta in zip(sigmas,betas):
          for src_ext in src_exts:
            test_exts = [src_ext]# if src_ext == 0 else [0, src_ext]
            #test_exts = [0]
            for test_ext in test_exts:
                for (src_gamma, src_cutoff) in spectra:
                    if self.opts.just_diff:
                        dec_degs = [150]
                    else:
                        #dec_degs = [16]
                        dec_degs = np.arange (-90, 90, 3)
                        dec_degs = np.concatenate([np.linspace(-90.,-6.,29), np.linspace(-5.,5.,11), np.linspace(6,90,29)])
                    for dec_deg in dec_degs:
                        dec_deg = max (-89, min (89, dec_deg))
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
                                    sigma,
                                    beta,
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
    def submit_do_n_sig_stacking (self):
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
        spectra = [(2, np.inf)]
        i_job = 0
        cat = self.opts.cat
        betas = [0.9,0.5]
        sigmas = [0,5]

        #Weighting
        weights = self.opts.weights
        for sigma, beta in zip(sigmas,betas):
            if not weights:
                raise ValueError('Must provide a weighting scheme for a stacking analysis')

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
                            sigma,
                            beta,
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
            if self.opts.weights:
                if self.opts.weights == 'equal':
                    print ("No weighting given - use equal weights")
                    weights = np.ones_like(params['dec'])
                else:
                    weights = params[str(self.opts.weights)]
            else:
                raise ValueError('Must specify a weighting scheme for stacking')
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            test_ext_deg = [ext for ext in test_ext]
            src_ext_deg = [ext for ext in src_ext]
        else:
            dec_deg = self.opts.test_dec
            test_ext_deg = self.opts.test_ext
            ra_deg = 0
            weights = 1. #put in to eliminate redundancy for ps injection 
            src_ext_deg = self.opts.src_ext
        src_ext = np.radians(src_ext_deg)
        test_ext = np.radians(test_ext_deg)
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
          self.analysis.ps_sindec_width = .02
          bg_tsd = self.bg_tsds
        else:
          self.analysis.ps_sindec_width  = .2 if dec_deg < 60 else 0.05
          bg_tsd = bk.get_best (self.bg_tsds, test_ext_deg, dec_deg)

        ts = bg_tsd.isf (stats.norm.sf (sigma) + eps, fit=use_fit)

        prush ('Setting up injection...')
        data = hyp.DataHypothesis (self.analysis)
        gammas = np.linspace(1.,4,25)
        ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, gamma, weights = weights,
                                        extensions = src_ext, energy_range=(10**thresh, 10**cutoff),
                                        sigsub = self.conf.root.sigsub)
        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas, extensions=test_ext))

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
    def do_one_TS (self):
        """
        Do a single unscrambled test for TS, fitted ns.
        """

        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        n_year = self.conf.root.n_year
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
        if not self.opts.weights:
          self.opts.weights = 'equal'
          weights = np.ones_like(params['dec'])
          print ("Weights = equal")
        else:
          weights = params[str(self.opts.weights)]
        if 'extension' in params.keys():
            test_ext = params['extension']
        else:
            test_ext = np.zeros_like(params['dec'])
        dec_deg = np.degrees(params['dec'])
        ra_deg = np.degrees(params['ra'])
        test_ext_deg = [ext for ext in test_ext]
        test_ext = np.radians(test_ext_deg)
        dec = np.radians(dec_deg)
        ra = np.radians(ra_deg)

        if self.opts.cat or self.conf.root.n_year==1:
          bg_tsd = self.bg_tsds
        else:
          bg_tsd = bk.get_best (self.bg_tsds)

        np.random.seed (seed)
        data = hyp.DataHypothesis (self.analysis)
        gammas = self.analysis.pdf_ratios_energy['IC40'].gammas
        ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, 2, extensions = test_ext, weights = weights, sigsub = self.conf.root.sigsub)

        tr = trial.PSTrialRunner (data, ps, ps.tests (gammas), cut_n_sigma=np.inf)

        c = tr.get_one_ts_ns_gamma (TRUTH=True)

        pvalue = bg_tsd.sf(c[0])
        
        result_dict = {'pvalue':pvalue,'TS':c[0],'ns':c[1],'gamma':c[2]}

        print(result_dict)
        sm = bk.SavingModel (
                'ts',
                'one_ts.dict',
                )
        filename = sm.save (
                 result_dict, '{0}/cats/{1}/{2}/one_ts_truth'.format (self.mode_dir,self.opts.cat,self.opts.weights))

        prush ('->', filename)

    @command
    def get_all_pdfs (self):
        """
        Do a single unscrambled test, get assoc. pdfs (bg and sig, energy and space)
        """

        seed = self.opts.seed
        n_year = 1

        data = hyp.DataHypothesis (self.analysis)
        ps = hyp.PointSourceHypothesis (self.analysis, 0, 0, 2, sigsub = False)

        tr = trial.PSTrialRunner (data, ps, ps.tests (self.analysis.gammas))

        calc = tr.get_one_calc (TRUTH=True)
        L_2_00 = calc.Ls[np.argmin (np.abs (2 - self.analysis.gammas))]

        sm = bk.SavingModel (
                'ts',
                'one_ts.dict',
                )
        filename = sm.save (
                 result_dict, '{0}/cats/{1}/{2}/one_ts_truth'.format (self.mode_dir,cat,weight))

        prush ('->', filename)

    @command
    def collect_sig_int (self):
        """Collect signal injection results and cache."""
        prush ('Collecting signal results...')
        if self.opts.cat:
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

        return d1

    @command
    def sig_plots (self):
        """
        Make plots for TS distribution for injected trials. If allsky, also plot eta and dof across sky.
        """
        if self.opts.cat:
            bg_tsd = loading ('{0}/cats/{1}/{2}/bg_tsds.dict'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            sens_ts = loading ('{0}/cats/{1}/{2}/sig_int/0/0.9/2.00/2.00/+000inf/tsdist.pickle'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            disc_ts = loading ('{0}/cats/{1}/{2}/sig_int/5/0.5/2.00/2.00/+000inf/tsdist.pickle'.format (self.mode_dir,self.opts.cat,self.opts.weights))
            plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}/plots/sig_int'.format (self.mode_dir,self.opts.cat,self.opts.weights))
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
            plt.close(fig)
        else:
            bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))[0]
            plot_dir  = misc.ensure_dir ('{0}/plots/bg_tsds'.format (self.mode_dir))
            sindec = np.sin (np.radians (sorted (bg_tsds)))

            sensdir = self.mode_dir+'/sig_int/'
            sig_int = bk.get_all(sensdir,'tsdist.pickle',lambda x:x)
            xs = sig_int[0][.9][0][0][2][2][np.inf]
            xd = sig_int[5][.5][0][0][2][2][np.inf]
             
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
              sens_ts = xs[dec][0] ; disc_ts = xd[dec][0]
              h_sens = histlite.hist(sens_ts,bins=100).normalize()
              h_disc = histlite.hist(disc_ts,bins=100).normalize()
              histlite.plot1d (ax, h, label = '{} background trials'.format(tsd.n_total)) ; plt.plot (x,tsd.eta*tsd.chi2.pdf(x), color= 'k', ls='--',
                            label = r'$\tilde{\chi}^2 (n_{dof} = $' + str(np.round(tsd.chi2.args[0],2)) + r', $\eta = $' + str(np.round(tsd.eta,2)) + ')')
              histlite.plot1d (ax, h_sens, color = 'red', label = '{} sens injection trials'.format(len(sens_ts)))
              histlite.plot1d (ax, h_disc, color = 'green', label = '{} disc injection trials'.format(len(disc_ts)))
              ax.semilogy()
              ax.set_title("TS distribution - {0}yr - dec {1}".format(self.conf.root.n_year, dec))
              ax.set_ylim (1e-4, 10)
              ax.set_xlim (0, 40)
              plt.legend(loc='upper right', prop=propsmall)
              ax.set_xlabel('TS')
              savingfig (fig, plot_dir, 'tsdist')  
              plt.close(fig)

    @command
    def plot_cat (self):
        """Plot a catalog of sources on a healpy skymap."""
        cat = self.opts.cat
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
        weights = self.opts.weights
        if weights:
          weight = params[weights]
        else:
          weight = np.ones_like(params['dec'])
          weights = 'equal'
        if not cat:
          raise ValueError("Catalog must be specified.")
        #catalog info:        
        weightplot=[(i*100)/max(weight) for i in weight]
        ra,dec = np.degrees(params['ra']),np.degrees(params['dec'])
        
        plot_dir  = misc.ensure_dir (self.root_dir+'/catalogs/{}'.format(cat))
        fig = plt.plot()
        #Rot = 180 shifts the axis from the center to the edge
        #centered by default
        healpy.mollview(title = 'Equatorial map of {}'.format(cat), cbar = False,
                    rot = 180, notext = True, cmap=None, coord='C')
        healpy.graticule(coord='C', color='DimGrey')
        py.title("{} - {} - Equatorial".format(cat,weights), fontsize = 25, fontweight='bold')
        #healpy.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
        healpy.projscatter(ra,dec,coord='C',lonlat=True, s = weightplot)
        healpy.projtext(185, 2,'180', lonlat=True, fontweight = 'bold')
        healpy.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
        healpy.projtext(275, 2, '270', lonlat=True, fontweight = 'bold')
        healpy.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
        healpy.projtext(359, 2, '360', lonlat=True, fontweight = 'bold')
        healpy.projtext(193, -8, 'RA (deg)', lonlat=True, fontweight = 'bold')
        healpy.projtext(350, 30, '30', lonlat=True, fontweight = 'bold')
        healpy.projtext(340, 60, '60', lonlat=True, fontweight = 'bold')
        #healpy.projtext(5, -5, 'Dec (deg)', lonlat=True)
        healpy.projtext(358, -33.5, '-30', lonlat=True, fontweight = 'bold')
        healpy.projtext(358, -63.5, '-60', lonlat=True, fontweight = 'bold')
        ra_gplane = np.arange(0.,361.,1.)
        dec_gplane = np.zeros(len(ra_gplane))
        
        healpy.projplot(ra_gplane, dec_gplane, coord='G', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)
        #savingfig (fig, plot_dir, '{0}_{1}_equatorial'.format(cat,self.opts.weights)) 
        fname = '/{0}_{1}_equatorial.png'.format(cat,weights)
        py.savefig(plot_dir+fname)
        #py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/equatorialmap.png')

        public_dir = '/home/brelethford/public_html/csky/' + plot_dir.split('brelethford/')[-1]
        print (public_dir+fname)
        subprocess.call(['ssh', 'cobalt06', 'mkdir', '-p', '{0}/'.format(public_dir)])
        subprocess.call(["scp", '-r','{0}/'.format (plot_dir), 'brelethford@cobalt06:{}'.format(public_dir)])
        '''
        fig_gal, ax_gal = pfig()
        #Galactic plane requires no rotation to put it in the form we want.
        healpy.mollview(title = 'Galactic Map of SwiftBAT AGNs', cbar = False,
                    cmap = None, notext = 'True', coord = 'G')
        healpy.graticule(coord='G',color='DimGrey')
        plt.title("{} {} - Galactic".format(cat,weight), fontsize = 25, fontweight = 'bold')
        #healpy.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
        healpy.projscatter(ra,dec, coord='C', lonlat=True, s = weightplot)
        healpy.projtext(175, 2,'180', lonlat=True, fontweight = 'bold')
        healpy.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
        healpy.projtext(275, 2, '-90', lonlat=True, fontweight = 'bold')
        healpy.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
        healpy.projtext(-160, 2, '-180', lonlat=True, fontweight = 'bold')
        healpy.projtext(0, -8, 'Latitude (deg)', lonlat=True, fontweight = 'bold')
        healpy.projtext(175, 30, '30', lonlat=True, fontweight = 'bold')
        healpy.projtext(165, 60, '60', lonlat=True, fontweight = 'bold')
        #healpy.projtext(5, -5, 'Longitude (deg)', lonlat=True)
        healpy.projtext(178, -33.5, '-30', lonlat=True, fontweight = 'bold')
        healpy.projtext(175, -63.5, '-60', lonlat=True, fontweight = 'bold')
        ra_cplane = np.arange(0.,361.,1.)
        dec_cplane = np.zeros(len(ra_cplane))
        healpy.projplot(ra_cplane, dec_cplane, coord='C', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)

        #save        

        py.savefig('/home/relethford/Documents/IceCube_Research/Plots/AGNCore/X-Ray_Catalogue/SwiftBAT_70M/galacticmap.png')
        savingfig (fig, plot_dir, '{0}_{1}'.format(cat,weight)) 
        
        print ( 'number of sources = ' + str(len(ra)))


        prush ('->', filename)
        '''

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
    
    @command
    def sens_ext_onedec (self):
        """Plot extended sensitivity and ."""

        dec = 16.0

        sig_int = self.sig_int

        misc.tex_mpl_rc (False)

        fig = plt.figure(100)
        fig.clf()
        gs = mpl.gridspec.GridSpec (1, 1)
        ax = plt.subplot (gs[0])
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
        label_ext = 'Extended Source llh'
        label_ps = 'Point Source llh'
        sensdir = self.mode_dir+'/sig_int/'
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        exts = np.arange(6)
        #x = [sig_int[0][.9][i][i][2][2][np.inf][16] for i in exts]
        #sens2 = np.array ([
        #    x[k][-1] for k in sorted (x)])
        #mask2 = np.r_[False,sens2[1:] / sens2[:-1] > .5]
        #ax.semilogy (exts, sens2,
        #        label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #disc - gamma=2
        x = [sig_int[5][.5][i][i][2][2][np.inf][16] for i in exts]
        y = [sig_int[5][.5][0][i][2][2][np.inf][16] for i in exts]
        discx = np.array ([
            k[-1] for k in sorted (x)])
        discy = np.array ([
            k[-1] for k in sorted (y)])
        ax.semilogy (exts, discx,
                label=label_ext, lw=lw, ls='--', color='grey')
        ax.semilogy (exts, discy,
                label=label_ps, lw=lw, ls=ls, color='k')
        ax.set_ylim (10**-12, 10**-9.)
        
        #ax.set_yticks (np.logspace (-12, -9, 7))
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=2, handlelength=2.2,
        )

        ax.grid ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=0, y=1.2 * 10**-9., ha='left', va='bottom')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        savingfig (fig, plot_dir, 'sens_ext_dec_{}'.format(dec))

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
        label2 = 'csky / $E^{{-{0}}}${1}'.format (
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
        ax.semilogy (sindec2, sens2,
                label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #ax.semilogy (sindec2[mask2], sens2[mask2],
        #        label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #sens - gamma = 3
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label3 = 'csky / $E^{{-{0}}}${1}'.format (
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
        #ax.semilogy (sindec3[mask3], sens3[mask3],
        #        label=label3 + ' - sens', lw=lw, ls=ls, color=colors[2])
        ax.semilogy (sindec3, sens3,
                label=label3 + ' - sens', lw=lw, ls=ls, color=colors[2])
        #disc - gamma=2
        x = sig_int[5][.5][0][0][2][2][np.inf]
        sindec2_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        mask2d = np.r_[False,disc2[1:] / disc2[:-1] > .5]
        #ax.semilogy (sindec2_disc[mask2d], disc2[mask2d],
        #        label=label2 + ' - disc', lw=lw, ls=ls, color=colors[1])
        ax.semilogy (sindec2_disc, disc2,
                label=label2 + ' - disc', lw=lw, ls=ls, color=colors[1])
        #disc - gamma=
        x = sig_int[5][.5][0][0][3][2][np.inf]
        sindec3_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc3 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        mask3d = np.r_[False,disc3[1:] / disc3[:-1] > .5]
        #ax.semilogy (sindec3_disc[mask3d], disc3[mask3d],
        #        label=label3 + ' - disc', lw=lw, ls=ls, color=colors[3])
        ax.semilogy (sindec3_disc, disc3,
                label=label3 + ' - disc', lw=lw, ls=ls, color=colors[3])
        #ax.semilogy (sindec_w, sens_w,
        #        lw=1, alpha=.5, color=color)

        #Stefan's results
        lw = 2
        thing_sens = np.load ('{0}/data/track_sens.npy'.format (self.root_dir))
        thing_disc = np.load ('{0}/data/track_disc.npy'.format (self.root_dir))
        sd, s2, s3 = thing_sens['dec'], sys*thing_sens['2'], sys*thing_sens['3']
        ax.semilogy (np.sin (sd), s2,
                   color=colors[0], lw=lw, ls=':', alpha = 0.9,
                   label=r'skylab / $E^{-2}$ - sens')
        ax.semilogy (np.sin (sd), s3 / 100,
                   color=colors[2], lw=lw, ls=':', alpha = 0.9,
                   label=r'skylab / $E^{-3}$ - sens')

        dd, d2, d3 = thing_disc['dec'], sys*thing_disc['2'], sys*thing_disc['3']
        ax.semilogy (np.sin (dd), d2,
                   color=colors[1], lw=lw, ls=':', alpha = 0.9,
                   label=r'skylab / $E^{-2}$ - disc')
        ax.semilogy (np.sin (dd), d3 / 100,
                   color=colors[3], lw=lw, ls=':', alpha = 0.9,
                   label=r'skylab / $E^{-3}$ - disc')
        #ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13, 10**-9.)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-13, -9, 5))
        #plt.yticks (np.logspace (-13, -7, 7))

        ax.text (.95, 10**-12.9, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-12.9, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=2, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        )
        #leg.get_frame ().set_linewidth (0)

        #some Ratio plots
        rat2 = np.interp (sd, np.arcsin(sindec2), sens2) / s2 
        rat3 = np.interp (sd, np.arcsin(sindec3), sens3) / (s3 / 100)
        rat2disc = np.interp (dd, np.arcsin(sindec2_disc), disc2) / d2
        rat3disc = np.interp (dd, np.arcsin(sindec3_disc), disc3) / (d3 / 100)
        #rat2 = np.interp (sd, np.arcsin(sindec2[mask2]), sens2[mask2]) / s2 
        #rat3 = np.interp (sd, np.arcsin(sindec3[mask3]), sens3[mask3]) / (s3 / 100)
        #rat2disc = np.interp (dd, np.arcsin(sindec2_disc[mask2d]), disc2[mask2d]) / d2
        #rat3disc = np.interp (dd, np.arcsin(sindec3_disc[mask3d]), disc3[mask3d]) / (d3 / 100)
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
        savingfig (fig, plot_dir, 'npz_vs_track_sensitivity_sq')

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
        ps = hyp.PointSourceHypothesis (self.analysis, multi_dec, multi_ra, 2, sigsub=False)
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
                                        energy_range=(10**thresh, 10**cutoff), sigsub=False)
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

        sens = tr.get_sens (ts, beta, n_batch=200, tol=.02,full_output=True)
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
        filename_ts = sm.save (
                sens['tss'],
                '{0}/tests/{1}/{2}/sig_int/'.format (self.mode_dir,test,n_multi),
                sigma, beta, test_ext_deg, src_ext_deg,
                gamma, thresh, cutoff, dec_deg,
                outname_ts)
        prush ('->', filename)
        prush ('->', filename_ts)

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
    def angres_kent (self):
        """Make Kent angular resolution plots."""
        misc.tex_mpl_rc (True)
        plt.rc ('font', size=16)
        apr = self.analysis.apr['IC86']
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
        savingfig (fig, plot_dir, 'my7yr_stacking_tests')


    @command
    def sens (self):
        """Plot integrated sensitivity."""

        #sig_int = self.sig_int
        sig_int = self.sig_int

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
        #sens - gamma = 
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label2 = 'csky / $E^{{-{0}}}${1}'.format (
                #int (self.conf.root.n_year), gamma, cutoff_str)
                2, cutoff_str)
        sensdir = self.mode_dir + '/sig_int/'
        #sig_int = bk.get_all(
        #        sensdir,'sens.pickle', lambda x:x[0])
        ext = self.opts.test_ext
        print('Using extension of {}'.format(ext))
        x = sig_int[0][.9][ext][0][2][2][np.inf]
        sindec2 = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        sys = self.opts.sys_total
        sens2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec2, sens2,
                label=label2 + ' - sens', lw=lw, ls=ls, color=colors[0])
        #disc - gamma=2

        x = sig_int[5][.5][ext][0][2][2][np.inf]
        sindec2_disc = np.array ([
            np.sin (np.radians(k)) for k in sorted (x)])
        disc2 = sys*np.array ([
            x[k][-1] for k in sorted (x)])
        ax.semilogy (sindec2_disc, disc2,
                label=label2 + ' - disc', lw=lw, ls=ls, color=colors[1])
        #Stefan's results
        if 'oneyear' in self.mode.words:
          thing_sens = np.genfromtxt('{0}/data/sens_86.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/disc_86.csv'.format (self.root_dir), delimiter = ',')
        elif 'IC86II' in self.mode.words or 'IC86II_corrected' in self.mode.words:
          thing_sens = np.genfromtxt('{0}/data/IC86II/IC86II_sens.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/IC86II/IC86II_disc.csv'.format (self.root_dir), delimiter = ',')
        elif 'yrs4' in self.mode.words[0]:
          thing_sens = np.genfromtxt('{0}/data/sens4yr.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/disc4yr.csv'.format (self.root_dir), delimiter = ',')
        #elif 'yrs7' or 'SNR' in self.mode.words:
        else:
          thing_sens = np.load ('{0}/data/track_sens.npy'.format (self.root_dir))
          thing_disc = np.load ('{0}/data/track_disc.npy'.format (self.root_dir))
        lw = 2
        if 'yrs7' in self.mode.words[0] or 'SNR' in self.mode.words:
          sys = self.opts.sys_total
          sd, s2, s3 = np.sin(thing_sens['dec']), sys*thing_sens['2'], sys*thing_sens['3']
          dd, d2, d3 = np.sin(thing_disc['dec']), sys*thing_disc['2'], sys*thing_disc['3']
        else:
          sd, s2 = thing_sens.T
          dd, d2 = thing_disc.T
        ax.semilogy (sd, s2,
                   color=colors[0], lw=lw, ls=':', alpha = 0.5,
                   label=r'official results / $E^{-2}$ - sens')

        ax.semilogy (dd, d2,
                   color=colors[1], lw=lw, ls=':', alpha = 0.5,
                   label=r'official results / $E^{-2}$ - disc')
        #ax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (10**-13., 10**-9.)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-13, -9, 5))
        ## optional catalog for src-specific weights
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
            if 'extension' in params.keys():
                test_ext = params['extension']
                src_ext = params['extension']
            if self.opts.weights:
                try:
                  weights = params[str(self.opts.weights)]
                except:
                  raise ValueError('{0} not a parameter in {1}'.format(self.opts.weights,self.opts.cat))
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            dec = dec_deg / 180.*pi
            ra = ra_deg / 180.*pi
            gamma = self.opts.src_gamma
            thresh = self.opts.src_thresh
            cutoff = self.opts.src_cutoff

            # get the params
            data = hyp.DataHypothesis (self.analysis)
            gammas = np.linspace(1.,4,25)
            ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, gamma, weights = weights,
                                            energy_range=(10**thresh, 10**cutoff),
                                            sigsub = False)#self.conf.root.sigsub)
            acc_w = ps.p_sources/np.sum(ps.p_sources)
            theo_w = weights/np.sum(weights)
            #load the sensdisc and divvy it up per source in each way.
            sensdir = self.mode_dir + '/cats/{0}/{1}/'.format(self.opts.cat,self.opts.weights)
            sig_int_cat = cache.load(sensdir+'sig_int.dict')
            sens_cat = sig_int_cat[0][.9][2][2]['+000inf'][2]
            disc_cat = sig_int_cat[5][.5][2][2]['+000inf'][2]
            #plot!
            ax.scatter (np.sin(params['dec']), sens_cat*theo_w,
                label='sens per source', marker = 'o',color='blue')
            ax.scatter (np.sin(params['dec']), disc_cat*theo_w,
                label='disc per source', marker = 'o',color='red')
            #also plot the sens across as a hline
            ax.hlines(sens_cat,-1.,1.,colors='blue',alpha=.5)
            ax.hlines(disc_cat,-1.,1.,colors='red',alpha=.5)
            ax.set_title('{0} - per source contribution'.format(self.opts.cat))
        #plt.yticks (np.logspace (-13, -7, 7))

        ax.text (.95, 10**-12.9, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-12.9, 'South', size='small', ha='left', va='bottom')
        leg = ax.legend (
            loc='upper right',
            prop=propsmall, ncol=2, handlelength=2.2,
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
        if 'oneyear' or 'IC86II' or 'SNR' in self.mode.words:
            rax.set_ylabel ('Drexel / Stefan')
        elif 'yrs4' in self.mode.words:
            rax.set_ylabel ('Drexel / Jake')
              #https://wiki.icecube.wisc.edu/index.php/IC86_I_Point_Source_Analysis/Performance_Plots
        rax.set_xlim (-1, 1)
        rax.set_ylim (0.5,1.5)

        ax.grid ()
        rax.grid ()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')
        if self.opts.cat:
          plot_dir = misc.ensure_dir ('{0}/{1}/{2}/plots'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        else:
          plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        if ext !=0.0:
          savingfig (fig, plot_dir, '{}_sensitivity_ext_{}'.format(self.mode.words[0],ext))
        else:
          savingfig (fig, plot_dir, '{}_sensitivity'.format(self.mode.words[0]))

    @command
    def sens_ext (self):
        """Plot integrated sensitivity."""

        #sig_int = self.sig_int
        sig_int = self.sig_int

        misc.tex_mpl_rc (False)

        fig = plt.figure(100)
        fig.clf()
        gs = mpl.gridspec.GridSpec (1, 1)#, height_ratios=[3,1], hspace=.15)
        ax = plt.subplot (gs[0])
        #rax = plt.subplot (gs[1], sharex=ax)
        colors = soft_colors
        lw = 2
        cutoff = np.inf
        if cutoff == 5:
            ls = '--'
        else:
            ls = '-'
        #Stefan's results
        if 'oneyear' in self.mode.words:
          thing_sens = np.genfromtxt('{0}/data/sens_86.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/disc_86.csv'.format (self.root_dir), delimiter = ',')
        elif 'IC86II' in self.mode.words or 'IC86II_corrected' in self.mode.words:
          thing_sens = np.genfromtxt('{0}/data/IC86II/IC86II_sens.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/IC86II/IC86II_disc.csv'.format (self.root_dir), delimiter = ',')
        elif 'yrs4' in self.mode.words:
          thing_sens = np.genfromtxt('{0}/data/sens4yr.csv'.format (self.root_dir), delimiter = ',')
          thing_disc = np.genfromtxt('{0}/data/disc4yr.csv'.format (self.root_dir), delimiter = ',')
        #elif 'yrs7' or 'SNR' in self.mode.words:
        else:
          thing_sens = np.load ('{0}/data/track_sens.npy'.format (self.root_dir))
          thing_disc = np.load ('{0}/data/track_disc.npy'.format (self.root_dir))
        lw = 2
        if 'yrs7' in self.mode.words[0] or 'SNR' in self.mode.words:
          sys = self.opts.sys_total
          sd, s2, s3 = np.sin(thing_sens['dec']), sys*thing_sens['2'], sys*thing_sens['3']
          dd, d2, d3 = np.sin(thing_disc['dec']), sys*thing_disc['2'], sys*thing_disc['3']
        else:
          sd, s2 = thing_sens.T
          dd, d2 = thing_disc.T
        ax.semilogy (sd, s2,
                   color=colors[0], lw=lw, ls=':', alpha = 0.5,
                   label=r'official results / $E^{-2}$ - sens')

        ## optional catalog for src-specific weights
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
            if 'extension' in params.keys():
                test_ext = params['extension']
                src_ext = params['extension']
            if self.opts.weights:
                try:
                  weights = params[str(self.opts.weights)]
                except:
                  raise ValueError('{0} not a parameter in {1}'.format(self.opts.weights,self.opts.cat))
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            dec = dec_deg / 180.*pi
            ra = ra_deg / 180.*pi
            gamma = self.opts.src_gamma
            thresh = self.opts.src_thresh
            cutoff = self.opts.src_cutoff

            # get the params
            data = hyp.DataHypothesis (self.analysis)
            gammas = np.linspace(1.,4,25)
            ps = hyp.PointSourceHypothesis (self.analysis, dec, ra, gamma, weights = weights,
                                            energy_range=(10**thresh, 10**cutoff),
                                            sigsub = False)#self.conf.root.sigsub)
            acc_w = ps.p_sources/np.sum(ps.p_sources)
            theo_w = weights/np.sum(weights)
            #load the sensdisc and divvy it up per source in each way.
            sensdir = self.mode_dir + '/cats/{0}/{1}/'.format(self.opts.cat,self.opts.weights)
            sig_int_cat = cache.load(sensdir+'sig_int.dict')
            sens_cat = sig_int_cat[0][.9][2][2]['+000inf'][2]
            disc_cat = sig_int_cat[5][.5][2][2]['+000inf'][2]
            #plot!
            ax.scatter (np.sin(params['dec']), sens_cat*theo_w,
                label='sens per source', marker = 'o',color='blue')
            ax.scatter (np.sin(params['dec']), disc_cat*theo_w,
                label='disc per source', marker = 'o',color='red')
            #also plot the sens across as a hline
            ax.hlines(sens_cat,-1.,1.,colors='blue',alpha=.5)
            ax.hlines(disc_cat,-1.,1.,colors='red',alpha=.5)
            ax.set_title('{0} - per source contribution'.format(self.opts.cat))
        #plt.yticks (np.logspace (-13, -7, 7))
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        ax.set_ylim (3*10**-13, 9*10**-11)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-12, -10, 5))

        ax.text (.95, 10**-12.9, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-12.9, 'South', size='small', ha='left', va='bottom')


        #sens - gamma = 
        cutoff_str = '' if cutoff == np.inf \
                else ' cutoff at ${0:.0f}$ TeV'.format (10**cutoff/1e3)
        label2 = 'csky / $E^{{-{0}}}${1}'.format (
                #int (self.conf.root.n_year), gamma, cutoff_str)
                2, cutoff_str)
        sensdir = self.mode_dir + '/sig_int/'
        sig_int = bk.get_all(
                sensdir,'sens.pickle', lambda x:x[0])
        ext = self.opts.test_ext
        print('Using extension of {}'.format(ext))
        colors = ['k','b','g','r','c','m','y','chartreuse','burlywood','brown','violet']
        for ext,color in zip(np.arange(11),colors):
          x = sig_int[0][.9][ext][ext][2][2][np.inf]
          sindec2 = np.array ([
              np.sin (np.radians(k)) for k in sorted (x)])
          sens2 = np.array ([
              x[k][-1] for k in sorted (x)])
          ax.semilogy (sindec2, sens2,
                  label=label2 + ' - sens - ext {}'.format(ext), lw=lw, ls=ls, color=color)
          #disc - gamma=2
          # x = sig_int[5][.5][ext][0][2][2][np.inf]
          # sindec2_disc = np.array ([
          #   np.sin (np.radians(k)) for k in sorted (x)])
          # disc2 = sys*np.array ([
          #   x[k][-1] for k in sorted (x)])
          # ax.semilogy (sindec2_disc, disc2,
          #       label=label2 + ' - disc', lw=lw, ls=ls, color=colors)


          #some Ratio plots
          rat2 = np.interp (sd, sindec2, sens2) / s2 
          #rat2disc = np.interp (dd, sindec2_disc, disc2) / d2
          #rax.plot (sd, rat2, color=color)
          #rax.plot (dd, rat2disc, color=colors[1])
        
        #leg = ax.legend (
        #    loc='upper right',
        #    prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        #)
        #rax.set_xlabel (r'$\sin(\delta)$')
        #if 'oneyear' or 'IC86II' or 'SNR' in self.mode.words:
        #    rax.set_ylabel ('Drexel / Stefan')
        #elif 'yrs4' in self.mode.words:
        #    rax.set_ylabel ('Drexel / Jake')
              #https://wiki.icecube.wisc.edu/index.php/IC86_I_Point_Source_Analysis/Performance_Plots
        #rax.set_xlim (-1, 1)
        #rax.set_ylim (0.5,1.5)

        ax.grid ()
        #rax.grid ()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')
        if self.opts.cat:
          plot_dir = misc.ensure_dir ('{0}/{1}/{2}/plots'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        else:
          plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        if ext !=0.0:
          savingfig (fig, plot_dir, '{}_sensitivity_ext'.format(self.mode.words[0],ext))
        else:
          savingfig (fig, plot_dir, '{}_sensitivity'.format(self.mode.words[0]))

        #now disc
        fig = plt.figure(100)
        fig.clf()
        gs = mpl.gridspec.GridSpec (1, 1)#, height_ratios=[3,1], hspace=.15)
        ax = plt.subplot (gs[0])
        #rax = plt.subplot (gs[1], sharex=ax)

        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')

        ax.set_ylim (10**-12, 2*10**-10)
        ymin, ymax = ax.get_ylim ()
        ax.set_yticks (np.logspace (-12, -10, 5))
        ax.semilogy (dd, d2,
                   color=colors[1], lw=lw, ls=':', alpha = 0.5,
                   label=r'official results / $E^{-2}$ - disc')
        ax.set_xlabel (r'$\sin(\delta)$')

        ax.text (.95, 10**-11.9, 'North', size='small', ha='right', va='bottom')
        ax.text (-.95, 10**-11.9, 'South', size='small', ha='left', va='bottom')
        for ext,color in zip(np.arange(11),colors):
          x = sig_int[5][.5][ext][ext][2][2][np.inf]
          sindec2 = np.array ([
              np.sin (np.radians(k)) for k in sorted (x)])
          sens2 = np.array ([
              x[k][-1] for k in sorted (x)])
          ax.semilogy (sindec2, sens2,
                  label=label2 + ' - disc - ext {}'.format(ext), lw=lw, ls=ls, color=color)
          #disc - gamma=2
          # x = sig_int[5][.5][ext][0][2][2][np.inf]
          # sindec2_disc = np.array ([
          #   np.sin (np.radians(k)) for k in sorted (x)])
          # disc2 = sys*np.array ([
          #   x[k][-1] for k in sorted (x)])
          # ax.semilogy (sindec2_disc, disc2,
          #       label=label2 + ' - disc', lw=lw, ls=ls, color=colors)


          #some Ratio plots
          rat2 = np.interp (dd, sindec2, sens2) / d2 
          #rat2disc = np.interp (dd, sindec2_disc, disc2) / d2
          #rax.plot (sd, rat2, color=color)
          #rax.plot (dd, rat2disc, color=colors[1])
        
        #leg = ax.legend (
        #    loc='upper right',
        #    prop=propsmall, ncol=1, handlelength=2.2,
            #bbox_to_anchor=(0.02, .9, .96, .10),
            #mode='expand'
        #)
        #rax.set_xlabel (r'$\sin(\delta)$')
        #if 'oneyear' or 'IC86II' or 'SNR' in self.mode.words:
        #    rax.set_ylabel ('Drexel / Stefan')
        #elif 'yrs4' in self.mode.words:
        #    rax.set_ylabel ('Drexel / Jake')
        #      #https://wiki.icecube.wisc.edu/index.php/IC86_I_Point_Source_Analysis/Performance_Plots
        #rax.set_xlim (-1, 1)
        #rax.set_ylim (0.5,1.5)

        ax.grid ()
        #rax.grid ()
        #fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        #plt.tight_layout ()
        fig.subplots_adjust (top=.93)
        icprelim (ax, x=-1, y=1.2 * ymax, ha='left', va='bottom')
        if self.opts.cat:
          plot_dir = misc.ensure_dir ('{0}/{1}/{2}/plots'.format (self.mode_dir,self.opts.cat,self.opts.weights))
        else:
          plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        #savingfig (fig, plot_dir,
        #        'int_sens_src_{0:04.1f}_test_{1:04.1f}'.format (
        #            self.opts.src_ext,
        #            self.opts.test_ext))
        #savingfig (fig, plot_dir, 'MESC_vs_track_sensitivity')
        if ext !=0.0:
          savingfig (fig, plot_dir, '{}_discovery_ext'.format(self.mode.words[0],ext))
        else:
          savingfig (fig, plot_dir, '{}_discovery'.format(self.mode.words[0]))

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
    def submit_compareTS (self):
        """Submit compare_TS jobs to cluster."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/compareTS/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))


        #Have to compare a range of decs, ras.
        dec_degs = np.arange (-80, 80, 10)
        ra_degs = np.arange (0, 360, 20)

        #for test_ext in np.arange(11):
        for ra_deg in ra_degs:
            for dec_deg in dec_degs:
                  command = '{} {} compareTS --conf-dirs={} ' \
                            ' --test-dec={:+07.2f}' \
                            ' --test-ra={:04.1f}' \
                            ' --seed=0'.format (
                                env_shell,
                                this_script,
                                confs,
                                dec_deg,
                                ra_deg,
                                )
                  label = 'do_compareTS__dec_{:+07.2f}__ra__{:04.1f}'.format (
                            dec_deg, ra_deg)
                  commands.append (command)
                  labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def compareTS (self):
        """
        Compare TS results between csky and skylab
        """

        seed = self.opts.seed
        n_trials = int (self.opts.n_trials)
        n_year = self.conf.root.n_year
        #Different parameters for a catalog vs. a test at 0 deg
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))

            if 'extension' in params.keys():
                test_ext = params['extension']
                src_ext = params['extension']
            else:
                test_ext = np.zeros_like(params['dec'])
                src_ext = np.zeros_like(params['dec'])
            if self.opts.weights:
                if self.opts.weights == 'equal':
                    print ("No weighting given - use equal weights")
                    weights = np.ones_like(params['dec'])
                else:
                    weights = params[str(self.opts.weights)]
            else:
                raise ValueError('Must specify a weighting scheme for stacking')
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            test_ext_deg = [ext for ext in test_ext]
            src_ext_deg = [ext for ext in src_ext]
        else:
            dec_deg = self.opts.test_dec
            ra_deg = self.opts.test_ra
            weights = np.atleast_1d(1.) #put in to eliminate redundancy for ps injection 
            test_ext_deg = 0.
            src_ext_deg = 0.

        #Here I'll include an optional mask to cut the number of sources down to the 10 largest ones
        #par_min_arg = np.argsort (params['flux2m'])[-10]
        #par_mask = params['flux2m'] <  params['flux2m'][par_min_arg]

        dec = np.atleast_1d(np.radians(dec_deg))
        ra = np.atleast_1d(np.radians(ra_deg))
        weights = np.atleast_1d(weights)
        test_ext = np.radians(test_ext_deg)
        src_ext = np.radians(src_ext_deg)

        seed = self.opts.seed
        np.random.seed (seed)
        analysis = self.analysis
        samples = analysis.keys

        try:
          gammas = analysis.pdf_ratios_energy['IC86'].gammas
        except:
          gammas = analysis.pdf_ratios_energy['IC40'].gammas

        #analysis.set_energy_pdfs (disable=True) #Include to test space-only distributions.
        data = hyp.DataHypothesis (analysis)
        ps = hyp.PointSourceHypothesis (analysis, dec, ra, 2, weights = weights, sigsub=False)
        tr = trial.PSTrialRunner (data, ps, ps.tests(gammas), cut_n_sigma=np.inf)
        L = tr.get_one_calc (TRUTH=True).Ls[0]


        #now let's import the llhmodels.
        llhs = self.llhs

        #Initialize the arrays to hold the pdfs.
        skylab_spacearray = []
        skylab_energyarray = []
        skylab_array = []


        best_gamma =4. 
        print(samples)

        for sample in samples:
          #now we replace the pdfs with that of skylab... then calculate again.

          ps_llh = llhs[sample]
          exp = ps_llh.exp 
          mc = ps_llh.mc 
          #ps_llh = PointSourceLLH(exp, mc, livetime=livetimes['{}'.format(sample)], llh_model=llhmodels[sample], mode = 'all', seed=0, scramble=False)
          #llhs.append(ps_llh)

          #Now that we have a llh object, I want to find which gamma is closest to the best_fit gamma.
          #loop through sources for effA * theo_weight
          m = []
          for i,src in zip(range(len(dec)),dec):
            m.append(ps_llh.llh_model.effA(src,**{'gamma':best_gamma})[0]*weights[i])

          norm = np.sum(m)
	  m    = np.asarray(m, dtype=float) / norm

          SoB_skylab_space = 0 # this will be added to per src.
          for i,src_dec,src_ra in zip(range(len(dec)),dec,ra):
            S_space = ps_llh.llh_model.signal(src_ra,src_dec,exp)
            B_space = ps_llh.llh_model.background(exp)
            #For now, I'll convolve the space pdf with the effA*theo_w values to make it similar to csky's space pdf.
            SoB_skylab_space += m[i]*(S_space / B_space)

          #whereas this spline is for energy.
          #spline = ps_llh.llh_model._w_spline_dict[(('gamma', 2.0),)]
          spline = ps_llh.llh_model._w_spline_dict[(('gamma', best_gamma),)]
          SoB_skylab_energy = np.exp(ps_llh.llh_model._spline_eval(spline,exp))

          SoB_skylab = SoB_skylab_space*SoB_skylab_energy 
          
          #append these to the arrays:
          skylab_spacearray.append(SoB_skylab_space)
          skylab_energyarray.append(SoB_skylab_energy)
          skylab_array.append(SoB_skylab)

        #I'll also create a llh for a multiyear sample.
        ps_llh=MultiPointSourceLLH()
        for sample,i in zip(samples,xrange(len(llhs))):
           ps_llh.add_sample(str(i), llhs[sample])
           print ("Adding sample" + str(llhs[sample]))

        #csky_result = L.ts_ns()
        csky_result = tr.get_one_ts_ns_gamma(TRUTH=True)
        print ('csky results (csky pdfs): '+ str(csky_result))

        #and find the ts_ns for that llh.
        skylab_result = ps_llh.fit_source(ra,dec,src_w=weights,scramble=False)#, approx_grad=True) - doesn't work currently?
        print ('skylab results: '+ str((skylab_result[0], skylab_result[1]['nsources'],skylab_result[1]['gamma'])))

        #now let's put these arrays into dictionaries so that csky recognizes them in the backend.
        skylab_space_dict = dict(zip(samples,skylab_spacearray))
        skylab_energy_dict = dict(zip(samples,skylab_energyarray))
        skylab_dict = dict(zip(samples,skylab_array))

        #Now replace each pdf with the skylab ones.
        calc = tr.get_one_calc (TRUTH=True)
        L = calc.Ls[np.argmin (np.abs (self.analysis.gammas - best_gamma))]

        L.backend = llh.csky_ext.PSLogLCalcBackend (L.pdf_ratio_space_keys, skylab_energy_dict,
                                        L.pdf_ratio_space_sigsub_keys, L.n_excluded_keys, L.ns_ratio_keys)

        print ('csky results (skylab pdfs): '+ str(L.ts_ns()))
        #print ('csky results (skylab pdfs): '+ str(tr.get_one_ts_ns_gamma(TRUTH=True)))

        if self.opts.cat:
          out_dir  = misc.ensure_dir ('{0}/pdfs/skylab/{1}'.format (self.mode_dir,self.opts.cat))
        else:
          out_dir  = misc.ensure_dir ('{0}/pdfs/skylab/dec_{1}/ra_{2}'.format (self.mode_dir,dec_deg,ra_deg))
        outfile = open(out_dir+"/compareTS.txt","w")
        outfile.write('csky results (csky pdfs): TS = {0:.2f}, fitted ns = {1:.2f}, fitted gamma = {2:.2f}'.format(csky_result[0],csky_result[1],csky_result[2]))
        outfile.write('\nskylab results : TS = {0:.2f}, fitted ns = {1:.2f}, fitted gamma = {2:.2f}'.format(skylab_result[0],skylab_result[1]['nsources'],skylab_result[1]['gamma']))
        outfile.write('\ncsky results (skylab pdfs) : TS = {0:.2f}, fitted ns = {1:.2f}'.format(L.ts_ns()[0],L.ts_ns()[1]))
        outfile.close()


        #Finally, let's save the skylab space and energy ratio stuff for mike to play with.
        saving (skylab_space_dict, out_dir+'/space_SoB.array') 
        saving (skylab_energy_dict, out_dir+'/energy_SoB.pickle') 
        saving (skylab_result, out_dir+'/results_skylab.pickle') 
        saving (csky_result, out_dir+'/results_csky.pickle')

    @command
    def submit_unblinding_test (self):
        """Submit unblinding job to cluster, for a growing catalog."""
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/unblinding_test/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir, memory=4)
        commands, labels = [], []
        #get catalog size for the purposes of doing batches of individual unblinding tests.
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
        N = len(params['dec'])
        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))
        if self.opts.dec_order:
          command1 = '{} {} unblind_sources_test --conf-dirs={} ' \
                            ' --cat={}' \
                            ' --weights={}' \
                            ' --dec-order'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.cat,
                                self.opts.weights)
          label1 = 'do_unblind_test_{}_{}_sn'.format (
                  self.opts.cat, self.opts.weights)

          command2 = '{} {} unblind_sources_test --conf-dirs={} ' \
                            ' --cat={}' \
                            ' --weights={}' \
                            ' --dec-order' \
                            ' --ns'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.cat,
                                self.opts.weights)
          label2 = 'do_unblind_test_{}_{}_ns'.format (
                  self.opts.cat, self.opts.weights)
          commands.append (command1)
          commands.append (command2)
          labels.append (label1)
          labels.append (label2)

        elif self.opts.weight_order:
          command = '{} {} unblind_sources_test --conf-dirs={} ' \
                            ' --cat={}' \
                            ' --weights={}' \
                            ' --weight-order'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.cat,
                                self.opts.weights)
          label = 'do_unblind_test_{}_{}_weights'.format (
                  self.opts.cat, self.opts.weights)
          commands.append (command)
          labels.append (label)
        #for individual sources in the catalog:
        else:
          for i in range(N):
            command = '{} {} unblind_sources_test --conf-dirs={} ' \
                            ' --cat={}' \
                            ' --weights={}' \
                            ' --i-source={}'.format (
                                env_shell,
                                this_script,
                                confs,
                                self.opts.cat,
                                self.opts.weights,
                                i)
            label = 'do_unblinding_test_{}'.format (i)
            commands.append (command)
            labels.append (label)

        submitter.submit_condor00 (commands, labels,
                blacklist=self.blacklist,
                )
        pjobdir (job_dir)

    @command
    def unblind_sources_test (self):
        """
        Using for determining how many sources change the unblinded TS_ns of a catalog
        """
        job_root = self.conf.root.cluster_job_dir
        job_dir = '{0}/do_bg_trials/{1}'.format (job_root, job_id)
        submitter = Submitter (job_dir=job_dir,max_jobs = 100, memory=4)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (
            os.path.abspath, self.conf.root_dirs))

        seed = self.opts.seed
        n_year = self.conf.root.n_year
        #Different parameters for a catalog vs. a test at 0 deg
        params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))

        if 'extension' in params.keys():
            test_ext = params['extension']
            src_ext = params['extension']
        else:
            test_ext = np.zeros_like(params['dec'])
            src_ext = np.zeros_like(params['dec'])
        if self.opts.weights:
            if self.opts.weights == 'equal':
                print ("No weighting given - use equal weights")
                weights = np.ones_like(params['dec'])
            else:
                weights = params[str(self.opts.weights)]
        else:
            raise ValueError('Must specify a weighting scheme for stacking')
        dec_deg = np.degrees(params['dec']) 
        ra_deg = np.degrees(params['ra']) 
        test_ext_deg = [ext for ext in test_ext]
        src_ext_deg = [ext for ext in src_ext]

        dec = np.atleast_1d(np.radians(dec_deg))
        ra = np.atleast_1d(np.radians(ra_deg))
        weights = np.atleast_1d(weights)
        test_ext = np.radians(test_ext_deg)
        src_ext = np.radians(src_ext_deg)

        seed = self.opts.seed
        np.random.seed (seed)
        analysis = self.analysis
        #To include sources in dec order, not weighted order
        dec_order = self.opts.dec_order
        i_source = self.opts.i_source
        try:
          gammas = analysis.pdf_ratios_energy['IC86'].gammas
        except:
          gammas = analysis.pdf_ratios_energy['IC40'].gammas

        #analysis.set_energy_pdfs (disable=True) #Include to test space-only distributions.
        data = hyp.DataHypothesis (analysis)
        N_sources = []
        ts = []
        ns = []
        gamma_fit = []
        newdec = []
        best_gamma = 2.

        if self.opts.i_source: #if I'm just doing one source at a time:
          par_mask = i_source #this is my clever way of choosing the right source.
          newdec.append(dec[par_mask])
          ps = hyp.PointSourceHypothesis (analysis, dec[par_mask], ra[par_mask], 2, weights = weights[par_mask], sigsub=False, can_inject=False)
          tr = trial.PSTrialRunner (data, ps, ps.tests(gammas), cut_n_sigma=np.inf)
          L = tr.get_one_calc (TRUTH=True).Ls[0]
          result = tr.get_one_ts_ns_gamma(TRUTH=True)
          N_sources.append(1) #only using one at a time
          ts.append(result[0])
          ns.append(result[1])
          gamma_fit.append(result[2])
          csky_result = {'ts':result[0],'ns':result[1],'gamma':result[2], 'newdec':newdec[0]}
          prush ('{0} unblinded, {1}th source: '.format(self.opts.cat,i_source) + str(csky_result))

        else:
          for N in range(1,len(params['dec'])):
            if self.opts.weights == 'equal' or dec_order == True:
              if self.opts.ns:
                par_mask = np.argsort (dec)[-N:]
                newdec.append(dec[par_mask][0])
              else:
                par_mask = np.argsort (dec)[:N]
                newdec.append(dec[par_mask][-1])
            else:
              par_mask = np.argsort (weights)[-N:]
              newdec.append(dec[par_mask][0])
            ps = hyp.PointSourceHypothesis (analysis, dec[par_mask], ra[par_mask], 2, weights = weights[par_mask], sigsub=False, can_inject=False)
            tr = trial.PSTrialRunner (data, ps, ps.tests(gammas), cut_n_sigma=np.inf)
            L = tr.get_one_calc (TRUTH=True).Ls[0]
            result = tr.get_one_ts_ns_gamma(TRUTH=True)
            N_sources.append(N)
            ts.append(result[0])
            ns.append(result[1])
            gamma_fit.append(result[2])
            csky_result = {'ts':result[0],'ns':result[1],'gamma':result[2], 'newdec':newdec[0]}
            prush ('{0} unblinded, {1} top sources: '.format(self.opts.cat,N) + str(csky_result))

        #now let's make it a dict.
        results = {'N':tuple(N_sources),'ts':tuple(ts),'ns':tuple(ns),'gamma':tuple(gamma_fit), 'newdec':tuple(newdec)} 

        print ('Complete!')
        out_dir  = misc.ensure_dir ('{0}/{1}/unblinding_test'.format (self.mode_dir,self.opts.cat))
        if dec_order:
          if self.opts.ns:
            saving (results, out_dir+'/results_decorder_ns.pickle')
          else:
            saving (results, out_dir+'/results_decorder_sn.pickle')
        elif self.opts.weight_order:
          saving (results, out_dir+'/results_weightorder.pickle')
        else:
          saving (results, out_dir+'/batch_{}_results.pickle'.format(self.opts.i_source,self.opts.weights))

    @command
    def plot_unblinding_test (self): #needs updating
        """
        Plot unblinding test with increasing numbers of sources. 
        """

        result_dir  = misc.ensure_dir ('{0}/{1}/unblinding_test'.format (self.mode_dir,self.opts.cat))
        dec_order = self.opts.dec_order

        #gather data:
        if self.opts.weight_order: 
          results = cache.load(result_dir+'/results_weightorder.pickle')
        elif dec_order:
          if self.opts.ns:
              results = cache.load(result_dir+'/results_decorder_ns.pickle')
          else:
              results = cache.load(result_dir+'/results_decorder_sn.pickle')
        else: #for indiv sources
          files = [cache.load(result_dir + '/' + file) for file in os.listdir(result_dir + '/') if file.startswith('batch')]
          results = {}
          for key in files[0].iterkeys():
            results[key] = tuple(file[key][0] for file in files)

        N,ts,ns,gamma,newdec = results['N'],results['ts'],results['ns'],results['gamma'],results['newdec']

        if dec_order:
          if self.opts.ns:
            plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}_dec_order_ns'.format (self.mode_dir,self.opts.cat,self.opts.weights))
          else:
            plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}_dec_order_sn'.format (self.mode_dir,self.opts.cat,self.opts.weights))
          xlabel = r'$\sin(\delta)$'
          x_axis = newdec
          xmin = -1.0
          xmax = 1.0
        elif self.opts.weight_order:
          plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}_weight_order'.format (self.mode_dir,self.opts.cat,self.opts.weights))
          xlabel = 'N sources (highest weight first)'
          x_axis = N
          xmin = -10
          xmax = len(N)+10
        else: #indiv sources
          plot_dir  = misc.ensure_dir ('{0}/cats/{1}/{2}'.format (self.mode_dir,self.opts.cat,self.opts.weights))
          xlabel = r'$\sin(\delta)$'
          x_axis = newdec
          xmin = -1.0
          xmax = 1.0

        ts_label = 'TS'
        ns_label = 'signal-like events (ns)'
        gamma_label = r'fitted $\gamma$'

        fig, ax = pfig()
        ax.cla()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ts_label)
        ax.set_ylim(0,12.)
        ax.grid()
        plt.xlim(xmin,xmax)
        ax.scatter (x_axis, ts, color = 'k', s=0.7,alpha=0.7)
        ax.set_title('Unblinding Test for Growing Catalog - {}'.format(self.opts.cat))
        savingfig (fig, plot_dir, 'TS') 
        plt.close(fig)

        fig, ax = pfig()
        ax.cla()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ns_label)
        plt.xlim(xmin,xmax)
        plt.ylim(0,60)
        ax.scatter (x_axis, ns, color = 'k', s=0.7,alpha=0.7)
        ax.set_title('Unblinding Test for Growing Catalog - {}'.format(self.opts.cat))
        savingfig (fig, plot_dir, 'ns') 
        plt.close(fig)

        fig, ax = pfig()
        ax.cla()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(gamma_label)
        plt.xlim(xmin,xmax)
        ax.scatter (x_axis, gamma, color = 'k', s=0.7,alpha=0.7)
        ax.set_title('Unblinding Test for Growing Catalog - {}'.format(self.opts.cat))
        ax.set_ylim(0.,5.)
        savingfig (fig, plot_dir, 'gamma') 
        plt.close(fig)

    @command
    def plot_ts_ns_ratio (self):
        """
        Plot ratio between ts_ns for csky and skylab across the sky
        """

        if self.opts.cat:
          result_dir  = misc.ensure_dir ('{0}/pdfs/skylab/{1}'.format (self.mode_dir.opts.cat))
        else:
          result_dir  = misc.ensure_dir ('{0}/pdfs/skylab'.format (self.mode_dir))

        #gather data
        csky_files = [cache.load(file) for file in glob(result_dir+'/*/*/*csky.pickle')]
        skylab_files = [cache.load(file) for file in glob(result_dir+'/*/*/*skylab.pickle')]

        TS_csky=[]
        TS_skylab=[]
        n_csky=[]
        n_skylab=[]

        for file in csky_files:
          TS_csky.append(file[0])
          n_csky.append(file[1])

        for file in skylab_files:
          TS_skylab.append(file[0])
          n_skylab.append(file[1]['nsources'])

        if self.opts.cat:
          plot_dir  = misc.ensure_dir ('{0}/plots/ts_ns_skylab/{1}'.format (self.mode_dir,self.opts.cat))
        else:
          plot_dir  = misc.ensure_dir ('{0}/plots/ts_ns'.format (self.mode_dir))

          xlabel = 'csky'
          ylabel = 'skylab'

          fig, ax = pfig()
          ax.cla()

          ax.set_xlabel(xlabel)
          ax.set_ylabel(ylabel)
          ax.set_xlim(0.,10.)
          ax.set_ylim(0.,10.)

          ax.scatter (TS_csky, TS_skylab, color = 'k', s=0.5,alpha=0.5)

          ax.set_title('TS Ratio (csky / skylab)')
          #ax.loglog()
          ax.set_aspect('equal')
          ax.plot ([0,10], [0,10], color='.5', ls='--')
          savingfig (fig, plot_dir, 'TS_ratio') 
          plt.close(fig)

          fig, ax = pfig()
          ax.cla()

          ax.set_xlabel(xlabel)
          ax.set_ylabel(ylabel)
          ax.set_xlim(0.,10.)
          ax.set_ylim(0.,10.)

          ax.scatter (n_csky, n_skylab, color = 'k', s=0.5,alpha=0.5)

          ax.set_title('ns Ratio (csky / skylab)')
          #ax.loglog()
          ax.set_aspect('equal')
          ax.plot ([0,10], [0,10], color='.5', ls='--')
          savingfig (fig, plot_dir, 'nsources_ratio') 
          plt.close(fig)

          #Alright, now I'd like to plot some hists showing the distribution of TS and ns values for each of these. Maybe 2d hist for each one?
          #h_csky = histlite.hist((TS_csky,n_csky), bins = (21,21), range = ((0,10),(0,10)))
          #h_skylab = histlite.hist((TS_skylab,n_skylab), bins = (21,21), range = ((0,10),(0,10)))
          csky_label = 'csky'
          skylab_label = 'skylab'
          bins = 41
          plotrange = (0,10)
          h_csky_TS = histlite.hist(TS_csky, bins = bins, range = plotrange)
          h_skylab_TS = histlite.hist(TS_skylab, bins = bins, range = plotrange)
          h_csky_ns = histlite.hist(n_csky, bins = bins, range = plotrange)
          h_skylab_ns = histlite.hist(n_skylab, bins = bins, range = plotrange)

          fig, ax = pfig()
          ax.cla()
          ax.set_title('TS distribution over whole sky')

          ax.set_xlabel('TS')
          ax.semilogy()
          ax.set_xlim(0,10.)
          ax.set_ylim(1e-1,1e3)
          histlite.plot1d (ax, h_csky_TS, label = csky_label)
          histlite.plot1d (ax, h_skylab_TS, label = skylab_label)
          ax.legend(loc='upper right', prop=propsmall)
          savingfig (fig, plot_dir, 'TS_dist') 
          plt.close(fig)

          fig, ax = pfig()
          ax.cla()
          ax.set_title('ns distribution over whole sky')

          ax.set_xlabel('ns')
          ax.semilogy()
          ax.set_xlim(0,10.)
          ax.set_ylim(1e-1,1e3)
          histlite.plot1d (ax, h_csky_ns, label = csky_label)
          histlite.plot1d (ax, h_skylab_ns, label = skylab_label)
          ax.legend(loc='upper right', prop=propsmall)
          savingfig (fig, plot_dir, 'ns_dist') 
          plt.close(fig)

    @command
    def plot_pdf_ratio (self):
        """
        Plot ratio between SoB for csky and skylab
        """

        #Different parameters for a catalog vs. a test at 0 deg
        if self.opts.cat:
            params = cache.load('/data/condor_builds/users/brelethford/Data/{}/pickle/params.pickle'.format(self.opts.cat))
            if 'extension' in params.keys():
                test_ext = params['extension']
                src_ext = params['extension']
            else:
                test_ext = np.zeros_like(params['dec'])
                src_ext = np.zeros_like(params['dec'])
            if self.opts.weights:
                if self.opts.weights == 'equal':
                    print ("No weighting given - use equal weights")
                    weights = np.ones_like(params['dec'])
                else:
                    weights = params[str(self.opts.weights)]
            else:
                raise ValueError('Must specify a weighting scheme for stacking')
            dec_deg = np.degrees(params['dec']) 
            ra_deg = np.degrees(params['ra']) 
            test_ext_deg = [ext for ext in test_ext]
            src_ext_deg = [ext for ext in src_ext]
        else:
            dec_deg = self.opts.test_dec
            ra_deg = self.opts.test_ra
            weights = np.atleast_1d(1.) #put in to eliminate redundancy for ps injection 
            test_ext_deg = 0.
            src_ext_deg = 0.

        dec = np.atleast_1d(np.radians(dec_deg))
        ra = np.atleast_1d(np.radians(ra_deg))
        weights = np.atleast_1d(weights)
        test_ext = np.radians(test_ext_deg)
        src_ext = np.radians(src_ext_deg)

        #csky setup
        analysis = self.analysis
        data = hyp.DataHypothesis (analysis)
        try:
          gammas = analysis.pdf_ratios_energy['IC86'].gammas
        except:
          gammas = analysis.pdf_ratios_energy['IC40'].gammas

        #skylab setup
        llhs = self.llhs
        gammas = np.arange(1.,4.5,1.0)
        #gammas = np.array([4.])

        #In order to correctly find the right fit_gamma, we change the analysis gammas:
        analysis.gammas = gammas
        samples = analysis.keys

        for gamma in gammas:
          ps = hyp.PointSourceHypothesis (analysis, dec, ra, gamma, weights = weights, sigsub=False)
          tr = trial.PSTrialRunner (data, ps, ps.tests (gammas), cut_n_sigma=np.inf)
          calc = tr.get_one_calc (TRUTH=True)
          L_2_00 = calc.Ls[np.argmin (np.abs (gamma - analysis.gammas))]

          for sample in samples:
            #in csky, space SoB has theo_weights and effA included. No need to explicitly calculate.
            SoB_csky_space = L_2_00.pdf_ratio_space_keys[sample]
            SoB_csky_energy = L_2_00.pdf_ratio_energy_keys[sample]
            SoB_csky = SoB_csky_space * SoB_csky_energy

            #now skylab: Use the premade llh from setup_ana
            llh = llhs[sample]
	    exp = llh.exp
            mc = llh.mc
            livetime = llh.livetime

            #loop through sources for effA * theo_weight
            m = []
            for i,src in zip(range(len(dec)),dec):
              m.append(llh.llh_model.effA(src,**{'gamma':gamma})[0]*weights[i])

            norm = np.sum(m)
	    m    = np.asarray(m, dtype=float) / norm

            SoB_skylab_space = 0 # this will be added to per src.
            for i,src_dec,src_ra in zip(range(len(dec)),dec,ra):
              S_space = llh.llh_model.signal(src_ra,src_dec,exp)
              B_space = llh.llh_model.background(exp)
              #For now, I'll convolve the space pdf with the effA*theo_w values to make it similar to csky's space pdf.
              SoB_skylab_space += m[i]*(S_space / B_space)

            #Simplify the plots by removing the very low values for SoB
            spacemask = SoB_csky_space > 1e-5

            #whereas this spline is for energy.
            try:
              spline = llh.llh_model._w_spline_dict[(('gamma', gamma),)]
            except:
              spline = llh.llh_model._w_spline_dict[(('gamma', 2.),)]
            SoB_skylab_energy = np.exp(llh.llh_model._spline_eval(spline,exp))

            SoB_skylab = SoB_skylab_space*SoB_skylab_energy 
            
            if self.opts.cat:
              plot_dir  = misc.ensure_dir ('{0}/plots/SoB_vs_skylab/{1}/{2}/{3}/{4}'.format (self.mode_dir,sample,self.opts.cat,self.opts.weights,gamma))
            else:
              plot_dir  = misc.ensure_dir ('{0}/plots/SoB_vs_skylab/{1}/dec_{2}/{3}'.format (self.mode_dir,sample,dec_deg,gamma))

            #Let's do a plot for effA*theo_w if we're dealing with a catalog. 

            if self.opts.cat:
              fig, ax = pfig()
              fig.clf()
              gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
              ax = plt.subplot (gs[0])
              rax = plt.subplot (gs[1], sharex=ax)
              ax.scatter (np.sin(dec), m, color = 'blue', marker='d', label = 'skylab')
              ax.scatter (np.sin(dec), np.concatenate(tr.get_test_gamma(gamma).p_keys_sources[sample],axis=0), color = 'red', label = 'csky')
              #ax.scatter (np.sin(dec), ps.p_keys_sources[sample], color = 'red', label = 'csky')
              ax.set_title('acceptance weights')
              ax.set_ylabel('effA * theo_w')
              ax.set_xlim(-1.,1.)
              ax.set_ylim(bottom=0.)
              ax.legend(loc='upper right', prop=propsmall)
              rax.scatter (np.sin(dec), np.concatenate(tr.get_test_gamma(gamma).p_keys_sources[sample],axis=0)/m)
        
              rax.set_xlabel (r'$\sin(\delta)$')
              rax.set_ylabel ('csky / skylab')
              rax.set_xlim (-1, 1)

              rax.grid ()
              savingfig (fig, plot_dir, 'acceptance') 
              plt.close(fig) 

            #Plot for just detector acceptance, no weights (effA)

            test_decs = np.radians(np.arange(-90.,91.,1))
            effA_skylab = [llh.llh_model.effA(src,**{'gamma':gamma})[0]  for src in test_decs]
            #hist it so we can view both at once
            h_skylab = histlite.hist(effA_skylab, bins=100, range = (-1.,1.))  

            effA_hist = analysis.dists_dec_sig[sample].dpdfs[gamma].h.normalize(integrate=True)
            maxval = max(effA_hist.values)

            #use maxval to get effA skylab to be at the same relative level.
            fig, ax = pfig()
            gs = mpl.gridspec.GridSpec (2, 1, height_ratios=[3,1], hspace=.15)
            ax = plt.subplot (gs[0])
            rax = plt.subplot (gs[1], sharex=ax)

            ax.plot (np.sin(test_decs), effA_skylab / (max(effA_skylab) / maxval), color = 'blue', label = 'skylab')
             
            #histlite.plot1d (ax, h_skylab.normalize(integrate = True), color = 'blue', label = 'skylab')
            #histlite.plot1d (ax, histlite.hist_from_function(effA_hist.bins,effA_hist.spline_fit().spline), color = 'red', label = 'csky')
            ax.plot (np.sin(test_decs), effA_hist.spline_fit().spline(np.sin(test_decs)), color = 'red', label = 'csky')
            ax.set_title('Detector Acceptance - $\gamma = ${}'.format(gamma))
            ax.set_ylabel('effA - normalized')
            ax.set_xlim(-1.,1.)
            ax.set_ylim(bottom=0.)
            ax.legend(loc='upper right', prop=propsmall)
            rax.scatter (np.sin(test_decs), effA_hist.spline_fit().spline(np.sin(test_decs))/(effA_skylab / (max(effA_skylab) / maxval)))
        
            rax.set_xlabel (r'$\sin(\delta)$')
            rax.set_ylabel ('csky / skylab')
            rax.set_xlim (-1, 1)

            rax.grid ()
            savingfig (fig, plot_dir, 'detector_only') 
            plt.close(fig)
            #Plots for space and energy hists for each of csky and skylab.

            #SoB - Space
            h_skylab =  histlite.hist(SoB_skylab_space, bins=101, range=(1e-8,1e4), log=True) 
            h_csky =  histlite.hist(SoB_csky_space, bins=101, range=(1e-8,1e4), log=True) 

            fig_space, ax_space = pfig()
            ax_space.cla()
            histlite.plot1d (ax_space, h_skylab, log=True, color = 'blue', alpha = 0.5, label = 'skylab')
            histlite.plot1d (ax_space, h_csky, log=True, color = 'red', alpha = 0.5, label = 'csky')
            ax_space.set_title('SoB - Space')
            ax_space.semilogx()
            ax_space.legend(loc='upper right', prop=propsmall)
            savingfig (fig_space, plot_dir, 'SoB_space_both')   
            plt.close(fig_space)

            #SoB - Energy
            h_skylab =  histlite.hist(SoB_skylab_energy, bins=101, range=(1e-4,1e4), log=True) 
            h_csky =  histlite.hist(SoB_csky_energy, bins=101, range=(1e-4,1e4), log=True) 

            fig_energy, ax_energy = pfig()
            ax_energy.cla()
            histlite.plot1d (ax_energy, h_skylab, log=True, color = 'blue', alpha = 0.5, label = 'skylab')
            histlite.plot1d (ax_energy, h_csky, log=True, color = 'red', alpha = 0.5, label = 'csky')
            ax_energy.set_title('SoB - Energy')
            ax_energy.semilogx()
            ax_energy.legend(loc='upper right', prop=propsmall)
            savingfig (fig_energy, plot_dir, 'SoB_energy_both')   
            plt.close(fig_energy)

            #Plots for the hist ratios of csky vs skylab
            h_space =  histlite.hist(SoB_csky_space[spacemask]/ SoB_skylab_space[spacemask], bins=101, range=(1e-1,1e1), log=True) 
            h_energy =  histlite.hist(SoB_csky_energy / SoB_skylab_energy, bins=101, range=(1e-1,1e1), log=True) 
            h_ratio =  histlite.hist(SoB_csky / SoB_skylab, bins=101, range=(1e-1,1e1), log=True) 


            #now the llh plots
            fig_space, ax_space = pfig()
            fig_energy, ax_energy = pfig()
            fig, ax = pfig()
            ax_space.cla()
            ax_energy.cla()
            ax.cla()

            histlite.plot1d (ax_space, h_space, range=(1e-2,1e2), log=True)
            histlite.plot1d (ax_energy, h_energy, range=(1e-2,1e2), log=True)
            histlite.plot1d (ax, h_ratio, range=(1e-2,1e2), log=True)
            ax_space.set_title('SoB Ratio - Space Only (csky / skylab)')
            ax_space.semilogx()
            plt.axvline (1, color='.5', ls='--')
            savingfig (fig_space, plot_dir, 'SoB_space_ratio_log')  
            plt.close(fig_space)

            ax_energy.set_title('SoB Ratio - Energy Only (csky / skylab)')
            plt.axvline (1, color='.5', ls='--')
            ax_energy.semilogx()
            savingfig (fig_energy, plot_dir, 'SoB_energy_ratio_log')  
            plt.close(fig_energy)

            ax.set_title('SoB Ratio (csky / skylab)')
            ax.semilogx()
            plt.axvline (1, color='.5', ls='--')
            savingfig (fig, plot_dir, 'SoB_ratio_log') 
            plt.close(fig)

            #now some scatter plots.
            fig_space, ax_space = pfig()
            fig_energy, ax_energy = pfig()
            fig, ax = pfig()
            ax_space.cla()
            ax_energy.cla()
            ax.cla()

            xlabel = 'csky'
            ylabel = 'skylab'

            ax_space.set_xlabel(xlabel)
            ax_space.set_ylabel(ylabel)
            ax_energy.set_xlabel(xlabel)
            ax_energy.set_ylabel(ylabel)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax_space.scatter (SoB_csky_space, SoB_skylab_space, color = 'blue', s=0.5,alpha=0.1)
            ax_energy.scatter (SoB_csky_energy, SoB_skylab_energy, color = 'red', s=0.5,alpha=0.1)
            ax.scatter (SoB_csky, SoB_skylab, color = 'k', s=0.5,alpha=0.1)

            ax_space.set_title('SoB Ratio - Space Only')
            ax_space.loglog()
	    min_space = min(min(SoB_csky_space),min(SoB_skylab_space))
            max_space = max(max(SoB_csky_space),max(SoB_skylab_space))

            ax_space.set_xlim(1e-1,max_space)
            ax_space.set_ylim(1e-1,max_space)
            ax_space.set_aspect('equal')
            ax_space.plot ([1e-1,max_space], [1e-1,max_space], color='.5', ls='--')
            savingfig (fig_space, plot_dir, 'SoB_space_scatter')  
            plt.close(fig_space)

            ax_energy.set_title('SoB Ratio- Energy Only')
            ax_energy.loglog()
	    min_energy = min(min(SoB_csky_energy),min(SoB_skylab_energy))
            max_energy = max(max(SoB_csky_energy),max(SoB_skylab_energy))
            ax_energy.set_xlim(min_energy,max_energy)
            ax_energy.set_ylim(min_energy,max_energy)
            ax_energy.set_aspect('equal')
            ax_energy.plot ([min_energy,max_energy], [min_energy,max_energy], color='.5', ls='--')
            savingfig (fig_energy, plot_dir, 'SoB_energy_ratio_scatter')  
            plt.close(fig_energy)

            ax.set_title('SoB Ratio (csky / skylab)')
            ax.loglog()
	    min_tot = min(min(SoB_csky),min(SoB_skylab))
            max_tot = max(max(SoB_csky),max(SoB_skylab))
            ax.set_xlim(1e-4,1e4)
            ax.set_ylim(1e-4,1e4)
            #ax.set_xlim(min_tot,max_tot)
            #ax.set_ylim(min_tot,max_tot)
            ax.set_aspect('equal')
            ax.plot ([min_tot,max_tot], [min_tot,max_tot], color='.5', ls='--')
            savingfig (fig, plot_dir, 'SoB_ratio_scatter') 
            plt.close(fig)

            #plot the top 10 SoB space ratio sources
            N = 10
            SoB_space_ratio = SoB_csky_space[spacemask] /  SoB_skylab_space[spacemask]
            SoB_energy_ratio = SoB_csky_energy / SoB_skylab_energy 
            high_args = (-1*SoB_space_ratio).argsort()[:N]
            ev_dec, ev_ra = np.degrees(exp['dec'][high_args]), np.degrees(exp['ra'][high_args])
            cat = self.opts.cat
            #catalog info:        
            weightplot=[(i*100)/max(weights) for i in weights]
             
            plot_dir  = misc.ensure_dir (self.root_dir+'/catalogs/{}'.format(cat))
            fig = plt.plot()
            #Rot = 180 shifts the axis from the center to the edge
            #centered by default
            healpy.mollview(title = 'Equatorial map of {} - top ten event weights'.format(cat), cbar = False,
                        rot = 180, notext = True, cmap=None, coord='C')
            healpy.graticule(coord='C', color='DimGrey')
            py.title("{} - {} - Equatorial".format(cat,self.opts.weights), fontsize = 25, fontweight='bold')
            #healpy.projscatter(0,0,coord='G',lonlat=True) # This one used to test GC
            healpy.projscatter(ra_deg,dec_deg,coord='C',lonlat=True, s = weightplot)
            healpy.projscatter(ev_ra,ev_dec,coord='C',lonlat=True, color = 'red')
            healpy.projtext(185, 2,'180', lonlat=True, fontweight = 'bold')
            healpy.projtext(95, 2, '90', lonlat=True, fontweight = 'bold')
            healpy.projtext(275, 2, '270', lonlat=True, fontweight = 'bold')
            healpy.projtext(8, 2,'0', lonlat=True, fontweight = 'bold')
            healpy.projtext(359, 2, '360', lonlat=True, fontweight = 'bold')
            healpy.projtext(193, -8, 'RA (deg)', lonlat=True, fontweight = 'bold')
            healpy.projtext(350, 30, '30', lonlat=True, fontweight = 'bold')
            healpy.projtext(340, 60, '60', lonlat=True, fontweight = 'bold')
            #healpy.projtext(5, -5, 'Dec (deg)', lonlat=True)
            healpy.projtext(358, -33.5, '-30', lonlat=True, fontweight = 'bold')
            healpy.projtext(358, -63.5, '-60', lonlat=True, fontweight = 'bold')
            ra_gplane = np.arange(0.,361.,1.)
            dec_gplane = np.zeros(len(ra_gplane))
            
            healpy.projplot(ra_gplane, dec_gplane, coord='G', lonlat=True, color='DimGrey', linewidth=2., alpha=0.5)
            fname = '/{0}_{1}_equatorial_with_events'.format(cat,self.opts.weights)
            py.savefig(plot_dir+fname)

            public_dir = '/home/brelethford/public_html/csky/' + plot_dir.split('brelethford/')[-1]
            print (public_dir+fname)
            subprocess.call(['ssh', 'cobalt06', 'mkdir', '-p', '{0}/'.format(public_dir)])
            subprocess.call(["scp", '-r','{0}/'.format (plot_dir), 'brelethford@cobalt06:{}'.format(public_dir)])

            ev_space =  SoB_space_ratio[high_args]
            ev_energy = SoB_energy_ratio[high_args]
            print ("SoB space ratios:")
            print (ev_space)
            print ("SoB energy ratios:")
            print (ev_energy)
            print ("SoB csky space:")
            print(SoB_csky_space[spacemask][high_args])
            print ("SoB csky energy:")
            print(SoB_csky_energy[spacemask][high_args])
            print ("SoB skylab space:")
            print(SoB_skylab_space[spacemask][high_args])
            print ("SoB skylab energy:")
            print(SoB_skylab_energy[spacemask][high_args])


if __name__ == '__main__':
    self = Csky6yr ()
    self.run ()
    try:
        import pandas as pd
    except:
        pass


hl = histlite


