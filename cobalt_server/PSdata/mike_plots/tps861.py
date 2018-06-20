#!/usr/bin/env python

from __future__ import print_function

import copy
from glob import glob
from itertools import izip
import os
import re
import socket
import sys
import time

import matplotlib as mpl
if socket.gethostname () not in ('ordi', 'zgaskami'):
    mpl.use ('Agg')
import healpy
import matplotlib.pyplot as plt
import numpy as np
pi = np.pi
from optparse import OptionParser
import scipy.optimize, scipy.stats
import tables

from optparse import OptionParser

from icecube.umdtools import arrays, cache, misc, submitter
from icecube.umdtools.apps import Timed, Driver, command
from icecube.umdtools.conf import Configurator
from icecube.umdtools.vars_class import Vars
from icecube.csky import bk, cat, coord, ens, pdf, trial
from icecube import histlite as hl
ensure_dir = misc.ensure_dir


no_emojify = lambda *a: a[0]

try:
    from emojify import emojify
except:
    emojify = no_emojify

if socket.gethostname () not in ('ordi', 'zgaskami', 'condor00'):
    emojify = no_emojify

job_id = '{}_T{:.0f}_J{}'.format (
        socket.gethostname (), time.time (), os.getpid ())

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

def pjobdir (job_dir):
    prush ('\nJob dir is {0} .\n'.format (job_dir))


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


class TPS861 (Timed, Driver):

    def __init__ (self):
        Timed.__init__ (self)
        Driver.__init__ (self)

    def run (self, arglist=[]):
        usage = '%prog {[options] [commands]}\n' + self._command_help
        self.parser = parser = OptionParser (usage=usage)

        parser.add_option ('--sigma', dest='sigma',
                default=0, type=int, metavar='N',
                help='handle N-sigma calculations')

        parser.add_option ('--beta', dest='beta',
                default=None, type=float, metavar='BETA',
                help='must surpass threshold in BETA fraction of trials')

        parser.add_option ('--zenith', dest='zenith',
                default=90, type=float, metavar='ZENITH',
                help='test and injection point is at ZENITH (CC, deg)')

        parser.add_option ('--gamma', dest='gamma',
                default=2, type=float, metavar='GAMMA',
                help='source has spectral index GAMMA')

        parser.add_option ('--n-trials', dest='n_trials',
                default=100, type=float, metavar='N',
                help='perform N trials')

        parser.add_option ('--n-jobs', dest='n_jobs',
                default=1., type=float, metavar='N',
                help='perform N jobs (with --n-trials each)')

        parser.add_option ('--seed', dest='seed',
                default=0, type=int, metavar='SEED',
                help='initialize RNG with SEED')

        parser.add_option ('-c', '--conf-dirs', dest='conf_dirs',
                default=os.path.abspath ('conf'), metavar='DIRS',
                help='load configuration from comma-separated list of DIRS')

        parser.add_option ('--blacklist', dest='blacklist',
                default='cobol61,cobol65', help='nodes to avoid')

        self.opts, self.commands = opts, commands = \
                parser.parse_args (arglist if arglist else sys.argv[1:])

        self.conf = Configurator (*self.opts.conf_dirs.split (','))
        self.mode = Vars ()
        self.mode.str = self.conf.root.mode
        self.mode.words = self.mode.str.split ('/')

        self.go (commands, announcement=
                 lambda s :emojify (
                     ':penguin: :penguin: :penguin: '
                     '{0} :penguin: :penguin: :penguin:'.format (s),
                     False))

    @property
    def blacklist (self):
        return [s+'.private.pa.umd.edu'
                for s in self.opts.blacklist.split (',')]


    @property
    def root_dir (self):
        return self.conf.root.root_dir

    @property
    def mode_dir (self):
        return ensure_dir ('{0}/{1}'.format (self.root_dir, self.mode.str))

    @property
    def data (self):
        try:
            return self._data
        except:
            self._data = loading (
                '{}/data/data.arrays'.format (self.mode_dir))
            return self._data

    @property
    def nu (self):
        try:
            return self._nu
        except:
            self._nu = loading (
                '{}/data/nu.arrays'.format (self.mode_dir))
            return self._nu

    @property
    def pdfs (self):
        try:
            return self._pdfs
        except:
            self._pdfs = loading (
                '{}/data/pdfs.pdfs'.format (self.mode_dir))
            return self._pdfs

    def bg (self):
        data = self.data
        bg = ens.DataInjector (data.sigma, np.arccos(data.sinDec), data.logE,
                               capsize=np.radians (3))
        return bg

    def sig (self, CC_zenith, gamma=None):
        if gamma is None:
            gamma = self.opts.gamma
        nu = self.nu
        sig = ens.PointInjector (
            nu.sigma, CC_zenith, pi/2, nu.trueDec+pi/2,
            nu.xaxis_zenith, nu.xaxis_azimuth, nu.logE,
            nu.ow*nu.trueE**-gamma,
            primary_logenergys=np.log10 (nu.trueE),
            logenergy_min=2, logenergy_max=np.inf, extension=0
        )
        return sig


    @property
    def bg_tsds (self):
        if hasattr (self, '_bg_tsds'):
            return self._bg_tsds
        try:
            self._bg_tsds = loading ('{0}/bg_tsds.dict'.format (self.mode_dir))
        except:
            self._bg_tsds = self.collect_bg_ts ()
        return self._bg_tsds

    @command
    def collect_bg_ts (self):
        """Collect bg_trials dict and cache it."""
        prush ('Collecting bg trials...')
        bg_tsds = bk.get_all (
                '{0}/bg_tsds'.format (self.mode_dir),
                '*.tsdist')
        saving (bg_tsds, '{0}/bg_tsds.dict'.format (self.mode_dir))
        return bg_tsds


    @property
    def sig_info (self):
        try:
            return self._sig_info
        except:
            try:
                return loading ('{}/sig_info.dict'.format (self.mode_dir))
            except:
                return self.collect_sig_info ()

    @command
    def collect_sig_info (self):
        """Collect signal injection results and cache."""
        prush ('Collecting signal results...')
        d = bk.get_all (
                '{}/n_sig'.format (self.mode_dir),
                'sens.pickle',
                lambda x: x[0])
        saving (d, '{}/sig_info.dict'.format (self.mode_dir))
        return d


    @command
    def setup_data (self):
        """Set up `data`, `nu`, and `pdfs`."""
        if 'txt' in self.mode.words:
            xd = np.genfromtxt (
                '{}/data/IC86-I_data.txt'.format (self.root_dir), names=True)
            data = arrays.Arrays (dict (
                (k,xd[k])
                for k in ('ra', 'dec', 'sigma', 'logE')))
            data.sinDec = np.sin (data.dec)
        else:
            xd = loading ('{}/data/exp.pickle'.format (self.root_dir))
            data = arrays.Arrays (dict (
                (k,xd[k])
                for k in ('ra', 'sinDec', 'sigma', 'logE')))

        if 'txt' in self.mode.words:
            xn = np.genfromtxt (
                '{}/data/IC86-I_MC.txt'.format (self.root_dir), names=True)
            nu = arrays.Arrays (dict (
                (k,xn[k])
                for k in ('ra', 'dec', 'sigma', 'logE',
                        'trueRa', 'trueDec', 'trueE', 'ow')))
            nu.sinDec = np.sin (nu.dec)
        else:
            xn = loading ('{}/data/MC.pickle'.format (self.root_dir))
            nu = arrays.Arrays (dict (
                (k,xn[k])
                for k in ('ra', 'sinDec', 'sigma', 'logE',
                        'trueRa', 'trueDec', 'trueE', 'ow')))
        xzenith, xazimuth = coord.rotate_source_to_xaxis (
            nu.trueDec+pi/2, nu.trueRa, np.arccos (-nu.sinDec), nu.ra)
        xazimuth[xazimuth < pi] += 2*pi
        xazimuth[xazimuth > pi] -= 2*pi
        nu.xaxis_zenith = xzenith
        nu.xaxis_azimuth = xazimuth

        data.apply_cut ((1 <= data.logE) & (data.logE < 10))
        nu.apply_cut ((1 <= nu.logE) & (nu.logE < 10))

        data_dir = misc.ensure_dir ('{}/data'.format (self.mode_dir))
        saving (data, '{}/data.arrays'.format (data_dir))
        saving (nu, '{}/nu.arrays'.format (data_dir))

        pdfs = pdf.PDFs (
            pdf.BgSpacePDF (
                -data.sinDec,
                bins=10, range=(-1,1), fit=False
            ),
            pdf.EnergyPDF (
                -data.sinDec, data.logE, weights=None,
                bins=20, range=((-1,1),(1,10)), fit=False
            ),
            pdf.EnergyPDFs (
                -nu.sinDec, nu.logE, np.log10 (nu.trueE), nu.ow * nu.trueE**-2,
                np.arange (1, 4.1, .25),
                bins=20, range=((-1,1),(1,10)), fit=False
            )
        )
        saving (pdfs, '{}/pdfs.pdfs'.format (data_dir))


    @command
    def submit_do_bg_ts (self):
        """Submit jobs for bg-only TSDists."""
        job_root = self.conf.root.job_dir
        job_dir = '{}/do_bg_ts/{}'.format (job_root, job_id)
        s = submitter.Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (os.path.abspath, self.conf.root_dirs))

        for zenith_deg in np.arange (1, 180, 2.):
            for i_job in xrange (int (self.opts.n_jobs)):
                command = '$IR4 {} do_bg_ts --conf-dirs={}' \
                    ' --n-trials={}' \
                    ' --zenith={:07.3f}' \
                    ' --seed={}'.format (
                        this_script, confs,
                        self.opts.n_trials, zenith_deg, i_job
                    )
                label = 'do_bg_ts__zen_{:07.3f}__seed_{}'.format (
                    zenith_deg, i_job
                )
                commands.append (command)
                labels.append (label)

        s.submit_condor00 (commands, labels, blacklist=self.blacklist)
        pjobdir (job_dir)


    @command
    def do_bg_ts (self):
        """Build a bg-only TSDist."""

        seed = self.opts.seed
        np.random.seed (seed)

        n_trials = self.opts.n_trials
        zenith_deg = self.opts.zenith
        zenith = np.radians (zenith_deg)

        tsd = trial.get_tsdist (n_trials, zenith, 0,
                                self.pdfs, self.bg(),
                                log_frac=100 / n_trials)

        sm = bk.SavingModel (
                'zenith_deg/tsd',
                '{:07.3f}/{:08d}.tsdist',
                )
        filename = sm.save (
                tsd, '{}/bg_tsds'.format (self.mode_dir),
                zenith_deg, seed)
        prush ('->', filename)


    @command
    def submit_do_n_sig (self):
        """Submit jobs for bg-only TSDists."""
        job_root = self.conf.root.job_dir
        job_dir = '{}/do_n_sig/{}'.format (job_root, job_id)
        s = submitter.Submitter (job_dir=job_dir)
        commands, labels = [], []

        this_script = os.path.abspath (__file__)
        confs = ','.join (map (os.path.abspath, self.conf.root_dirs))

        for zenith_deg in np.arange (0, 181, 7.5):
            command = '$IR4 {} do_n_sig --conf-dirs={}' \
                ' --n-trials={}' \
                ' --zenith={:07.3f}' \
                ' --gamma=2' \
                ' --sigma={}' \
                ' --beta={}' \
                ' --seed=0'.format (
                    this_script, confs,
                    self.opts.n_trials, zenith_deg,
                    self.opts.sigma, self.opts.beta
                )
            label = 'do_n_sig__zen_{:07.3f}__sigma_{:.1f}__beta_{:.2f}'.format (
                zenith_deg, self.opts.sigma, self.opts.beta
            )
            commands.append (command)
            labels.append (label)

        s.submit_condor00 (commands, labels, blacklist=self.blacklist)
        pjobdir (job_dir)


    @command
    def do_n_sig (self):
        """Calculate sensitivity and discovery potential fluxes."""

        seed = self.opts.seed
        np.random.seed (seed)

        n_trials = self.opts.n_trials
        zenith_deg = self.opts.zenith
        zenith = np.radians (zenith_deg)
        sigma = self.opts.sigma
        beta = self.opts.beta
        gamma = self.opts.gamma

        b = self.bg ()
        s = self.sig (zenith)
        pdfs = self.pdfs

        bg_tsd = bk.get_best (self.bg_tsds, zenith_deg)
        fit = sigma > 2
        ts = bg_tsd.sigma_thresh (sigma, fit=fit)

        prush ()
        prush ('Calculating n_sig...')
        prush ('- zenith = {0:.3f} deg'.format (zenith_deg))
        prush ('- spectrum ~ E ^ (- {0:.3f} )'.format (gamma))
        prush ('- n_sigma = {0:.2f} %'.format (sigma))
        prush ('- beta = {0:.3f} %'.format (beta * 100))
        prush ('- ts > {0:.5f}'.format (ts))
        prush ()

        result = trial.get_n_sig(ts, beta,
                                s.source_zenith, s.source_azimuth,
                                pdfs, b, s,
                                n_trials=n_trials, tol=0.01,
                                log=True, full_output=True)
        n_sig = result['n_sig'][-1]
        sens = s.to_flux (n_sig, 1e-3 / (321 * 86400))

        prush ('Obtained Phi0 = {:.4e} (n_sig = {:.4f})'.format (sens, n_sig))

        sm = bk.SavingModel (
            'sigma/beta/gamma/zenith',
            '{:.1f}/{:3.1f}/{:4.2f}/{:07.3f}/sens.pickle'
        )
        filename = sm.save (
            (sens, n_sig, s.n_exp),
            '{}/n_sig'.format (self.mode_dir),
            sigma, beta, gamma, zenith_deg,
        )
        prush ('->', filename)


    @command
    def sd (self):
        """Plot sensitivity and discovery potential."""

        misc.tex_mpl_rc ()

        sig_info = self.sig_info

        fig = getfig (aspect=16/10., width=6)
        ax = fig.add_subplot (111)

        nfig = getfig (aspect=16/10., width=6)
        nax = nfig.add_subplot (111)

        rfig = getfig (aspect=16/10., width=6)
        rax = rfig.add_subplot (111)

        curves = {}
        for n_sigma in (0, 5):
            if n_sigma == 0:
                thing = 'Sensitivity'
                CL = 0.9
                ls = '--'
            else:
                thing = 'Disc. Pot.'
                CL = 0.5
                ls = '-'

            color='b'
            alpha=.8

            label = r'$E^{{-2}}$ {}'.format (thing)
            x = bk.get_best (
                sig_info, n_sigma, CL, 2
            )
            sin_dec = np.array ([
                np.cos (-k/180.*pi) for k in sorted (x)]
            )
            curve = np.array ([
                x[k][0] for k in sorted (x)
            ])
            ncurve = np.array ([
                x[k][1] for k in sorted (x)
            ])

            ax.semilogy (sin_dec, curve,
                         label=label, ls=ls, color=color, alpha=alpha, lw=2)
            nax.plot (sin_dec, ncurve,
                      label=label, ls=ls, color=color, alpha=alpha, lw=2)
            curves[n_sigma] = sin_dec, curve

        x, y = np.genfromtxt ('{}/etc/orig_sens.txt'.format (self.root_dir)).T
        label = r'$E^{{-2}}$ Sensitivity (original)'.format ()
        ax.semilogy (x, 1e-3 * y,
                     label=label, ls='--', color='k', alpha=.8, lw=1)
        rax.plot (curves[0][0],
                  curves[0][1] / np.interp (curves[0][0], x, 1e-3 * y),
                  ls='--', color='k', lw=1, label='Sensitivity ratio')

        x, y = np.genfromtxt ('{}/etc/orig_disc.txt'.format (self.root_dir)).T
        label = r'$E^{{-2}}$ Disc. Pot. (original)'.format ()
        ax.semilogy (x, 1e-3 * y,
                     label=label, ls='-', color='k', alpha=.8, lw=1)
        rax.plot (curves[5][0],
                  curves[5][1] / np.interp (curves[5][0], x, 1e-3 * y),
                  ls='-', color='k', lw=1, label='Disc. Pot. ratio')

        ax.set_xlabel (r'$\sin(\delta)$')
        nax.set_xlabel (r'$\sin(\delta)$')
        rax.set_xlabel (r'$\sin(\delta)$')
        ax.set_ylabel (r'$E^2 '
                '\cdot (E/100\,\mathrm{TeV})^{\gamma-2}'
                '\cdot dN/dE\,\,\,'
                '[\mathrm{TeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
        nax.set_ylabel (r'$n_\mathrm{inj}$')
        nax.set_ylabel (r'ratio')

        ax.grid ()
        legend = ax.legend (loc='upper right', prop=propsmall,
                            handlelength=4, ncol=2)
        frame = legend.get_frame ()
        frame.set_linewidth (0)

        nax.grid ()
        legend = nax.legend (loc='upper right', prop=propsmall,
                            handlelength=4, ncol=2)
        frame = legend.get_frame ()
        frame.set_linewidth (0)

        rax.set_ylim (0.5, 1.5)
        rax.grid ()
        legend = rax.legend (loc='best', prop=propsmall,
                             handlelength=4, ncol=2)
        legend.get_frame ().set_linewidth (0)

        fig.subplots_adjust (bottom=.14, top=.91, right=.97)
        nfig.subplots_adjust (bottom=.14, top=.91, right=.97)
        rfig.subplots_adjust (bottom=.14, top=.91, right=.97)

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        savingfig (fig, plot_dir, 'sensdisc')
        savingfig (nfig, plot_dir, 'sensdisc_ninj')
        savingfig (rfig, plot_dir, 'sensdisc_ratio')


    @command
    def tsds (self):
        """Plot TSDist's."""

        plot_dir = misc.ensure_dir ('{0}/plots/tsds'.format (self.mode_dir))

        tsds = self.bg_tsds

        for zenith_deg in sorted (tsds):
            tsd = tsds[zenith_deg]

            fig = getfig (aspect=16/10., width=6)
            ax = fig.add_subplot (111)
            ax.semilogy ()
            hl.plot1d (ax, tsd.get_hist (bins=40).normalize (integrate=True),
                       color='b')
            ts = np.linspace (1e-3, tsd.ts_values.max (), 100)
            ax.plot (ts, tsd.chi2.pdf (ts), color='.8', ls='--')
            ax.set_xlabel ('TS')
            ax.set_ylabel ('probability density')
            fig.subplots_adjust (bottom=.14, top=.91, right=.97)
            savingfig (fig, plot_dir,
                       'fscale_bg_tsd_zen_{:07.3f}'.format (zenith_deg))
            plt.close (fig)


    @command
    def aeff (self):
        """
        Plot and tabulate Aeff.
        """

        import colormaps as cmaps
        plt.register_cmap (name='viridis', cmap=cmaps.viridis)
        plt.set_cmap (cmaps.viridis)

        logEmin, logEmax = 2., 9.
        dlogE = 0.1
        n_bins_E = (logEmax - logEmin) / dlogE
        dcz = 0.01
        dOmega = 2 * pi * dcz
        n_bins_cz = 2 / dcz

        nu = self.nu
        nu.cz = -np.sin (nu.trueDec)
        w_aeff = 1 / (1e4 * np.log (10)) * nu.ow / nu.trueE / dOmega / dlogE

        h_aeff = hl.hist (
            (nu.trueE, nu.cz), w_aeff,
            bins=(n_bins_E, n_bins_cz),
            range=((10**logEmin, 10**logEmax), (-1, 1)),
            log=(True, False),
        )

        misc.tex_mpl_rc (True)
        fig = getfig (aspect=4/3., width=5)
        ax = fig.add_subplot (111)
        fig.subplots_adjust (bottom=.15, left=.15)
        result = hl.plot2d (ax, h_aeff, cbar=True, log=True,
                            vmin=5e-6, vmax=1e4, zmin=5e-6)
        result['colorbar'].set_label (r'effective area $[\text{m}^2]$')
        ax.set_xlabel ('neutrino energy [GeV]')
        ax.set_ylabel (r'$\cos(\text{zenith})$')

        plot_dir = misc.ensure_dir ('{0}/plots'.format (self.mode_dir))
        savingfig (fig, plot_dir, 'aeff')

        bins = h_aeff.bins
        filename = '{}/aeff.txt'.format (plot_dir)
        prush ('-> {} ...'.format (filename))
        with open (filename, 'w') as f:
            pr = lambda *a, **kw: print (*a, file=f, **kw)
            pr ('# {:>11}{:>13}{:>16}{:>16}{:>16}'.format (
                'E_min[GeV]', 'E_max[GeV]',
                'cos(zenith)_min', 'cos(zenith)_max',
                'Aeff[m^2]'
            ))
            for (Emin, Emax) in izip (bins[0][:-1], bins[0][1:]):
                for (czmin, czmax) in izip (bins[1][:-1], bins[1][1:]):
                    pr ('{:13.3e}{:13.3e}{:+16.2f}{:+16.2f}{:16.3e}'.format (
                        Emin, Emax,
                        czmin, czmax,
                        h_aeff.get_value (1.001 * Emin, 1e-3 + czmin)
                    ))



if __name__ == '__main__':
    app = TPS861 ()
    app.run ()
    try:
        __IPYTHON__
    except:
        pass
    else:
        try:
            data = app.data
            nu = app.nu
            pdfs = app.pdfs
            sig = app.sig
            bg = app.bg
        except:
            pass


