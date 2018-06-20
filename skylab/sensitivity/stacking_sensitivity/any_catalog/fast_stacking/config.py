###############################################################################
# @file   config.py
# @author Josh Wood
# @date   Oct 13 2017
# @brief  Sets up skylab likelihood and injector for 7yr point source anlaysis
#         done by Stefan Coenders.
###############################################################################

import os
import numpy as np 

from   skylab.ps_llh            import StackingPointSourceLLH, PointSourceLLH, MultiPointSourceLLH
from   skylab.ps_injector       import PointSourceInjector

from   skylab.llh_models        import ClassicLLH, EnergyLLH
from   skylab.datasets          import Datasets
from   skylab.sensitivity_utils import DeltaChiSquare

import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

if int(mpl.__version__[0])>1:
  mpl.style.use('classic')

###############################################################################

###########
# GLOBALS #
###########

# energy units
GeV = 1
TeV = 1000*GeV

def TXS_location():
    src_ra  = np.radians(77.358)
    src_dec = np.radians(5.693)
    return(src_ra, src_dec)

###############################################################################

def config(seasons, seed = 1, scramble = True, e_range=(0,np.inf),
           verbose = True, gamma = 2.0, dec = 0., remove = False, src_w = None):
  r""" Configure multi season point source likelihood and injector. 

  Parameters
  ----------
  seasons : list
    List of season names
  seed : int
    Seed for random number generator

  Returns
  -------
  multillh : MultiPointSourceLLH
    Multi year point source likelihood object
  inj : PointSourceInjector
     Point source injector object
  """

  print("Scramble is %s" % str(scramble))

  # store individual llh as lists to prevent pointer over-writing
  llh = []

  multillh = MultiPointSourceLLH(seed = seed, ncpu = 25)

  # setup likelihoods
  if verbose: print("\n seasons:")
  for season in np.atleast_1d(seasons):

    sample = season[0]
    name   = season[1]

    exp, mc, livetime = Datasets[sample].season(name)
    sinDec_bins = Datasets[sample].sinDec_bins(name)
    energy_bins = Datasets[sample].energy_bins(name)

    if sample == "GFU" and remove:
      exp = Datasets['GFU'].remove_ev(exp, 58018.87118560489) # remove EHE 170922A
      

    mc = mc[mc["logE"] > 1.]

    if verbose:
      print("   - % 15s : % 15s" % (sample, name))
      vals = (livetime, exp.size, min(exp['time']), max(exp['time']))
      print("     (livetime %7.2f days, %6d events, mjd0 %.2f, mjd1 %.2f)" % vals)

    llh_model = EnergyLLH(twodim_bins = [energy_bins, sinDec_bins], allow_empty=True,
                          bounds = [1., 4.], seed = 2., kernel=1)

    llh.append( StackingPointSourceLLH(exp, mc, livetime, mode = "box", 
                               scramble = scramble, llh_model = llh_model,
                               nsource_bounds=(0., 1e3), nsource=15.) )

    multillh.add_sample(sample+" : "+name, llh[-1])

    # save a little RAM by removing items copied into LLHs
    del exp, mc

  # END for (season)

  #######
  # LLH #
  #######

  #############################################################################

  ############
  # INJECTOR #
  ############

  inj = PointSourceInjector(gamma = gamma, E0 = 1*TeV, seed = seed, e_range=e_range)
  inj.fill(dec, multillh.exp, multillh.mc, multillh.livetime, src_w)

  ############
  # INJECTOR #
  ############

  #############################################################################

  if verbose:
      print("\n fitted spectrum:")
      vals = (inj.E0/TeV)
      print("   - dN/dE = A (E / %.1f TeV)^-index TeV^-1cm^-2s^-1" % vals)
      print("   - index is *fit*")

      print("\n injected spectrum:")
      vals = (inj.E0/TeV, inj.gamma)
      print("   - dN/dE = A (E / %.1f TeV)^-%.2f TeV^-1cm^-2s^-1" % vals)

  return (multillh, inj)

###############################################################################

def fit_background(trials, image):
  r""" Fit background only TS distribution using delta function
       at ts = 0 and chi squared function for ts > 0.

  Parameters
  ----------
  trials : numpy.ndarray
    trials array from BaseLLH.do_trials() 
  image : str
    path to output image

  Return
  ------
  median_ts : float
    median of TS distribution
  eta : float
    fraction TS > 0
  ndf : float
    effective chi squared degrees of freedom ts > 0
  scale : float
    effective scaling of TS > 0
  """

  ################################
  # HISTOGRAM BACKGROUND ONLY TS #
  ################################

  # fix rounding errors
  trials['TS'][ trials['TS'] < 0 ] = 0.

  # histogram
  y, bins = np.histogram(trials['TS'], np.arange(0., 20.01, .25))
  x  = (bins[1:] + bins[:-1])/2
  w  = (bins[1:] - bins[:-1])

  # see sensitivity_utils.py for definition of DeltaChiSquare()
  func = DeltaChiSquare()
  func.fit(x,y,w,verbose=True)

  # parameters describing TS > 0
  ndf   = func.ndf
  eta   = func.eta
  scale = func.scale

  # median TS for background only
  median_ts = np.percentile(trials['TS'], 50)

  ###########################
  # PLOT BACKGROUND ONLY TS #
  ###########################

  plt.bar(x-0.5*w, y, w, color='b', alpha=0.5, linewidth=0, label=("%d scrambles" % trials.size))
  plt.plot(x, trials.size*func.binval(x,w), color='r', linestyle='--', label="$\chi^{2}$ fit")

  ax = plt.gca()
  ax.set_yscale("log")
  ax.set_xlabel("TS", horizontalalignment='right', x=1.0)
  ax.set_ylabel("Number of Scrambles", horizontalalignment='right', y=1.0)
  ax.set_xlim(0, 20)
  ax.set_ylim(0.5, 1.1*trials.size)

  plt.legend(loc=2)
  plt.savefig(image)
  plt.clf()

  print ("\nSaved " + image)

  return [median_ts, eta, ndf, scale]

