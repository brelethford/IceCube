###############################################################################
# @file   config.py
# @author Josh Wood
# @date   Oct 13 2017
# @brief  Sets up skylab likelihood and injector for 7yr point source anlaysis
#         done by Stefan Coenders.
###############################################################################

import os
import numpy as np 
from optparse import OptionParser
import argparse
from   skylab.ps_llh            import PointSourceLLH, MultiPointSourceLLH, FastMultiPointSourceLLH
#from   skylab_master.psLLH            import PointSourceLLH, MultiPointSourceLLH
from   skylab.ps_injector       import PointSourceInjector

from   skylab.llh_models        import ClassicLLH, EnergyLLH
#from   skylab_master.ps_model        import ClassicLLH, EnergyLLH
from   skylab.datasets          import Datasets
from   skylab.sensitivity_utils import DoubleChiSquare

###############################################################################

###########
# GLOBALS #
###########

# energy units
GeV = 1
TeV = 1000*GeV

###############################################################################




def config(seasons, seed = 1, scramble = True, verbose = True, dec = 0.):
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

  # store individual llh as lists to prevent pointer over-writing
  llh = []

  multillh = FastMultiPointSourceLLH(seed = seed, ncpu = 25)

  # setup likelihoods
  if verbose: print("\n seasons:")
  for season in np.atleast_1d(seasons):

    #sample = 'PointSourceTracks7yr_galactic_plane'
    #sample = 'PointSourceTracks7yr'
    #sample = 'PointSourceTracks_v002p01b'
    sample = 'PointSourceTracks'

    exp, mc, livetime = Datasets[sample].season(season)
    sinDec_bins = Datasets[sample].sinDec_bins(season)
    energy_bins = Datasets[sample].energy_bins(season)

    mc = mc[mc["logE"] > 1.]

    if verbose: print("   - %-15s (livetime %7.2f days, %6d events)" % (season, livetime, exp.size))

    llh_model = EnergyLLH(twodim_bins = [energy_bins, sinDec_bins], allow_empty=True,
                          bounds = [1., 4.], seed = 2., kernel=1)

    llh.append( PointSourceLLH(exp, mc, livetime, mode = "box", seed = seed, 
                               scramble = scramble, llh_model = llh_model,
                               nsource_bounds=(0., 1e3), delta_ang=np.deg2rad(1.), nsource=15.) )

    multillh.add_sample(season, llh[-1])

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

  '''
  inj = PointSourceInjector(gamma = 2., E0 = 1*TeV, seed = seed)
  inj.fill(dec, multillh.exp, multillh.mc, multillh.livetime)
  '''
  inj = None
  ############
  # INJECTOR #
  ############

  #############################################################################

  '''
  if verbose:
      print("\n fitted spectrum:")
      vals = (inj.E0/TeV)
      print("   - dN/dE = A (E / %.1f TeV)^-index TeV^-1cm^-2s^-1" % vals)
      print("   - index is *fit*")

      print("\n injected spectrum:")
      vals = (inj.E0/TeV, inj.gamma)
      print("   - dN/dE = A (E / %.1f TeV)^-%.2f TeV^-1cm^-2s^-1" % vals)
  '''

  return (multillh, inj)

###############################################################################
