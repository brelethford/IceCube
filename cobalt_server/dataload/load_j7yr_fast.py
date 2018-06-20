# -*-coding:utf8-*-
from __future__ import print_function
import sys
from skylab.datasets import Datasets
r"""Load different years of data in and load it correctly to the LLH classes. Specifically, I'm going to try to use the Dataset object."""

#Loading Zone#
sys.path.append("/home/brelethford/Documents/IceCube_Research/Scripts/AGN_Core/dataload/coenders_config")
import config

#Note - no sirin, no mese.

def load1yr_86():
    seasons = ['IC86, 2011']
    llh,inj  = config.config(seasons=seasons)
    return (llh,inj)

def load1yr():
    seasons = ['IC40']
    llh,inj  = config.config(seasons=seasons)
    return (llh,inj)

def load4yr():
    seasons = ['IC40','IC59','IC79','IC86, 2011']
    llh,inj  = config.config(seasons=seasons)
    return (llh,inj)

def load7yr():
    seasons = ['IC40','IC59','IC79','IC86, 2011','IC86, 2012-2014']
    llh,inj  = config.config(seasons=seasons)
    return (llh,inj)

