#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##
import os
import socket
import time

import numpy as np

from icecube.umdtools.submitter import Submitter

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())

# assign decs for the purposes of naming directories. I'll have to generalize this bin number to match the number of single source background files... #

bins = 11
degdecrange=np.linspace(-85.0,85.0,bins)
sindecrange=np.sin(np.radians(degdecrange))
#sindecrange=np.linspace(np.sin(np.radians(-85.0)),np.sin(np.radians(85.0)))
#degdecrange = np.degrees(np.arcsin(sindecrange))
#Note - I used to do it with equal spacing in degrees, but I think I want equal spacing in sin_dec, for the purpose of graphing.

# set up submitter
submitter = Submitter (job_dir='/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/jobs/{0}'.format (job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/calculate_sensitivity_single_source.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []

#Time to get each dec command to go after a separate dec folder of background trials.
for sin_dec, deg_dec in zip(sindecrange,degdecrange):
    commands.append('{0} {1} --dec {2}'.format (env_shell, job_script, sin_dec))
    labels.append ('handle_dec_{0:+010.5}'.format(deg_dec))


##Sumbit!
submitter.submit_npx4 (commands, labels)
