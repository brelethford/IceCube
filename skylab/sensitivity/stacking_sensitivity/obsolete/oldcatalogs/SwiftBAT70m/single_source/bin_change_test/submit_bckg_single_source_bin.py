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

# set up submitter
submitter = Submitter (job_dir='/data/user/brelethford/Output/all_sky_sensitivity/results/single_stacked/bin_test/jobs/{0}'.format (job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/bckg_single_source_bin.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []
## For this we'll need to do background trials in batches of 1000. So...
batchsize=1000
total_bckg_trials = 30000
batches = int(np.round(total_bckg_trials/batchsize))
bins = 11
degdecrange=np.linspace(-85.0,85.0,bins)
sindecrange=np.sin(np.radians(degdecrange))
#sindecrange=np.linspace(np.sin(np.radians(-85.0)),np.sin(np.radians(85.0)),bins)
#degdecrange = np.degrees(np.arcsin(sindecrange))

for sin_dec, deg_dec in zip(sindecrange,degdecrange):
  for i in range(batches):
    commands.append('{0} {1} --dec {2} --batch {3} --batchsize {4}'.format (env_shell, job_script, sin_dec, i, batchsize))
    labels.append ('handle_dec_{0:+010.5}_batch'.format(deg_dec)+str(i))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
