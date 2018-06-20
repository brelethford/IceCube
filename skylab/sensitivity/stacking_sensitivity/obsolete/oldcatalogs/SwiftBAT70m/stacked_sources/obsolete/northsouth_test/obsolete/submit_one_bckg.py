#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##

#I want to test how well I've reassembled a buncha background trials, so this one is just gonna do 1 batch of 30000.
import os
import socket
import time
from optparse import OptionParser
import numpy as np

from icecube.umdtools.submitter import Submitter

# establish weighting scheme for searching (only llh_model needs weighting for bckg)

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--sky', dest = 'sky', type = str,
                default = None, metavar = 'SKY',
                help = 'tells which sources to use (and which folders to reference)')

opts, args = parser.parse_args ()
sky = opts.sky

if sky == None:
  raise NameError('Error - did not establish which part of the sky to test.')

print ("The sky is: {}".format (sky))

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())

# set up submitter
submitter = Submitter (job_dir='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/jobs/{}'.format (job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/one_bckg.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []
## For this we'll need to do background trials in batches of 1000. So...
batchsize=30000
total_bckg_trials = 30000

batches = int(np.round(total_bckg_trials/batchsize))

for i in range(batches):
    commands.append('{0} {1} --batch {2} --batchsize {3} --sky {4}'.format (env_shell, job_script, i, batchsize, sky))
    labels.append ('handle_batch_'+str(i))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
