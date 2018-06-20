#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##
import os
import socket
import time
import numpy as np
from optparse import OptionParser 

from icecube.umdtools.submitter import Submitter

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--sky', dest = 'sky', type = str,
                default = 'None', metavar = 'SKY',
                help = 'tells which sources to use (and which folders to reference)')

opts,args = parser.parse_args ()
sky = opts.sky

if sky == None:
  raise NameError('Error - did not establish which part of the sky to test.')



print ("The skytest is: {}".format (sky))

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}_{3}'.format (
        socket.gethostname (), os.getpid (), time.time (), sky)

# set up submitter
submitter = Submitter (job_dir='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/northsouth/jobs/{0}_box'.format ( job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/northsouth_test_box.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []

commands.append('{0} {1} --sky {2}'.format (env_shell, job_script, sky))
labels.append ('handle_{}'.format(sky))

submitter.submit_npx4 (commands, labels)
