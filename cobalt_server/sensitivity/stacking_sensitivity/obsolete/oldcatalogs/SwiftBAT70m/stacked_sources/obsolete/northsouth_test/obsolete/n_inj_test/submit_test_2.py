#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##
import os
import socket
import time
import numpy as np
from optparse import OptionParser 

from icecube.umdtools.submitter import Submitter

parser = OptionParser (usage = '%prog [options]')
parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = 'uniform', metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = 'uniform', metavar = 'INJWEIGHT',
                help = 'Sets the weighting used for signal injection.')

opts, args = parser.parse_args ()
llhweight = opts.llhweight
injweight = opts.injweight

print ("The llhweight is: {}".format (llhweight))
print ("The injweight is: {}".format (injweight))

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}_{3}inj'.format (
        socket.gethostname (), os.getpid (), time.time (), injweight)

# set up submitter
submitter = Submitter (job_dir='/data/user/brelethford/Output/stacking_sensitivity/SwiftBAT70m/{0}/jobs/{1}'.format (llhweight, job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/n_inj_test_2.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []
## However, I don't know how I'd split this up, so I'm just gonna hope that I can do it all at once...

bins=1
gammarange = np.linspace(2.0,2.0,bins)

for gamma in gammarange:
    commands.append('{0} {1} --llhweight {2} --injweight {3} --gamma {4}'.format (env_shell, job_script, llhweight, injweight, gamma))
    labels.append ('handle_sensitivity_gamma_{}'.format(gamma))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
