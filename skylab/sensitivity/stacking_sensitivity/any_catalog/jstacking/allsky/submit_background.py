#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##
import os
import socket
import time
from optparse import OptionParser
import numpy as np
from icecube.umdtools import misc
from icecube.umdtools.submitter import Submitter

# establish weighting scheme for searching (only llh_model needs weighting for bckg)

parser = OptionParser (usage = '%prog [options]')

parser.add_option ('--sirin', dest = 'sirin', 
                default = 0, action='store_true', metavar = 'SIRIN',
                help = 'Determines if unsplined IC79 is used')

parser.add_option ('--n', dest = 'n', type = int,
                default = '4', metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--mese', dest = 'mese',
                default = 0, action='store_true', metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
n = opts.n
sirin = opts.sirin
mese = opts.mese

if mese:
  print ("Using MESE data")
  mese=1
if sirin:
  print ("Using sirin's IC79")
  sirin=1

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())

jobdir = misc.ensure_dir ('/data/user/brelethford/Output/all_sky_sensitivity/mstacking_{0}yr/jobs/'.format (str(n)))
# set up submitter
submitter = Submitter (job_dir=jobdir+'{0}'.format (job_id))
if n>4:
  submitter.memory = 4

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/background_trials.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []
## For this we'll need to do background trials in batches of 1000. So...
batchsize=3000
total_bckg_trials = 30000

batches = int(np.round(total_bckg_trials/batchsize))
bins=11
degdecrange=np.linspace(-85.0,85.0,bins)
sindecrange=np.sin(np.radians(degdecrange))
for i in range(batches):
  for sindec,deg_dec in zip(sindecrange,degdecrange):
    commands.append('{0} {1} --batch {2} --batchsize {3} --sindec {4} --years {5} --sirin {6} --mese {7}'.format (env_shell, job_script, i, batchsize, sindec, str(n), sirin, mese))
    labels.append ('handle_dec_{0:+010.5}_batch_{1}'.format(deg_dec,i))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
