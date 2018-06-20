#!/usr/bin/env python
##The purpose of this script is to use submitter to simultaneously find many background events for the same weighting algorithm: this one is for uniform weights.##
import os
import socket
import time
import numpy as np
from optparse import OptionParser 
from icecube.umdtools import misc
from icecube.umdtools.submitter import Submitter

parser = OptionParser (usage = '%prog [options]')

catalog = 'SwiftBAT70m'
weight = 'equal'
n=int(7)

if n == 1:
  sample = "PointSourceTracks1yr"
elif n == 4:
  sample = "PointSourceTracks4yr"
elif n == 7:
  sample = "PointSourceTracks7yr"

if catalog:
  if weight is None:
        if catalog == 'SwiftBAT70m':
                weight = 'uniform'
        elif catalog == '4yr_Starburst':
                weight = 'flux'
        elif catalog == 'teststack' or catalog == 'teststack50' or catalog == 'teststack300':
                weight = 'equal'
        elif catalog == 'Milagro17':
                weight = 'weight'
        elif catalog == 'blackhole':
                weight = 'flux2m'
        elif catalog == 'SNR_noPWN' or catalog == 'SNR_cloud' or catalog == 'SNR_PWN':
                weight = 'weight'
        elif catalog == 'WHSP_Blazars':
                weight = 'weights'
        elif catalog == '2LAC':
                weight = 'flux'
        else:
                raise ValueError('catalog name not found')

print ("The weight is: {}".format (weight))

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())
# set up submitter
if catalog:
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/fast/stacking_sensitivity/{0}yr/{1}/jobs/'.format (str(n), catalog))
else:
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/fast/allsky_sensitivity/{0}yr/jobs/'.format (str(n)))
submitter = Submitter (job_dir=jobdir+'{0}_{1}'.format (job_id, weight))


if catalog == 'SwiftBAT70m':
  if int(n) == 7:
    submitter.memory = 22 
  else:
    submitter.memory = 20
elif catalog == 'teststack300':
  submitter.memory = 20
elif catalog == 'blackhole':
  submitter.memory = 14
else:
  submitter.memory = 8

if not catalog:
  submitter.memory = 2


# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/background_test.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
batchsize=1000
total_bckg_trials = 1000

batches = int(np.round(total_bckg_trials/batchsize))

bins = 35
degdecrange=np.linspace(-85.0,85.0,bins)
sindecrange=np.sin(np.radians(degdecrange))

commands, labels, = [], []

if catalog:
  for batch in range(batches):
    opts_list = '{0} {1} --cat {2} --batch {3} --nscramble {4} --seed {5} --weight {6} --sample {7}'.format (env_shell, job_script, catalog, batch, batchsize, np.random.randint(1e8), weight, sample)
    commands.append(opts_list)
    labels.append ('handle_sensdisc_{}_{}'.format(catalog,batch))
else:
  for sindec,deg_dec in zip(sindecrange,degdecrange):
    opts_list = '{0} {1} --sindec {2} --sample {3}'.format (env_shell, job_script, sindec, sample)
    commands.append(opts_list)
    labels.append ('handle_sensdisc_{}'.format(np.int(deg_dec)))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
