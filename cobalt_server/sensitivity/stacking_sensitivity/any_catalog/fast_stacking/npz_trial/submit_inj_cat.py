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

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG', default = None,
                help = 'Selects the catalog of sources.')

parser.add_option ('--weight', dest = 'weight', type = str,
                default = None, metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--gamma', dest = 'gamma', type = float,
                default = 2.0, metavar = 'Gamma',
                help = 'Set injection spectrum.')

parser.add_option ('--n', dest = 'n',
                default = 4, metavar = 'N',
                help = 'Number of years of data')

opts, args = parser.parse_args ()
catalog = opts.catalog
weight = opts.weight
gamma = opts.gamma
n = int(opts.n)

if n == 1:
  sample = "PointSourceTracks1yr"
elif n == 4:
  sample = "PointSourceTracks4yr"
elif n == 7:
  sample = "PointSourceTracks7yr"

print ("Using sample {}".format(sample))

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
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/npz/stacking_sensitivity/{0}yr/{1}/jobs/'.format (str(n), catalog))
else:
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/npz/allsky_sensitivity/{0}yr/jobs/'.format (str(n)))
submitter = Submitter (job_dir=jobdir+'{0}_{1}'.format (job_id, weight))


if catalog == 'SwiftBAT70m':
  if int(n) == 7:
    submitter.memory = 30
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
job_script = this_dir + '/n_inj_catalog.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
batchsize=50
total_inj_trials = 1000
batches = int(np.round(total_inj_trials/batchsize))
#So, 100 trials at each injection spot - of which there will be 11.

#I'm eschewing estimating a sens and disc range (we'll do 0-50 sens, 0-200 disc...
if catalog in ['Milagro17']:
  n_sens_bounds = [0.,50.]
  n_disc_bounds = [50., 201.]
  nstep_sens = 2
  nstep_disc = 25 
elif catalog in ['SNR_noPWN','SNR_PWN','SNR_cloud']:
  n_sens_bounds = [0.,16.]
  n_disc_bounds = [16., 61.]
  nstep_sens = 2
  nstep_disc = 5
else:
  n_sens_bounds = [0.,100.]
  n_disc_bounds = [100., 401.]
  nstep_sens = 10
  nstep_disc = 25

ninjs_sens = np.arange(n_sens_bounds[0], n_sens_bounds[1], nstep_sens)
ninjs_disc = np.arange(n_disc_bounds[0], n_disc_bounds[1], nstep_disc)
ninjs = np.concatenate((ninjs_sens,ninjs_disc)) 

commands, labels, = [], []

for ninj in ninjs:
  for batch in range(batches):
    opts_list = '{0} {1} --cat {2} --batch {3} --nscramble {4} --ninj {5} --seed {6} --weight {7} --sample {8} --gamma {9}'.format (env_shell, job_script, catalog, batch, batchsize, ninj, np.random.randint(1e8), weight, sample, gamma)
    commands.append(opts_list)
    labels.append ('handle_sensdisc_{}_n_sig_{}_batch_{}'.format(catalog,ninj,batch))

# submit! This will return several batches of injection events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
