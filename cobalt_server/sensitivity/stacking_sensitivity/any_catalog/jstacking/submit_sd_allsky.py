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

parser.add_option ('--n', dest = 'n',
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--sirin', dest = 'sirin',
                default = 0, action='store_true', metavar = 'SIRIN',
                help = 'Determines if unsplined IC79 is used')

parser.add_option ('--datatype', dest = 'datatype',
                default = 'standard', metavar = 'DATA',
                help = 'Determines if using npz or mydata')

parser.add_option ('--mese', dest = 'mese',
                default = 0, action='store_true', metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
n = opts.n
datatype = opts.datatype
sirin = opts.sirin
mese = opts.mese


if mese:
  print ("Using MESE data")
  mese=1
if sirin:
  print ("Using sirin's IC79 dataset")
  sirin = 1
else:
  print ("Using kai's IC79 dataset")

# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())
# set up submitter
jobdir = misc.ensure_dir ('/data/user/brelethford/Output/{0}/allsky_sensitivity/{1}yr/jobs/'.format (datatype,str(n)))
submitter_sens = Submitter (job_dir=jobdir+'sens_{0}'.format (job_id))
submitter_disc = Submitter (job_dir=jobdir+'disc_{0}'.format (job_id))

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/calculate_sensdisc.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands_sens, commands_disc, labels_sens, labels_disc = [], [], [], []

submitter_sens.memory = 6
submitter_disc.memory = 6

bins=35
degdecrange=np.linspace(-85.0,85.0,bins)
sindecrange=np.sin(np.radians(degdecrange))

gammabins=1
gammarange = np.linspace(2.0,2.0,gammabins)

for gamma in gammarange:
    for sindec,deg_dec in zip(sindecrange,degdecrange):
        opts_list = '{0} {1} --n {2} --gamma {3} --sirin {4} --mese {5} --sindec {6} --datatype {7}'.format (env_shell, job_script, n, gamma, sirin, mese, sindec, datatype)
        commands_sens.append(opts_list)
        labels_sens.append ('handle_sensitivity_gamma_{0}_dec_{1:+010.5}'.format(gamma,deg_dec))
        commands_disc.append(opts_list + ' --disc 1')
        labels_disc.append ('handle_disc_gamma_{0}_dec_{1:+010.5}'.format(gamma,deg_dec))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter_sens.submit_npx4 (commands_sens, labels_sens)
submitter_disc.submit_npx4 (commands_disc, labels_disc)
