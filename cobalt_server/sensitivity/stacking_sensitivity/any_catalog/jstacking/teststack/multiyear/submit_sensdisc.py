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

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = None, metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = None, metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--n', dest = 'n', type = int,
                default = 4, metavar = 'N',
                help = 'Number of years tried.')


opts, args = parser.parse_args ()
catalog = opts.catalog
llhweight = opts.llhweight
injweight = opts.injweight
n_year = opts.n

try:
   print ("Catalog used: " + catalog)
except:
   print ('No catalog selected! Please select a catalog.')

if llhweight is None:
        if catalog == 'teststack':
                llhweight = 'equal'
        elif catalog == 'teststack50':
                llhweight = 'equal'
        elif catalog == 'teststack300':
                llhweight = 'equal'
        else:
                raise ValueError('catalog name not found')

if injweight is None:
        injweight = llhweight


print ("The llhweight is: {}".format (llhweight))
print ("The injweight is: {}".format (injweight))



# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())

jobdir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/jstacking_{1}yr/jobs/'.format (catalog,n_year))
# set up submitter
submitter_sens = Submitter (job_dir=jobdir+'sens_{0}_{1}'.format (job_id, llhweight))
submitter_disc = Submitter (job_dir=jobdir+'disc_{0}_{1}'.format (job_id, llhweight))

submitter_sens.memory = 6
submitter_disc.memory = 6

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/calculate_sensdisc.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands_sens, commands_disc, labels_sens, labels_disc = [], [], [], []

## For this we'll need to do background trials in batches of 1000. So...
bins=1
gammarange = np.linspace(2.0,2.0,bins)

for gamma in gammarange:
    opts_list ='{0} {1} --cat {2} --llhweight {3} --injweight {4} --gamma {5}'.format (env_shell, job_script, catalog, llhweight, injweight, gamma)
    commands_sens.append(opts_list)
    labels_sens.append ('handle_sensitivity_gamma_{}'.format(gamma))
    commands_disc.append(opts_list + ' --disc 1')
    labels_disc.append ('handle_disc_gamma_{}'.format(gamma))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter_sens.submit_npx4 (commands_sens, labels_sens)
submitter_disc.submit_npx4 (commands_disc, labels_disc)

