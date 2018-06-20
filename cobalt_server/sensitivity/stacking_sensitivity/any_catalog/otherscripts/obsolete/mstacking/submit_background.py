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

parser.add_option ('--fom', dest = 'fom', type = float,
                default = 0.0, metavar = 'FOM',
                help = 'Sets the figure of merit cut for WHSP blazars.')

parser.add_option ('--sirin', dest = 'sirin', 
                default = 0, action='store_true', metavar = 'SIRIN',
                help = 'Determines if unsplined IC79 is used')

parser.add_option ('--mhubertest', dest = 'mhubertest',
                default = 0, action='store_true',metavar = 'MHUBERTEST',
                help = 'toggles mhuber code priority')

parser.add_option ('--n', dest = 'n', type = int,
                default = '4', metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--mese', dest = 'mese',
                default = 0, action='store_true', metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
catalog = opts.catalog
llhweight = opts.llhweight
injweight = opts.injweight
n = opts.n
sirin = opts.sirin
fom = opts.fom
mhubertest = opts.mhubertest
mese = opts.mese

if mese:
  print ("Using MESE data")
  mese=1
if mhubertest:
  print ("Using mhubertest")
  mhubertest=1
if sirin:
  print ("Using sirin's IC79")
  sirin=1

try:
   print ("Catalog used: " + catalog)
except:
   print ('No catalog selected! Please select a catalog.')

if llhweight is None:
        if catalog == 'SwiftBAT70m':
                llhweight = 'uniform'
        elif catalog == '4yr_Starburst':
                llhweight = 'flux'
        elif catalog == 'SNR_cloud' or catalog == 'SNR_noPWN' or catalog == 'SNR_PWN':
                llhweight = 'weight'
        elif catalog == 'WHSP_Blazars':
                llhweight = 'weights'
        elif catalog == '2LAC':
                llhweight = 'flux'
        else:
                raise ValueError('catalog name not found')

if injweight is None:
        injweight = llhweight

print ("The llhweight is: {}".format (llhweight))
print ("The injweight is: {}".format (injweight))



# get unique job id
job_id = '{0}_nixtime_{2:.0f}_job_{1}'.format (
        socket.gethostname (), os.getpid (), time.time ())
if mese:
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr_mese/jobs/'.format (catalog, str(n)))
else:
  jobdir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/jobs/'.format (catalog, str(n)))
# set up submitter
submitter = Submitter (job_dir=jobdir+'{0}_{1}'.format (job_id, llhweight))
if n>4:
  submitter.memory = 6
else:
  if catalog == 'WHSP_Blazars':
    submitter.memory = 4
  else:
    submitter.memory = 3

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/background_trials.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands, labels = [], []
## For this we'll need to do background trials in batches of 1000. So...
batchsize=100
total_bckg_trials = 30000

batches = int(np.round(total_bckg_trials/batchsize))

for i in range(batches):
    commands.append('{0} {1} --batch {2} --batchsize {3} --cat {4} --llhweight {5} --injweight {5} --years {6} --fom {7} --sirin {8} --mhubertest {9} --mese {10}'.format (env_shell, job_script, i, batchsize, catalog, llhweight, str(n), fom, sirin, mhubertest, mese))
    labels.append ('handle_batch_'+str(i))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter.submit_npx4 (commands, labels)
