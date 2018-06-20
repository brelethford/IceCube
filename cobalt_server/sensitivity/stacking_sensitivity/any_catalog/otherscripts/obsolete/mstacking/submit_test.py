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
parser.add_option ('--llhweight', dest = 'llhweight', type = str,
                default = None, metavar = 'LLHWEIGHT',
                help = 'Sets the weighting used in the llh model.')

parser.add_option ('--cat', dest = 'catalog', type = str,
                metavar = 'CATALOG',
                help = 'Selects the catalog of sources.')

parser.add_option ('--injweight', dest = 'injweight', type = str,
                default = None, metavar = 'INJWEIGHT',
                help = 'Sets the weighting used in the injection model.')

parser.add_option ('--n', dest = 'n',
                default = 4, metavar = 'N',
                help = 'Number of years of data')

parser.add_option ('--fom', dest = 'fom', type = float,
                default = 0.0, metavar = 'FOM',
                help = 'Sets the figure of merit cut for WHSP blazars.')

parser.add_option ('--sirin', dest = 'sirin',
                default = 0, action='store_true', metavar = 'SIRIN',
                help = 'Determines if unsplined IC79 is used')

parser.add_option ('--mhubertest', dest = 'mhubertest',
                default = 0, action='store_true', metavar = 'MHUBERTEST',
                help = 'toggles mhuber code priority')

parser.add_option ('--mese', dest = 'mese',
                default = 0, action='store_true', metavar = 'MESE',
                help = 'Toggles inclusion of mese data for 7yr+mese things.')

opts, args = parser.parse_args ()
catalog = opts.catalog
llhweight = opts.llhweight
injweight = opts.injweight
n = opts.n
fom = opts.fom
sirin = opts.sirin
mhubertest = opts.mhubertest
mese = opts.mese


if mese:
  print ("Using MESE data")
  mese=1
if mhubertest:
  print ("Using mhubertest")
  mhubertest=1
if sirin:
  print ("Using sirin's IC79 dataset")
  sirin = 1
else:
  print ("Using kai's IC79 dataset")


try:
   print ("Catalog used: " + catalog)
except:
   print ('No catalog selected! Please select a catalog.')

if llhweight is None:
        if catalog == 'SwiftBAT70m':
                llhweight = 'uniform'
        elif catalog == '4yr_Starburst':
                llhweight = 'flux'
        elif catalog == 'SNR_noPWN' or catalog == 'SNR_cloud' or catalog == 'SNR_PWN':
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
job_id = '{0}_nixtime_{2:.0f}_job_{1}_{3}inj'.format (
        socket.gethostname (), os.getpid (), time.time (), injweight)
# set up submitter
jobdir = misc.ensure_dir ('/data/user/brelethford/Output/stacking_sensitivity/{0}/mstacking_{1}yr/jobs/'.format (catalog, str(n)))
submitter_sens = Submitter (job_dir=jobdir+'sens_{0}_{1}'.format (job_id, llhweight))
submitter_disc = Submitter (job_dir=jobdir+'disc_{0}_{1}'.format (job_id, llhweight))

if catalog == 'WHSP_Blazars':
  submitter_sens.memory = 25
  submitter_disc.memory = 25
elif catalog == 'SwiftBAT70m':
  submitter_sens.memory = 28
  submitter_disc.memory = 28
elif catalog == '2LAC':
  submitter_sens.memory = 22
  submitter_disc.memory = 22
else:
  submitter_sens.memory = 13
  submitter_disc.memory = 13
  

# figure out what dir we're in, and get the path to the actual job script
this_dir = os.path.dirname (os.path.abspath (__file__))
job_script = this_dir + '/calculate_sensdisc.py'

# get the metaproject path
env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'

# build the lists of commands and labels
commands_sens, commands_disc, labels_sens, labels_disc = [], [], [], []
bins=1
gammarange = np.linspace(2.0,2.0,bins)

for gamma in gammarange:
    opts_list = '{0} {1} --cat {2} --llhweight {3} --injweight {4} --n {5} --gamma {6} --fom {7} --sirin {8} --mhuber {9} --mese {10}'.format (env_shell, job_script, catalog, llhweight, injweight, n, gamma, fom, sirin, mhubertest, mese)
    commands_sens.append(opts_list)
    labels_sens.append ('handle_sensitivity_gamma_{}'.format(gamma))
    commands_disc.append(opts_list + ' --disc 1')
    labels_disc.append ('handle_disc_gamma_{}'.format(gamma))

# submit! This will return several batches of background events. I'll have to synthesize these together afterwards.
submitter_sens.submit_npx4 (commands_sens, labels_sens)
submitter_disc.submit_npx4 (commands_disc, labels_disc)
