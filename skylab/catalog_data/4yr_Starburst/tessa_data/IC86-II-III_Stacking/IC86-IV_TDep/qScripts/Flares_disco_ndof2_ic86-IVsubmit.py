#!/usr/bin/env python

import os
import sys
import glob
import socket

def call(args, cwd=None, stdout=None, stderr=None, env=None):
    """
    Wrap the call function of the subprocess module for older pythons.
    """

    if not cwd:
        cwd = os.getcwd()
    if not stdout:
        stdout = sys.stdout
    if not stderr:
        stderr = sys.stderr
    if not env:
        env = os.environ

    try:
        import subprocess
        return subprocess.call(args,cwd=cwd,stdout=stdout,stderr=stderr,env=env)
    except:
        import popen2
        if cwd:
            os.chdir(cwd)
        cmd = " ".join(args)
        proc = popen2.Popen4(cmd)
        for line in proc.fromchild:
            stdout.write(line)
        rc = proc.wait()
        return rc


if __name__ == "__main__":
    """
    Main function for the script.
    """

    host = socket.gethostname()
    print host

    print "Submitting job with input and output files:"
    if host=="lightcube":
        env="/home/christov/psLab/start.sh"
        files = glob.glob("/home/christov/Fermi_LAT/54971_to_57160/lcFinalFMonitored_from_54971_to_57160/*.rawBlock.dat")
    elif host=="submitter.icecube.wisc.edu":
        env = "/data/user/achristov/psLab/start.sh"
        files = glob.glob("/data/user/achristov/psLab/macro_llh/IC86-IV_TDep/lcAll/*.dat")

    gendir="flares/discoGenScripts/ic86-I-II-III-IV-79-59/"
    rezdir="flares/discoRez/ic86-I-II-III-IV-79-59/"

    if not os.path.exists(gendir):
        os.makedirs(gendir)
    if not os.path.exists(rezdir):
        os.makedirs(rezdir)
    
    for blockfilepath in files:
        
        blockfile = blockfilepath.split("/")[-1]
        
        ra = float(blockfile.split("_")[-3])
        
        declination = float(blockfile.split("_")[-2])
        
        source = blockfile.split("_")[1]+"_"+blockfile.split("_")[2]

        fmaster = open("Flares_disco_ndof2_ic86-IV_MASTER.C")
        
        outfile = gendir+"Flares_disco_ndof2_ic86-IV_"+source+".C"

        out = file(outfile,"w")
        for line in fmaster:
            line = line.replace("MASTER_BLOCKFILE_MASTER", blockfilepath)
            line = line.replace("MASTER_RA_MASTER", str(ra))
            line = line.replace("MASTER_DECLINATION_MASTER", str(declination))
            line = line.replace("MASTER_OUT_FILE", rezdir+blockfile.rstrip(".rawBlock.dat")+"_ndof2_disco")
            if host=="lightcube":
                line = line.replace("MASTER_SITEATT","_lc")
            elif host=="submitter.icecube.wisc.edu":
                line = line.replace("MASTER_SITEATT","_co")
            out.write(line)
            
        if host == "submitter.icecube.wisc.edu":
            command = "\'/data/user/achristov/psLab/macro_llh/IC86-IV_TDep/qScripts/"+outfile+"'"
            print command
            call(["/home/achristov/bin/submit_npx4_long.sh", env, "root", "-l", "-b", "-q", command])
            
        else:
            print "This is not the cluster!"
            command = "/home/christov/psLab/macro_llh/IC86-IV_TDep/qScripts/"+outfile
            print command
            #call([ "root", "-l", "-b", "-q", command])
        
