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

    print "Submitting job with input and output files:"
    
    env = "/atlas/users/christov/psLab/start.sh"

    files = glob.glob("/atlas/users/christov/psLab/macro_llh/IC86-IV_TDep/lcAll/*.dat")
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
            line = line.replace("MASTER_OUT_FILE", "/atlas/users/christov/psLab/macro_llh/IC86-IV_TDep/qScripts/"+rezdir+blockfile.rstrip(".rawBlock.dat")+"_ndof2_disco")
            if host=="atlas074":
                line = line.replace("MASTER_SITEATT","_at")

            out.write(line)
                        
        if host == "atlas074":
            command = "/atlas/users/christov/psLab/macro_llh/IC86-IV_TDep/qScripts/"+outfile
            print command
            if not os.path.exists(gendir+source):
                os.makedirs(gendir+source)

            script_file=open(gendir+source+"/submit_script.sh", "w+")
            script_file.write("#/bin/bash\n")
            script_file.write("export DISPLAY=\"\"\n")
            script_file.write("export OS_ARCH=RHEL_6_x86_64\n")
            script_file.write("tmpDir=`mktemp -d`\n")
            script_file.write("cd ${tmpDir}\n")
            script_file.write("echo 'we are in '${PWD}\n")
            script_file.write("mkdir rundir\n")
            script_file.write("ln -s /atlas/users/christov/data/IC59/IC59_GC/\n")
            script_file.write("cd rundir\n")
            script_file.write("ln -s /atlas/users/christov/data/IC79/BDT_GRLv2_5/\n")
            script_file.write("ln -s /atlas/users/christov/data/IC79/BDT_v4/\n")
            script_file.write("ln -s /atlas/users/christov/data/IC59/data/\n")
            script_file.write("eval `/cvmfs/icecube.opensciencegrid.org/py2-v1/setup.sh`\n")
            script_file.write(env+" \'root -qlb "+command+"\'\n")
            script_file.close()

            call(["qsub -l mem=5000mb -l vmem=5000mb -q verylong -j oe -o "+gendir+source+"/output.txt", gendir+source+"/submit_script.sh"])
            #call(["qsub -l mem=5000mb -l vmem=5000mb -q veryshort -j oe -o "+gendir+source+"/output.txt", gendir+source+"/submit_script.sh"])

        else:
            print "This is not the cluster!"
            command = "/atlas/users/christov/psLab/macro_llh/IC86-IV_TDep/qScripts/"+outfile
            print command
            #call([ "root", "-l", "-b", "-q", command])
        
