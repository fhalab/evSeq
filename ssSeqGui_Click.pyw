#!/usr/bin/env pythonw

# Import subprocess
import subprocess
import os

# Change cwd
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Run the sequence
subprocess.run('conda run -n ssSeq pythonw ssSeqGui', shell=True)