#!/usr/bin/env python
# Import necessary modules
import subprocess
import os

# Get the location of ssSeqGui (same directory as this file)
_ssSeqGuiLoc = os.path.dirname(os.path.realpath(__file__))

# Change cwd
os.chdir(_ssSeqGuiLoc)

# Activate from command line
subprocess.run("conda run -n ssSeq python ./ssSeqGui", shell = True)