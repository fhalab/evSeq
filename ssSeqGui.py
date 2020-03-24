#!/usr/bin/env python

# Import subprocess
import subprocess
import os

# Change cwd
os.chdir(os.path.dirname(__file__))

# Run the sequence
subprocess.run('conda run -n ssSeq python ssSeqSupport/GUI.py', shell=True)