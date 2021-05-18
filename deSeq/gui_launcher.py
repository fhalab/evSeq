#!/usr/bin/env python
# Import necessary modules
import subprocess
import os

# # Get the location of ssSeqGui (same directory as this file)
gui_loc = os.path.dirname(os.path.realpath(__file__))

# # Change cwd
os.chdir(gui_loc)

# Activate from command line
subprocess.run("conda activate deSeq_test", shell = True)
subprocess.run("deSeq-GUI", shell=True)