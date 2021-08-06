#!/usr/bin/env pythonw
# Import necessary modules
import subprocess

# Activate from command line
subprocess.run("conda activate evSeq", shell=True)
subprocess.run("evSeq-GUI", shell=True)