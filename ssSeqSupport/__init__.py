"""Imports relevant modules and variables into the global scope"""
# Import third party modules relied upon by custom modules
import numpy as np
import pandas as pd
import re
import csv
import gzip
import os
import os.path
from Bio import pairwise2
from tqdm import tqdm
from multiprocessing import Pool
from collections import Counter
from time import strftime
from glob import glob

# Import plotting tools
import holoviews as hv
import colorcet as cc
import bokeh.io
hv.extension('bokeh')

# Ignore divide by 0.
np.seterr(invalid="ignore")

# Import system info now
from ssSeqSupport.SystemInfo import NCpus, Homedir, Logfilename

# Load all predefined global variables
from ssSeqSupport.PredefinedGlobals import (ReverseCompDict, AllowedBases, 
                                            AllowedBasesNoDeg, AllowedWells,
                                            CodonTable, BpToInd, IndToBp,
                                            AaToInd, IndToAa, AdapterLengthF, 
                                            AdapterLengthR, BarcodeLength,
                                            IdParser, AaOpts, BpOpts)

# Load loggers
from ssSeqSupport.Logging import (LogError, LogInit, LogInputFiles,
                                  LogWarning, LogInfo)

# Load support functions
from ssSeqSupport.SupportFuncs import (GetBlockInfo, CreateID, FindNNN,
                                       ReverseComplement, Translate,
                                       BuildOutputDirs)
                                       
# Import reference sequence
from ssSeqSupport.RefSeq import RefSeq, BuildRefSeqs

# Load viz functions
from ssSeqSupport.VizFuncs import GenerateSequencingHeatmap, GenerateReadQualChart

# Now load DNAObjects
from ssSeqSupport.DNAObjects import Translation, SeqPair

# Now load checker functions
from ssSeqSupport.InputChecks import CheckArgs, CheckRefSeqs, CheckIndexMap

# Now run/load relevant variables from GlobalSetup. This will construct global
# variables from provided input files
from ssSeqSupport.LoadFiles import LoadRefSeq, LoadDualInds, ConstructBCsToRefSeq

# Load plate and well objects followed by multiprocessing support
from ssSeqSupport.ssSeqMultiprocessing import MultiprocessPlateAnalyzer, MultiprocessPlateAnalyzerTS
from ssSeqSupport.PlateObjects import Well, Plate

# Now load the single required function from RunssSeq
from ssSeqSupport.RunssSeq import RunssSeq