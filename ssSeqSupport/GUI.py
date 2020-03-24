#!/usr/bin/env python

# Import relevant modules
import argparse
import os
from time import strftime

# Import relevant modules from ssSeqSupport.
from ssSeqSupport import LogInit, LogInfo, LogError
from ssSeqSupport import RunssSeq
from ssSeqSupport import NCpus

# Import gooey
from gooey import Gooey, GooeyParser

# Create a "main" function
@Gooey(program_name = "ssSeq",
       required_cols = 1,
       optional_cols = 1)
def main():

    # Get the cwd 
    cwd = os.getcwd()

    # Instantiate GooeyParser
    parser = GooeyParser(description = "User Interface for ssSeq")

    # Create the first argument group
    required_args_group = parser.add_argument_group("Required Arguments", 
                                                    "Arguments required for each run")
    required_args_group.add_argument("refseq", 
                                     help = "csv containing reference sequences.",
                                     widget = "FileChooser")
    required_args_group.add_argument("folder",
                                     help = "Folder containing fastq or fastq.gz files",
                                     widget = "DirChooser")
    
    # Add an argument group for passing in the reverse file
    optional_args_group = parser.add_argument_group("Optional Arguments",
                                                    "Optional arguments for specific use cases.")
    optional_args_group.add_argument("--fastq_r", 
                                     help = "Reverse fastq or fastq.gz file. Usually not needed.",
                                     required = False,
                                     default = "",
                                     widget = "FileChooser")
    optional_args_group.add_argument("--output", 
                                     help = "Save location for run.",
                                     required = False,
                                     default = cwd,
                                     widget = "DirChooser")
    
    # Add read analysis parameters argument group
    params_args_group = parser.add_argument_group("Read Analysis Parameters",
                                                  "Parameters for filtering out poor quality reads.")
    params_args_group.add_argument("--q_cutoff", 
                                   help = "Reads below this Q-score are discarded.",
                                   required = False,
                                   default = 30,
                                   type = int)
    params_args_group.add_argument("--alignment_filter", 
                                   help = "Reads with normalized alignment score below this value are discarded.",
                                   required = False,
                                   default = 0.5,
                                   type = float)
    
    # Add an advanced argument group
    advanced_args_group = parser.add_argument_group("Advanced",
                                                    "Advanced use-case parameters")
    advanced_args_group.add_argument("--analysis_only", 
                                     help = "Check box to perform quality analysis only",
                                     required = False, 
                                     action = "store_true")
    advanced_args_group.add_argument("--detailed_refseq",
                                     help = "Check box if you are using different reference sequences by well",
                                     required = False, 
                                     action = "store_true")
    advanced_args_group.add_argument("--troubleshoot", 
                                     help = "Check box to run in troubleshoot mode",
                                     action = "store_true", 
                                     required = False)
    advanced_args_group.add_argument("--jobs", 
                                     help = "Computer processors used. Must be between 1 and {}.".format(NCpus),
                                     required = False,
                                     dest = "jobs",
                                     default = NCpus - 1, 
                                     type = int)
    advanced_args_group.add_argument("--read_length", 
                                     help = "Read length found in fastq files. Calculated if not specified.",
                                     required = False,
                                     type = int,
                                     default = None)
    

    # Parse the arguments
    CLArgs = vars(parser.parse_args())
        
    # Identify the cwd and start time and add to the "CLArgs" dict. Also create an
    # output directory from the two and add this to CLArgs as well.
    base_output = CLArgs["output"]
    datetime = strftime("%Y%m%d-%H%M%S")
    output_dir = os.path.join(base_output, "ssSeq_Output", datetime)
    CLArgs.update({"datetime": datetime, "output": output_dir})
    
    # Log CLArgs
    LogInit(CLArgs)    

    # Run ssSeq
    try:
        RunssSeq(CLArgs)
    except Exception as e:
        LogError("\nUnhandled exception encountered: '{}'".format(e))
    
    # Log that we have successfully completed the run
    LogInfo("\nRun completed. Log may contain warnings.")
    
# Run main()
if __name__ == "__main__":
    main()