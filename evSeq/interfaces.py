"""
Contains default arguments for evSeq. These are stored here so that we don't
need to change both cmd.py and gui.py each time we change arguments.

This file also contains a constructor function for the argparser components
as well as execution code that is shared between the GUI and CLI parsers.
"""
# Import 3rd party packages
import os
import argparse
import tqdm

from gooey import GooeyParser
from time import strftime

# Import evseq functions
from .util.globals import N_CPUS
from .util.input_processing import build_output_dirs
from .util.logging import log_init, log_info, log_error
from evSeq.run_evSeq import run_evSeq
    
# Get the working directory
CWD = os.getcwd()
    
# Define default args
AVERAGE_Q_CUTOFF = 25
BP_Q_CUTOFF = 30
LENGTH_CUTOFF = 0.9

MATCH_SCORE = 1
MISMATCH_PENALTY = 0
GAP_OPEN_PENALTY = 3
GAP_EXTENSION_PENALTY = 1

VARIABLE_THRESH = 0.2
VARIABLE_COUNT = 10

READ_LENGTH = None

# Build the shared parser constructor
def add_shared_args(parser):
    """
    Adds arguments to a GooeyParser or parser
    """    
    # Add an argument group for IO
    io_args_group = parser.add_argument_group("I/O", "Arguments specific to input/output of evSeq.")
    io_args_group.add_argument("--detailed_refseq",
                               help="Set if you are using different reference sequences by well",
                               required=False, 
                               action="store_true")
    io_args_group.add_argument("--analysis_only", 
                               help="Stop the run after analyzing read qualities",
                               required=False,
                               action="store_true")
    io_args_group.add_argument("--only_parse_fastqs",
                               help="Stop the run after generation of parsed, well-filtered fastq files",
                               required=False,
                               action="store_true")
    io_args_group.add_argument("--return_alignments", 
                               help="Save the constructed alignments to files",
                               required=False,
                               action="store_true")
    io_args_group.add_argument("--keep_parsed_fastqs",
                               help="Whether or not to keep the well-separated parsed filtered fastq files.",
                               required=False,
                               action="store_true")
    
    # Add read analysis parameters argument group
    params_args_group = parser.add_argument_group("Read Analysis",
                                                  "Parameters for filtering out poor quality reads.")
    params_args_group.add_argument("--average_q_cutoff", 
                                   help="Reads with an average Q-score below this are discarded.",
                                   required=False,
                                   default=AVERAGE_Q_CUTOFF,
                                   type=int)
    params_args_group.add_argument("--bp_q_cutoff", 
                                   help="Bases with a Q-score below this are not counted.",
                                   required=False,
                                   default=BP_Q_CUTOFF,
                                   type=int)
    params_args_group.add_argument("--length_cutoff", 
                                   help="Reads shorter than this fraction of `read_length` are discarded.",
                                   required=False,
                                   default=LENGTH_CUTOFF,
                                   type=float)
    params_args_group.add_argument("--match_score",
                                   help = "The bonus given to matching bases in an alignment score. Default = 1.",
                                   required = False,
                                   default = MATCH_SCORE,
                                   type = int)
    params_args_group.add_argument("--mismatch_penalty",
                                   help = "The penalty given to mismatching bases in an alignment score. Default = 0.",
                                   required = False,
                                   default = MISMATCH_PENALTY,
                                   type = int)
    params_args_group.add_argument("--gap_open_penalty",
                                   help = "The penalty given to opening a gap in an alignment score. Default = 3.",
                                   required = False,
                                   default = GAP_OPEN_PENALTY,
                                   type = int)
    params_args_group.add_argument("--gap_extension_penalty",
                                   help = "The penalty given to extending a gap in an alignment score. Default = 1.",
                                   required = False,
                                   default = GAP_EXTENSION_PENALTY,
                                   type = int)
    
    # Add a group containing arguments used for identifying variable positions
    variable_group = parser.add_argument_group("Position Identification",
                                                  "Parameters used for identifying variable positions.")
    variable_group.add_argument("--variable_thresh", 
                                help="Non-reference count frequency above which a position is considered to contain a variant.",
                                required=False,
                                default=VARIABLE_THRESH,
                                type=float)
    variable_group.add_argument("--variable_count", 
                                help="Wells with less reads than this are considered to contain nothing.",
                                required=False,
                                default=VARIABLE_COUNT,
                                type=int)
    
    # Add an advanced arguments group
    advanced_args_group = parser.add_argument_group("Advanced",
                                                    "Advanced use-case parameters")
    advanced_args_group.add_argument("--jobs", 
                                     help=f"Computer processors used. Must be between 1 and {N_CPUS}.",
                                     required=False,
                                     dest="jobs",
                                     default=N_CPUS - 1, 
                                     type=int)
    advanced_args_group.add_argument("--read_length", 
                                     help="Read length found in fastq files. Calculated if not specified.",
                                     required=False,
                                     type=int,
                                     default=READ_LENGTH)
    advanced_args_group.add_argument("--fancy_progress_bar", 
                                     help="EXPERIMENTAL, but usually works. Uses tqdm.gui to create a pop-out progress bar that shows instantaneous and average rates of processing.",
                                     required=False,
                                     action="store_true")
    
    return io_args_group

# Build the GUI parser
def build_gui_parser():
    
    # Instantiate GooeyParser
    parser = GooeyParser(description = "User Interface for evSeq")
        
    # Add required arguments
    required_args_group = parser.add_argument_group("Required Arguments", 
                                                    "Arguments required for each run")
    required_args_group.add_argument("refseq", 
                                     help="csv containing reference sequence information.",
                                     widget="FileChooser")
    required_args_group.add_argument("folder",
                                     help="Folder containing fastq or fastq.gz files. Can also be forward fastq or fastq.gz file.",
                                     widget="DirChooser")
    
    # Add arguments shared with CLI
    io_args_group = add_shared_args(parser)
    
    # Add additional IO args
    io_args_group.add_argument("--output", 
                               help="Save location for run. Defaults to same location as 'folder'.",
                               required=False,
                               default=None, # updated later
                               widget="DirChooser")
    io_args_group.add_argument("--fastq_r",
                               help="Reverse fastq or fastq.gz file. Usually not needed.",
                               required=False,
                               default="",
                               widget="FileChooser")
    
    return parser
    
# Build the CLI parser
def build_cli_parser():
    
    # Instantiate argparser
    parser = argparse.ArgumentParser()
    
    # Add required arguments
    required_args_group = parser.add_argument_group("Required Arguments",
                                                    "Arguments required for each run")
    required_args_group.add_argument("refseq",
                                     help="csv containing reference sequence information.")
    required_args_group.add_argument("folder",
                                     help="Folder containing fastq or fastq.gz files. Can also be forward fastq or fastq.gz file.")
    
    # Add arguments shared with CLI
    io_args_group = add_shared_args(parser)
    
    # Add additional IO args
    io_args_group.add_argument("--output",
                               help="Save location for run. Defaults to current working directory.",
                               required=False,
                               default=CWD)
    io_args_group.add_argument("--fastq_r",
                               help="Reverse fastq or fastq.gz file. Usually not needed.",
                               required=False,
                               default="")
    
    return parser

# Define shared execution code
def execute_evseq(gui = False):
    
    # Build the parser
    if gui:
        parser = build_gui_parser()
    else:
        parser = build_cli_parser()
        
    # Parse the arguments
    CL_ARGS = vars(parser.parse_args())
    
    # If GUI, determine output folder default
    if gui:
        if CL_ARGS['fastq_r'] == '':
            default_output = CL_ARGS['folder']
        else:
            default_output = os.path.dirname(CL_ARGS['folder'])
            
        if CL_ARGS['output'] is None:
            CL_ARGS['output'] = default_output
            
    # Identify the start time and add to the "CLArgs" dict. Also create
    # an output directory from the two and add this to CLArgs as well.
    base_output = CL_ARGS["output"]
    datetime = strftime("%Y%m%d-%H%M%S")
    output_dir = os.path.join(base_output, "evSeqOutput", datetime)
    CL_ARGS.update({"datetime": datetime, "output": output_dir})

    # Build all output directories
    build_output_dirs(CL_ARGS)

    # Log CL_ARGS
    log_init(CL_ARGS)
    
    # Set up progres bar
    if CL_ARGS['fancy_progress_bar']:
        tqdm_fn = tqdm.gui.tqdm
        
    elif gui:
        def blank(iterable, desc=None, total=None): return iterable
        tqdm_fn = blank # do nothing
        
    else:
        tqdm_fn = tqdm.tqdm
        
    # Run evSeq
    try:
        run_evSeq(CL_ARGS, tqdm_fn)
    except Exception as e:
       log_error(f"\nUnhandled exception encountered: '{e}'")

    # Log that we have successfully completed the run
    log_info("Run completed. Log may contain warnings.")
   