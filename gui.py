#!/usr/bin/env pythonw

# Import relevant modules
import os
from time import strftime

# Import relevant functions
import evSeq
from evSeq.util.logging import log_init, log_info, log_error
from evSeq.util.globals import N_CPUS
from evSeq.util.input_processing import build_output_dirs
from evSeq.run_evSeq import run_evSeq

# Import gooey
from gooey import Gooey, GooeyParser

# Get the path to evSeq repo
evSeq_path = os.path.dirname(os.path.dirname(os.path.abspath(evSeq.__file__)))

# Create a "main" function
@Gooey(program_name = "evSeq",
       required_cols = 1,
       optional_cols = 1,
       image_dir=os.path.join(evSeq_path, 'icons'))
def main():

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
    
    
    # Add an argument group for passing in the reverse file
    io_args_group = parser.add_argument_group("I/O",
                                                    "Arguments specific to input/output of evSeq.")
    io_args_group.add_argument("--output", 
                                     help="Save location for run. Defaults to same location as 'folder'.",
                                     required=False,
                                     default=None, # updated below
                                     widget="DirChooser")
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
    io_args_group.add_argument("--fastq_r",
                               help="Reverse fastq or fastq.gz file. Usually not needed.",
                               required=False,
                               default="",
                               widget="FileChooser")
    
    
    # Add read analysis parameters argument group
    params_args_group = parser.add_argument_group("Read Analysis",
                                                  "Parameters for filtering out poor quality reads.")
    params_args_group.add_argument("--average_q_cutoff", 
                                   help="Reads with an average Q-score below this are discarded.",
                                   required=False,
                                   default=25,
                                   type=int)
    params_args_group.add_argument("--bp_q_cutoff", 
                                   help="Bases with a Q-score below this are not counted.",
                                   required=False,
                                   default=30,
                                   type=int)
    params_args_group.add_argument("--length_cutoff", 
                                   help="Reads shorter than this fraction of `read_length` are discarded.",
                                   required=False,
                                   default=0.9,
                                   type=float)
    
    # Add a group containing arguments used for identifying variable positions
    variable_group = parser.add_argument_group("Position Identification",
                                                  "Parameters used for identifying variable positions.")
    variable_group.add_argument("--variable_thresh", 
                                   help="Non-reference count frequency above which a position is considered to contain a variant.",
                                   required=False,
                                   default=0.2,
                                   type=float)
    variable_group.add_argument("--variable_count", 
                                   help="Wells with less reads than this are considered to contain nothing.",
                                   required=False,
                                   default=10,
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
                                     default=None)


    # Parse the arguments
    CL_ARGS = vars(parser.parse_args())

    # Determine output folder default
    if CL_ARGS['fastq_r'] == '':
        default_output = CL_ARGS['folder']
    else:
        default_output = os.path.dirname(CL_ARGS['folder'])
    if CL_ARGS['output'] is None:
        CL_ARGS['output'] = default_output
    # Identify the cwd and start time and add to the "CLArgs" dict. Also create
    # an output directory from the two and add this to CLArgs as well.
    base_output = CL_ARGS["output"]
    datetime = strftime("%Y%m%d-%H%M%S")
    output_dir = os.path.join(base_output, "evSeqOutput", datetime)
    CL_ARGS.update({"datetime": datetime, "output": output_dir})

    # Build all output directories
    build_output_dirs(CL_ARGS)
    
    # Log CLArgs
    log_init(CL_ARGS)    

    # Run evSeq
    try:
        run_evSeq(CL_ARGS)
    except Exception as e:
        log_error(f"\nUnhandled exception encountered: '{e}'")
    
    # Log that we have successfully completed the run
    log_info("\nRun completed. Log may contain warnings.")

# Run main()
if __name__ == "__main__":
    main()
