# Import third party modules
import os.path

# # Import functions
from .globals import UTILDIR, LOG_FILENAME

# Write a function to log all cl_args passed in from each run
def log_init(cl_args):
    """
    Logs the initialization variables passed in through either the command line
    or the GUI submission sheet.
    
    Parameters
    ----------
    cl_args: dict: Dictionary resulting from calling 'parser.parse_args()' in the 
        initial setup of the run.
        
    Returns
    -------
    None
    """
    
    # Format a string for adding to the logfile
    logstr = f"""

{cl_args["datetime"]}
----------------------------------------------------------------------------
Required Arguments:
    refseq: {cl_args["refseq"]}
    folder: {cl_args["folder"]}
    
Input/Output Arguments:
    fastq_r: {cl_args["fastq_r"]}
    output: {cl_args["output"]}
    detailed_refseq: {cl_args["detailed_refseq"]}
    analysis_only: {cl_args["analysis_only"]}
    stop_after_fastq: {cl_args["stop_after_fastq"]}
    return_alignments: {cl_args["return_alignments"]}
    
Read Analysis:
    average_q_cutoff: {cl_args["average_q_cutoff"]}
    bp_q_cutoff: {cl_args["bp_q_cutoff"]}
    length_cutoff: {cl_args["length_cutoff"]}

Position Identification:
    variable_thresh: {cl_args["variable_thresh"]}
    variable_count: {cl_args["variable_count"]}
    
Advanced:
    jobs: {cl_args["jobs"]}
    read_length: {cl_args["read_length"]}

Logged Messages:
    """

    # Build the run-specific log file and make it global
    # Export RunSpecLog as a global variable
    global RUN_SPEC_LOG
    RUN_SPEC_LOG = os.path.join(cl_args["output"], "RunSpecificLog.txt")

    # Add to logs
    write_to_log(logstr)

# Write a function that writes a message to both the RunSpecific and continuous logs
def write_to_log(message):
    
    # Add to the run-specific log
    with open(RUN_SPEC_LOG, "a") as f:
        f.write(message)

    # Open the continuous log file and append the latest logging information
    with open(LOG_FILENAME, "a") as f:
        f.write(message)

# Write a function that logs a warning encountered during the rurn
def log_warning(w):
    """
    Writes an warning message to the log file. 
    
    Parameters
    ----------
    w: Warning string containing the captured warning to log.
    
    Returns
    -------
    None
    """
    
    # Define the message
    m = "\n\nWarning: {}".format(w)
    
    # Write the message
    write_to_log(m)
    
    # Print the warning
    print(w)

# Write a function that logs any critical error encountered during the run
def log_error(e):
    """
    Writes an error message to the log file. 
    
    Parameters
    ----------
    e: Exception class object containing the captured error to log or a string 
        containing the error message.
    
    Returns
    -------
    None
    """
    
    # Define the message
    m = "\n\nError Encountered: {}".format(e)
    
    # Write the message
    write_to_log(m)
    
    # Print the error
    print(e)
    
    # Terminate the program
    quit()
     
# Write a function to log the identified file matches in the target folder
def log_input_file(forward_file, reverse_file, unmatched_files):
    
    # Define the message
    message = f"""
Identified Seq File Pairs:
    Forward Reads: {forward_file}
    Reverse Reads: {reverse_file}
    
Unmatched Files in Folder:
    """
    
    # Add the unmatched files to the message
    for filename in unmatched_files:
        message += "\n\t {}".format(filename)
    
    # Write the message
    write_to_log(message)
                
# Write a generic log function for logging information
def log_info(m):
    
    # Define the message
    message = "\n\n" + m
    
    # Write the message
    write_to_log(message)