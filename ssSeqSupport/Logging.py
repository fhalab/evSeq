# Import third party modules
import os.path

# # Import ssSeqSupport variables
from . import Homedir, Logfilename

# Write a function that writes a message to both the RunSpecific and continuous logs
def WriteToLog(message):
    
    # Add to the run-specific log
    with open(RunSpecLog, "a") as f:
        f.write(message)

    # Open the continuous log file and append the latest logging information
    with open(Logfilename, "a") as f:
        f.write(message)

# Write a function to log all args passed in from each run
def LogInit(args):
    """
    Logs the initialization variables passed in through either the command line
    or the GUI submission sheet.
    
    Parameters
    ----------
    args: dict: Dictionary resulting from calling 'parser.parse_args()' in the 
        initial setup of the run.
        
    Returns
    -------
    None
    """
    
    # Format a string for adding to the logfile
    logstr = """

{}
----------------------------------------------------------------------------
refseq: {}
source folder/fastq_f: {}
fastq_r: {}
analysis only: {}
input_readlength: {}
troubleshoot mode: {}
n jobs: {}
Q-score cutoff: {}
alignment filter: {}
output folder: {} 
    """.format(args["datetime"], args["refseq"], args["folder"], args["fastq_r"],
               args["analysis_only"], args["read_length"],  args["troubleshoot"],
               args["jobs"], args["q_cutoff"], args["alignment_filter"],
               args["output"])

    # Build the run-specific log file and make it global
    # Export RunSpecLog as a global variable
    global RunSpecLog
    RunSpecLog = os.path.join(args["output"], "RunSpecificLog.txt")

    # Add to logs
    WriteToLog(logstr)

# Write a function that logs any critical error encountered during the run
def LogError(e):
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
    WriteToLog(m)
    
    # Print the error
    print(e)
    
    # Terminate the program
    quit()
    
# Write a function that logs a warning encountered during the rurn
def LogWarning(w):
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
    WriteToLog(m)
    
    # Print the warning
    print(w)
    
# Write a function to log the identified file matches in the target folder
def LogInputFiles(matched_files, unmatched_files):
    
    # Define the message
    message = """
Identified Seq File Pairs:
    Forward Reads: {}
    Reverse Reads: {}
    
Unmatched Files in Folder:
    """
    
    # Add the unmatched files to the message
    for filename in unmatched_files:
        message += "\n\t {}".format(filename)
    
    # Write the message
    WriteToLog(message)
                
# Write a generic log function for logging information
def LogInfo(m):
    
    # Define the message
    message = "\n" + m
    
    # Write the message
    WriteToLog(message)