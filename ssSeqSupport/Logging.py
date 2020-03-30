# Import third party modules
import os.path

# # Import ssSeqSupport variables
from . import Homedir, Logfilename


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
source folder/fastq_f: {}
fastq_r: {}
analysis only: {}
input_readlength: {}
troubleshoot mode: {}
n jobs: {}
Q-score cutoff: {}
alignment filter: {}
output folder: {} 
    """.format(args["datetime"], args["folder"], args["fastq_r"],
               args["analysis_only"], args["read_length"],  args["troubleshoot"],
               args["jobs"], args["q_cutoff"], args["alignment_filter"],
               args["output"])

    # Open the log file and append the latest logging information
    with open(Logfilename, "a") as f:
        f.write(logstr)
            
# Write a function that logs summary information. This function will only be
# triggered if a run completes successfully
def LogSummary(args):
    pass

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
    
    # Open the logfile
    with open(Logfilename, "a") as f:
        
        # Write the error
        f.write("\n\nError Encountered: {}".format(e))
        
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
    
    # Open the logfile
    with open(Logfilename, "a") as f:
        
        # Write the warning
        f.write("\n\nWarning: {}".format(w))

        # Print the error
        print(w)

    
# Write a function to log the identified file matches in the target folder
def LogInputFiles(matched_files, unmatched_files):
    
    # Open the logfilename and write the filenames
    with open(Logfilename, "a") as f:
        
        # Write a header
        f.write("\nIdentified Seq File Pairs:")
        
        # Write all filenames that matches
        for f_file, r_file in matched_files.items():
            f.write("\n\t- Forward Reads: {} \n\t- Reverse Reads: {}".format(f_file, r_file))
            
        # Write a header for the unmatched files
        f.write("\nUnmatched Files in Folder:")
        
        # Write all filenames that did not match
        for filename in unmatched_files:
            f.write(filename)
            
# Write a generic log function for logging information
def LogInfo(m):
    
    # Write the message to the log file
    with open(Logfilename, "a") as f:
        f.write("\n\n" + m)