# # Import third party modules/objects
from collections import Counter
from glob import glob
import numpy as np
import os.path

# # Load ssSeqSupport modules
from . import LogError, LogWarning
from . import AllowedBasesNoDeg, AllowedBases, AllowedWells
from . import NCpus

# Write a function that checks the validity of the index file
def CheckIndexMap(index_df):
    
    # Define the expected columns
    expected_cols = ("IndexPlate", "Well", "F-BC", "R-BC")
    
    # Pull the columns and make sure all expected ones exist
    index_cols = set(index_df.columns)
    missing_cols = [col for col in expected_cols if col not in index_cols]
    
    # If we are missing any columns, log an error and terminate
    if len(missing_cols) != 0:
        LogError("Expected columns missing from IndexMap.csv: {}".format(missing_cols))
        
    # Make sure all combinations of F-BC and R-BC are unique
    bc_combos = [(fbc, rbc) for fbc, rbc in 
                 index_df.loc[:, ["F-BC", "R-BC"]].itertuples(index=False)]
    bc_combo_counts = Counter(bc_combos)
    non_unique_combos = [combo for combo, count in bc_combo_counts.items()
                         if count > 1]
    
    # Report non-unique combos by throwing an error
    if len(non_unique_combos) > 0:
        LogError("""F-BC and R-BC combos must be unique in IndexMap.csv.
                  The following F-BC, R-BC combos are not unique: {}""".format(non_unique_combos))
        
    # Loop over each row in the input file
    for i, (_, row) in enumerate(index_df.iterrows()):
        
        # Check that each value is present. If not, throw an error.
        if any([row[col]==None for col in expected_cols]):
            LogError("Empty value in row {} of IndexMap.csv.".format(i))
            
        # Confirm that the forward barcode sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char not in AllowedBasesNoDeg for char in row["F-BC"]]):
            LogError("F-BC sequence in row {} has base other than 'A', 'C', 'T', or 'G'")
        
        # Confirm that the reverse barcode sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char not in AllowedBasesNoDeg for char in row["R-BC"]]):
            LogError("R-BC sequence in row {} has base other than 'A', 'C', 'T', or 'G'")

        # Confirm that the named well is a real well
        if row["Well"] not in AllowedWells:
            LogError("""Unexpected well in row {}: {} of reference sequence file. 
                        Well must take form 'A##'""".format(i, row["Well"]))
            
    # Confirm that all barcodes are the same length
    pairwise_check = np.equal(index_df["F-BC"].str.len().values,
                              index_df["R-BC"].str.len().values)
    if not all(pairwise_check):
        LogError("Barcodes must all be the same length. Check IndexMap.csv")

# Write a function that checks the validity of the reference sequence file
def CheckRefSeqs(ref_seqs_df, detailed_file):
    
    # Define expected_cols based on whether we are expecting a detailed reference
    # sequence file or not
    if detailed_file:
        expected_cols = ("PlateName", "IndexPlate", "Well", "ReferenceSequence")
    else:
        expected_cols = ("PlateName", "IndexPlate", "ReferenceSequence")
        
    # Identify columns missing from the reference sequence file
    ref_seq_cols = set(ref_seqs_df.columns)
    missing_cols = [col for col in expected_cols if col not in ref_seq_cols]
    
    # If we are missing any columns, log an error and terminate the program
    if len(missing_cols) != 0:
        LogError("Expected columns missing from refseq file: {}".format(missing_cols))
                
    # Get unique plate names. Create a dictionary for confirming that a unique
    # plate nickname always goes with a unique index plate
    index_to_nick_check = {}
    
    # Loop over each row in the input file and check the values
    for i, (_, row) in enumerate(ref_seqs_df.iterrows()):
        
        # Check that each value is present. If not, throw an error.
        if any([row[col]==None for col in expected_cols]):
            LogError("Empty value in row {} of refseqs file.".format(i))
            
        # Confirm that the reference sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char not in AllowedBases for char in row["ReferenceSequence"]]):
            LogError("Reference sequence in row {} has base other than 'A', 'C', 'T', 'G', or 'N'")
                        
        # Make sure that the same plate nickname and dual index plate name always 
        # go together.
        if row["PlateName"] not in index_to_nick_check:
            index_to_nick_check[row["PlateName"]] = row["IndexPlate"]
        else:
            if index_to_nick_check["PlateName"] != row["IndexPlate"]:
                LogError("""Each plate nickname must be associated with a single dual index plate.
                         In row {}, nickname {} is duplicated with index plate {}""".format(i, row["PlateName"],
                                                                                            row["IndexPlate"]))
        
        # Additional row checks for a detailed file:
        if detailed_file:
            
            # Confirm that the named well is a real well
            if row["Well"] not in AllowedWells:
                LogError("""Unexpected well in row {}: {} of reference sequence file. 
                          Well must take form 'A##'""".format(i, row["Well"]))
    
    # Additional wholistic checks for a detailed file: Make sure that there are 
    # no duplicate dual index plate-well combos
    if detailed_file:
        
        # Get all combos of dual index plate and well. 
        combos = [(row["IndexPlate"], row["IndexWell"]) 
                  for _, row in ref_seqs_df.iterrows()]
        
        # Count the number of occurences of each. Identify those combos with 
        # more than one occurence
        counts = Counter(combos)
        duplicates = [combo for combo, count in counts.items() if count > 1]
        
        # Throw an error if there are duplicates
        if len(duplicates) > 0:
            LogError("""Each sample must have a unique index plate and well.
                     The following combinations of plate and well occured more than once in the refseq file:
                     {}""".format(set(duplicates)))
           
# Write a function that checks the arguments passed into the command line prompt
def CheckArgs(args):
    
    # Confirm that the csv file with reference sequences exists
    if not os.path.exists(args["refseq"]):
        
        # Write the error and terminate the program
        LogError("The file or folder does not exist: {}".format(args["refseq"]))
        
    # Confirm that the reference sequence csv file is indeed a file
    if not os.path.isfile(args["refseq"]):
        
        # Write the error and terminate the pgoram
        LogError("Specified fastq-r file is not a file: {}".format(args["refseq"]))
    
    # Confirm that the specified folder containing sequencing results exists 
    if not os.path.exists(args["folder"]):
        
        # Write the error and terminate the program
        LogError("The file or folder does not exist: {}".format(args["folder"]))
    
    # Confirm that we have a folder or a file in the positional argument
    if os.path.isdir(args["folder"]):
        
        # Record that we are working with a folder input
        folder = True
        
        # If we have a folder, make sure there are fastq or fastq.gz files inside of it
        all_files = glob(os.path.join(args["folder"], "*.fastq*"))
        
        # If we can't find any files, report the error
        if len(all_files) == 0:
            
            # Write the error and terminate the program
            LogError("No fastq or fastq.gz files found in target folder.")
                
    # If we have a file, make sure the fastq_r flag is also thrown and another
    # file is given
    elif os.path.isfile(args["folder"]):
        
        # Record that we are working with a file input
        folder = False
        
        # Check to be sure that the fastq_r flag is thrown
        if args["fastq_r"] == "":
            
            # Write the error and terminate the program
            LogError("--fastq-r must also be specified.")
            
        # Check to be sure that the fastq_r flag exists
        if not os.path.exists(args["fastq_r"]):
            
            # Write the error and terminate the program
            LogError("Specified fastq-r file does not exist: {}".format(args["fastq_r"]))
            
        # Check to be sure that the fastq_r flag contains a file
        if not os.path.isfile(args["fastq_r"]):
            
            # Write the error and terminate the pgoram
            LogError("Specified fastq-r file is not a file: {}".format(args["fastq_r"]))
            
        # Check to be sure that the fastq_r flag contains a file and that it is
        # either a fastq or fastq.gz    
        if "fastq" not in args["fastq_r"] and "fastq.gz" not in args["fastq_r"]:
            
            # Write the error and terminate the program
            LogError("The reverse file is neither a fastq nor fastq.gz file.")
            
        # Check to be sure that the fastq_f flag contains a file and that it is
        # either a fastq or fastq.gz    
        if "fastq" not in args["folder"] and "fastq.gz" not in args["folder"]:
            
            # Write the error and terminate the program
            LogError("The forward file is neither a fastq nor fastq.gz file.")
            
        # Construct all_files if we have made it this far
        all_files = [args["folder"], args["fastq_r"]]
        
    # If we haven't found a file or a folder, throw an error
    else:
        
        # Log the error and terminate the program
        LogError("The positional argument is neither a file nor a folder.")
        
    # Confirm that the read length >0. Skip this if it is 'None' as this indicates
    # that we will be calculating it later on
    if args["read_length"] is not None:
        if not args["read_length"] > 0:
            
            # Write the error and terminate the program
            LogError("--read_length must be an integer greater than 0.")
    
    # Confirm that the Q-score cutoff >0
    if not args["q_cutoff"] > 0:
        
        # Write the error and terminate the program
        LogError("--q_cutoff must be an integer greater than 0.")
    
    # Confirm that the alignment filter is a float between 0 and 1
    if not (args["alignment_filter"] >=0 and args["alignment_filter"] <= 1):
        
        # Write the error and terminate the program
        LogError("--alignment_filter must be a float between 0 and 1.")
    
    # Confirm that the number of jobs is an integer between 1 and the number of 
    # cpu's available to the computer. Adjust the number to 1 if below 1 and to
    # the maximum number of cpus if above the maximum available.
    if not args["jobs"] >= 1:
        
        # Raise the number of jobs to 1 and record warning
        args["jobs"] = 1
        LogWarning("--jobs must be greater than or equal to 1. Defaulting to 1 job.")
    
    elif not args["jobs"] <= NCpus:
        
        # Lower the number of jobs to the maximum available and record warning
        args["jobs"] = NCpus
        LogWarning("--jobs must be less than or equal to the available number of cpus. Defaulting to max number of cpus.")
        
    # Return all_files
    return all_files, folder