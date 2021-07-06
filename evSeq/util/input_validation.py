# # Import third party modules/objects
from collections import Counter
from glob import glob
import numpy as np
import os.path

# # Load functions
from .logging import log_error, log_warning
from .globals import ALLOWED_BASES_NO_DEG, ALLOWED_BASES, ALLOWED_WELLS, N_CPUS

# Write a function that checks the validity of the index file
def check_index_map(index_df):
    
    # Define the expected columns
    expected_cols = ("IndexPlate", "Well", "FBC", "RBC")
    
    # Pull the columns and make sure all expected ones exist
    index_cols = set(index_df.columns)
    missing_cols = [col for col in expected_cols if col not in index_cols]
    
    # If we are missing any columns, log an error and terminate
    if len(missing_cols) != 0:
        log_error(f"Expected columns missing from index_map.csv: {missing_cols}")
        
    # Make sure all combinations of F-BC and R-BC are unique
    bc_combos = [(fbc, rbc) for fbc, rbc in 
                 index_df.loc[:, ["FBC", "RBC"]].itertuples(index=False)]
    bc_combo_counts = Counter(bc_combos)
    non_unique_combos = [combo for combo, count in bc_combo_counts.items()
                         if count > 1]
    
    # Report non-unique combos by throwing an error
    if len(non_unique_combos) > 0:
        log_error(f"""FBC and RBC combos must be unique in index_map.csv.
                  The following FBC, RBC combos are not unique: {non_unique_combos}""")
        
    # Loop over each row in the input file
    for i, (_, row) in enumerate(index_df.iterrows()):
        
        # Check that each value is present. If not, throw an error.
        if any([row[col]==None for col in expected_cols]):
            log_error(f"Empty value in row {i} of index_map.csv.")
            
        # Confirm that the forward barcode sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char not in ALLOWED_BASES_NO_DEG for char in row["FBC"]]):
            log_error(f"F-BC sequence in row {i} has base other than 'A', 'C', 'T', or 'G'")
        
        # Confirm that the reverse barcode sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char not in ALLOWED_BASES_NO_DEG for char in row["RBC"]]):
            log_error(f"R-BC sequence in row {i} has base other than 'A', 'C', 'T', or 'G'")

        # Confirm that the named well is a real well
        if row["Well"] not in ALLOWED_WELLS:
            log_error(f"""Unexpected well in row {i}: {row["Well"]} of reference sequence file. 
                        Well must take form 'A##'""")
            
    # Confirm that all barcodes are the same length
    pairwise_check = np.equal(index_df["FBC"].str.len().values,
                              index_df["RBC"].str.len().values)
    if not all(pairwise_check):
        log_error("Barcodes must all be the same length. Check index_map.csv")

# Write a function that checks the validity of the reference sequence file
def check_ref_seqs(ref_seqs_df, detailed_file):
    
    # Define expected_cols based on whether we are expecting a detailed reference
    # sequence file or not
    if detailed_file:
        expected_cols = ("PlateName", "IndexPlate", "Well", "ReferenceSequence",
                         "InFrameBase", "BpIndStart", "AaIndStart")
    else:
        expected_cols = ("PlateName", "IndexPlate", "ReferenceSequence",
                         "InFrameBase", "BpIndStart", "AaIndStart")
        
    # Identify columns missing from the reference sequence file
    ref_seq_cols = set(ref_seqs_df.columns)
    missing_cols = [col for col in expected_cols if col not in ref_seq_cols]
    
    # If this is not a detailed file but we have a 'Well' column, throw an error
    if not detailed_file and "Well" in ref_seq_cols:
        log_error("It looks like you're trying to pass in a detailed reference sequence file. "
                 "You must throw the 'detailed_refseq' flag if this is your intent.")
    
    # If we are missing any columns, log an error and terminate the program
    if len(missing_cols) != 0:
        log_error(f"Expected columns missing from refseq file: {missing_cols}")
                
    # Get unique plate names. Create a dictionary for confirming that a unique
    # plate nickname always goes with a unique index plate
    index_to_nick_check = {}
    
    # Loop over each row in the input file and check the values
    for i, (_, row) in enumerate(ref_seqs_df.iterrows()):
        
        # The InFrameBase must be between 1 and 3, inclusive.
        if not (1 <= row["InFrameBase"] <= 3):
            log_error(f"InFrameBase must be between 1 and 3. Check row {i}.")
            
        # The BpIndStart and AaIndStart must be positive or zero
        if row["BpIndStart"] < 0:
            log_error(f"BpIndStart must be greater than or equal to 0. Check row {i}.")
        if row["AaIndStart"] < 0:
            log_error(f"AaIndStart must be greater than or equal to 0. Check row {i}.")
        
        # Check that each value is present. If not, throw an error.
        if any([row[col]==None for col in expected_cols]):
            log_error(f"Empty value in row {i} of refseqs file.")
            
        # Confirm that the reference sequence is made up entirely of 'A', 'C', 
        # 'T', 'G', and 'N'.
        if any([char.upper() not in ALLOWED_BASES for char in row["ReferenceSequence"]]):
            log_error(f"Reference sequence in row {i} has base other than 'A', 'C', 'T', 'G', or 'N'")
                        
        # Make sure that the same plate nickname and dual index plate name always 
        # go together.
        if row["PlateName"] not in index_to_nick_check:
            index_to_nick_check[row["PlateName"]] = row["IndexPlate"]
        else:
            if index_to_nick_check[row["PlateName"]] != row["IndexPlate"]:
                log_error(f"""Each plate nickname must be associated with a single dual index plate. In row {i}, nickname {row["PlateName"]} is duplicated with index plate {row["IndexPlate"]}""")
        
        # Additional row checks for a detailed file:
        if detailed_file:
            
            # Confirm that the named well is a real well
            if row["Well"] not in ALLOWED_WELLS:
                log_error(f"""Unexpected well in row {i}: {row["Well"]} of reference sequence file. 
                          Well must take form 'A##'""")
          
    
    # Additional wholistic checks for a detailed file: Make sure that there are 
    # no duplicate dual index plate-well combos
    if detailed_file:
        
        # Get all combos of dual index plate and well. 
        combos = [(row["IndexPlate"], row["Well"]) 
                  for _, row in ref_seqs_df.iterrows()]
        
        # Count the number of occurences of each. Identify those combos with 
        # more than one occurence
        counts = Counter(combos)
        duplicates = [combo for combo, count in counts.items() if count > 1]
        
        # Throw an error if there are duplicates
        if len(duplicates) > 0:
            log_error(f"""Each sample must have a unique index plate and well.
                     The following combinations of plate and well occured more than once in the refseq file:
                     {set(duplicates)}""")
           
# Write a function that checks the arguments passed into the command line prompt
def check_args(cl_args):
    
    # Confirm that the csv file with reference sequences exists
    if not os.path.exists(cl_args["refseq"]):
        
        # Write the error and terminate the program
        log_error(f"The file or folder does not exist: {cl_args['refseq']}")
        
    # Confirm that the reference sequence csv file is indeed a file
    if not os.path.isfile(cl_args["refseq"]):
        
        # Write the error and terminate the pgoram
        log_error(f"Specified refseq file is not a file: {cl_args['refseq']}")
    
    # Confirm that the specified folder containing sequencing results exists 
    if not os.path.exists(cl_args["folder"]):
        
        # Write the error and terminate the program
        log_error(f"The file or folder does not exist: {cl_args['folder']}")
    
    # Confirm that we have a folder or a file in the positional argument
    if os.path.isdir(cl_args["folder"]):
        
        # If we have a folder, make sure there are fastq or fastq.gz files inside of it
        all_files = glob(os.path.join(cl_args["folder"], "*.fastq*"))
        
        # If we can't find any files, report the error
        if len(all_files) == 0:
            
            # Write the error and terminate the program
            log_error("No fastq or fastq.gz files found in target folder.")
                
    # If we have a file, make sure the fastq_r flag is also thrown and another
    # file is given
    elif os.path.isfile(cl_args["folder"]):
        
        # Check to be sure that the fastq_r flag is thrown
        if cl_args["fastq_r"] == "":
            
            # Write the error and terminate the program
            log_error("--fastq-r must also be specified.")
            
        # Check to be sure that the fastq_r flag exists
        if not os.path.exists(cl_args["fastq_r"]):
            
            # Write the error and terminate the program
            log_error(f"Specified fastq-r file does not exist: {cl_args['fastq_r']}")
            
        # Check to be sure that the fastq_r flag contains a file
        if not os.path.isfile(cl_args["fastq_r"]):
            
            # Write the error and terminate the pgoram
            log_error(
                f"Specified fastq-r file is not a file: {cl_args['fastq_r']}")
            
        # Check to be sure that the fastq_r flag contains a file and that it is
        # either a fastq or fastq.gz    
        if "fastq" not in cl_args["fastq_r"] and "fastq.gz" not in cl_args["fastq_r"]:
            
            # Write the error and terminate the program
            log_error("The reverse file is neither a fastq nor fastq.gz file.")
            
        # Check to be sure that the fastq_f flag contains a file and that it is
        # either a fastq or fastq.gz    
        if "fastq" not in cl_args["folder"] and "fastq.gz" not in cl_args["folder"]:
            
            # Write the error and terminate the program
            log_error("The forward file is neither a fastq nor fastq.gz file.")
        
    # If we haven't found a file or a folder, throw an error
    else:
        
        # Log the error and terminate the program
        log_error("The positional argument is neither a file nor a folder.")
        
    # Confirm that the read length >0. Skip this if it is 'None' as this indicates
    # that we will be calculating it later on
    if cl_args["read_length"] is not None:
        if not cl_args["read_length"] > 0:
            
            # Write the error and terminate the program
            log_error("--read_length must be an integer greater than 0.")
    
    # Confirm that the Q-score cutoff >0
    if not cl_args["average_q_cutoff"] > 0:
        
        # Write the error and terminate the program
        log_error("--average_q_cutoff must be an integer greater than 0.")
        
    if not cl_args["bp_q_cutoff"] > 0:
        
        # Write the error and terminate the program
        log_error("--bp_q_cutoff must be an integer greater than 0.")
        
    # Confirm that the average q cutoff is less than or equal to the basepair
    # one
    if cl_args["average_q_cutoff"] > cl_args["bp_q_cutoff"]:
        log_error("--average_q_cutoff must be less than or equal to --bp_q_cutoff")
    
    # Confirm that the length cutoff is a float between 0 and 1
    if not ( 0 <= cl_args["length_cutoff"] <= 1):
        
        # Write the error and terminate the program
        log_error("--length_cutoff must be a float between 0 and 1.")
        
    # Confirm that the variable threshold is a float between 0 and 1
    if not ( 0 <= cl_args["variable_thresh"] <= 1):
        
        # Write the error and terminate the program
        log_error("--variable_thresh must be a float between 0 and 1.")
        
    # Confirm that the variable count is greater than or equal to 1
    if cl_args["variable_count"] < 1:
        
        # Write the error and terminate the program
        log_error("--variable_count must be 1 or higher.")
    
    # Confirm that the number of jobs is an integer between 1 and the number of 
    # cpu's available to the computer. Adjust the number to 1 if below 1 and to
    # the maximum number of cpus if above the maximum available.
    if cl_args["jobs"] < 1:
        
        # Raise the number of jobs to 1 and record warning
        cl_args["jobs"] = 1
        log_warning("--jobs must be greater than or equal to 1. Defaulting to 1 job.")
    
    elif not cl_args["jobs"] <= N_CPUS:
        
        # Lower the number of jobs to the maximum available and record warning
        cl_args["jobs"] = N_CPUS
        log_warning("--jobs must be less than or equal to the available number of cpus. Defaulting to max number of cpus.")
