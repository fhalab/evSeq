# Import third party packages
import pandas as pd
import os
import os.path

# Import ssSeqSupport variables
from . import CheckRefSeqs, CheckIndexMap
from . import AllowedWells
from . import Homedir, LogError


# Write a function that loads the reference sequence file, checks the integrity
# of the input information, and generates an appropriate reference sequence dataframe
def LoadRefSeq(args):
    
    # Load the reference sequence using pandas
    ref_seq_df = pd.read_csv(args["refseq"])
    
    # Check the validity of the reference sequence file
    CheckRefSeqs(ref_seq_df, args["detailed_refseq"])
    
    # If this is not a detailed reference sequence form, expand the reference 
    # sequence form to mimic one
    if not args["detailed_refseq"]:
        
        # Create a list in which to store reference sequence information
        expanded_ref_seqs = []
        
        # Loop over each row and construct a plate
        for _, row in ref_seq_df.iterrows():
            
            # Make a new list with well information appended
            expanded_ref_seqs.extend([[row["PlateName"], row["IndexPlate"],
                                      row["ReferenceSequence"], well] for 
                                      well in AllowedWells])
            
        # Construct a new ref_seq_df
        ref_seq_df = pd.DataFrame(expanded_ref_seqs, 
                                  columns = ("PlateName", "IndexPlate", 
                                             "ReferenceSequence", "Well"))
    
    # Output the ref_seq_df
    return ref_seq_df
        
# Write a function that loads and checks the contents of IndexMap.csv
def LoadDualInds():
    
    # Identify the expected location of IndexMap.csv
    index_map_loc = os.path.join(Homedir, "IndexMap.csv")
    
    # Check to make sure that IndexMap.csv can be found as a file
    if not os.path.exists(index_map_loc):
        print(index_map_loc)
        # Log an error
        LogError("Cannot find 'IndexMap.csv' in ssSeq/ssSeqSupport/")
        
    # Load the file
    index_df = pd.read_csv(index_map_loc)
    
    # Check the index file integrity
    CheckIndexMap(index_df)
    
    # Return the index_df
    return index_df

# Write a function which constructs global objects from the index_df and 
# ref_seq_df
def ConstructBCsToRefSeq(ref_seq_df, index_df):
    
    # Create a dictionary mapping index plate and well combos to reference 
    # sequences and plate nicknames
    inds_to_ref_plate = {(row["IndexPlate"], row["Well"]): (row["PlateName"], row["ReferenceSequence"])
                         for _, row in ref_seq_df.iterrows()}
    
    # Create a dictionary mapping index plate and well combos to forward and 
    # reverse barcodes
    inds_to_bcs = {(row["IndexPlate"], row["Well"]): (row["F-BC"], row["R-BC"])
                   for _, row in index_df.iterrows()}
        
    # Construct a dictionary linking forward and reverse barcodes to the 
    # reference sequence, plate name, and well
    bcs_to_refseq = {}
    error_message = "The following dual index plate-well pairs do not exist in IndexMap.csv:"
    error_found = False
    for key, (plate_nick, refseq) in inds_to_ref_plate.items():
        
        # Check to make sure there is a matching plate-well location in 
        # the barcode plate. Record if there is not one.
        if key not in inds_to_bcs:
            
            # Record that an error was found
            error_found = True
            
            # Report the key that doesn't match
            error_message += "\n{}-{}".format(*key)
            
        # If an error has been identified, just continue
        if error_found:
            continue
            
        # Pull the barcodes associated with the plate-well pair
        fbc, rbc = inds_to_bcs[key]        
        
        # Deconstruct the key
        ind_plate, well = key
    
        # Add to the dictionary
        bcs_to_refseq[(fbc, rbc)] = (ind_plate, well, plate_nick, refseq)
            
    # Log the error if we couldn't find the dual index plate-well combo
    if error_found:
        LogError(error_message)
        
    # Return bcs_to_refseq and the barcode lengths
    return bcs_to_refseq