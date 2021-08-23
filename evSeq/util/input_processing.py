
# Import evSeq objects
from .globals import ALLOWED_WELLS, UTILDIR, ADAPTER_F, ADAPTER_R
from .logging import log_error, log_input_file
from .input_validation import check_ref_seqs, check_index_map

# Import other required modules
import re
import os
import pandas as pd
import numpy as np
from glob import glob
from Bio.Seq import Seq

def build_output_dirs(cl_args):
    """Builds the output directory structure."""
    
    # Build the folder structure if it does not exist
    if not os.path.exists(cl_args["output"]):
        os.makedirs(cl_args["output"])
    
    # Set up required folders based on type of run
    if cl_args["analysis_only"]:
        folders = ["Qualities"]

    elif cl_args["only_parse_fastqs"]:
        folders = [
            "Qualities",
            os.path.join("ParsedFilteredFastqs", "F"),
            os.path.join("ParsedFilteredFastqs", "R"),
        ]
    else:
        # Standard run
        folders = [
            "Qualities",
            "OutputCounts",
            "Platemaps"
        ]
    
    # Add folders for non-standard run
    if cl_args["keep_parsed_fastqs"]:
        folders = folders + [
            os.path.join("ParsedFilteredFastqs", "F"),
            os.path.join("ParsedFilteredFastqs", "R")
        ]
    if cl_args["return_alignments"]:
        folders.append("Alignments")

    # Build required folders
    for folder in folders:
        os.makedirs(os.path.join(cl_args["output"], folder), exist_ok=True)

def find_matches(input_folder):
    """Matches forward and reverse reads in a passed in folder."""
    
    # Glob all potential sequencing files in the directory
    seqfiles = glob(os.path.join(input_folder, "*fastq*"))
    
    # Define a generic regex to identify forward and reverse reads in the
    # folder
    f_regex = re.compile("(.*)_R1_.*")
    r_regex = re.compile("(.*)_R2_.*")
    
    # Extract all forward filenames from the seqfiles
    f_searches = [f_regex.search(filename) for filename in seqfiles]
    f_filenames = [filename for filename, test in zip(seqfiles, f_searches) if test]
    f_filenames_set = set(f_filenames)
    
    # Identify remaining files
    remaining_filenames = [filename for filename in seqfiles 
                           if filename not in f_filenames_set]
    
    # Create a dictionary storing all forward files. Make the values None
    # for now
    filepairs = {filename: None for filename in f_filenames}
    
    # Create a dictionary linking a root name to f_filename
    root_to_f = {f_regex.search(filename).group(1): filename for 
                 filename in f_filenames}
    
    # Loop over the remaining filenames and match to the forward files. Also
    # record those filenames with no match
    unmatched_files = []
    for filename in remaining_filenames:
        
        # Search to see if this is a reverse file
        reverse_search = r_regex.search(filename)
        
        # Get the root name check to see if there is a matching forward file
        if reverse_search:
                    
            # Get the rootname
            rootname = reverse_search.group(1)
            
            # Confirm that the rootname is present in root_to_f
            if rootname in root_to_f:
            
                # Match the reverse filename to the forward filename
                filepairs[root_to_f[rootname]] = filename
                
            # If we can't find a match, record this as an unmatched file
            else:
                unmatched_files.append(filename) 
            
        # If we can't find a match, record this as an unmatched file
        else:
            unmatched_files.append(filename) 
            
    # Check to see which filepairs still have 'None' as their value. Record 
    # these as unmatched files
    final_filepairs = {}
    for filekey, fileval in filepairs.items():
        
        # If fileval is None, record the filekey as an unmatched file
        if fileval is None:
            unmatched_files.append(filekey)
            
        # Otherwise add to the final filepair dictionary
        else:
            final_filepairs[filekey] = fileval
    
    # If final_filepairs is 0, throw and error and terminate the program
    n_pairs = len(final_filepairs)
    if n_pairs == 0:
        log_error("Could not match forward and reverse files. "
                  "Make sure they have appropriate names.")
    
    # If more than 1 pair is identified, throw and error and terminate the 
    # program
    if n_pairs > 1:
        log_error("More than 1 pair of sequencing files found in input directory.")
    
    # Return the dictionary matching files as well as the unmatched files
    forward_file, reverse_file = list(final_filepairs.items())[0]
    return forward_file, reverse_file, unmatched_files

def load_dual_inds():
    """Loads and checks the contents of index_map.csv."""
    
    # Identify the expected location of index_map.csv
    index_map_loc = os.path.join(UTILDIR, "index_map.csv")
    
    # Check to make sure that index_map.csv can be found as a file
    if not os.path.exists(index_map_loc):
        print(index_map_loc)
        # Log an error
        log_error("Cannot find 'index_map.csv'")
        
    # Load the file
    index_df = pd.read_csv(index_map_loc)
    
    # Check the index file integrity
    check_index_map(index_df)
    
    # Return the index_df
    return index_df

def construct_ref_seq(refseq_df):
    """
    Constructs a reference sequence based on the input sequences, and
    adjusts indexes and parameters accordingly.

    Returns four lists:
    -------------------
    refseqs: list of strings (DNA seqs)
        the new ReferenceSequence for each row in the refseq file.
    new_frame: list of ints, ranging [0, 2]
        the new FrameDistance, adjusted for the length of the FPrimer
        seed region.
    new_bp: list of ints
        new BpStartInd, adjusted for the length of the FPrimer seed 
        region.
    new_aa: list of ints; usually ~ (new_bp // 3)
        new AaStartInd, adjusted for the length of the FPrimer seed
        region.
    """
    # Get the forward primer and remove adapter
    f_primers = [seq.upper() for seq in refseq_df["FPrimer"]]
    adapterless_f_primers = [seq.replace(ADAPTER_F, '') for seq in f_primers]

    # Same for reverse
    r_primers = [seq.upper() for seq in refseq_df["RPrimer"]]
    adapterless_r_primers = [seq.replace(ADAPTER_R, '') for seq in r_primers]

    # Reverse complement the reverse primer
    r_seqs = [Seq(seq).reverse_complement() for seq in adapterless_r_primers]

    # Create a reference sequence from the primers and variable region
    var_regions = [seq.upper() for seq in refseq_df["VariableRegion"]]
    refseqs = [None for _ in var_regions]
    zipped = enumerate(zip(adapterless_f_primers, var_regions, r_seqs))
    for i, (f, var, r) in zipped:
        refseqs[i] = str(f + var + r)

    # Adjust frame distance relative to refseq
    bp_offsets = [len(seq) for seq in adapterless_f_primers]

    # For a given N (bp offset), the FrameDistance changes periodically
    # according to (I + N) % 3
    new_frame = [(I + N) % 3
                 for N, I in zip(bp_offsets, refseq_df['FrameDistance'])]

    # Now adjust the Bp and Aa IndStart to match refseq
    new_bp = [I - N
              for N, I in zip(bp_offsets, refseq_df['BpIndStart'])]
    new_aa = [I - ((N+F) // 3)
              for N, F, I, in zip(
                  bp_offsets,
                  refseq_df['FrameDistance'],
                  refseq_df['AaIndStart'],
              )]

    return refseqs, new_frame, new_bp, new_aa

def load_ref_seq(cl_args):
    """Loads the reference sequence file, checks the integrity of the
    input information, and generates an appropriate reference sequence
    dataframe.
    """
    
    # Load the reference sequence using pandas
    refseq_df = pd.read_csv(cl_args["refseq"])
    
    # Check the validity of the reference sequence file
    check_ref_seqs(refseq_df, cl_args["detailed_refseq"])

    # Construct reference sequence
    refseqs, new_frame, new_bp, new_aa = construct_ref_seq(refseq_df)

    # Save new values in dataframe
    refseq_df["ReferenceSequence"] = refseqs
    refseq_df['FrameDistance'] = new_frame
    refseq_df['BpIndStart'] = new_bp
    refseq_df['AaIndStart'] = new_aa
    
    # If this is not a detailed reference sequence, expand it
    if not cl_args["detailed_refseq"]:
    
        # Expand to a detailed reference sequence dataframe
        updated_ref_array = []
        for row in refseq_df.itertuples(index = False):
            updated_ref_array.extend([[row.PlateName, row.IndexPlate, well,
                                       row.FPrimer, row.RPrimer,
                                       row.VariableRegion,
                                       row.ReferenceSequence, row.FrameDistance,
                                       row.BpIndStart, row.AaIndStart]
                                      for well in ALLOWED_WELLS])

        # Define the complete reference sequence dataframe
        refseq_df = pd.DataFrame(updated_ref_array,
                                 columns = ("PlateName", "IndexPlate", "Well",
                                            "FPrimer", "RPrimer",
                                            "VariableRegion",
                                            "ReferenceSequence", "FrameDistance",
                                            "BpIndStart", "AaIndStart"))
    
    # Output the ref_seq_df
    return refseq_df

def find_codons_variable_positions(refseq, inframe_bp, plate, well):
    """Finds variable positions in a reference sequence."""
    
    # Get the number of codons and number of variable positions
    n_codons = refseq.count("NNN")
    n_positions = refseq.count("N")
    
    # Check to be sure the number of positions is divisible by 3 and divides
    # to the number of codons
    if (n_positions % 3) != 0:
        log_error(f"Error for {plate}-{well} refseq:"
                  "Must enter `N` in groups of 3 to signify codon.")
    assert (n_positions / 3) == n_codons, "Mismatch in number of codons and number of variable positions."

    # If no codons are found, return an empty array
    if n_codons == 0:
        return np.array([], dtype = int)
    
    # Split the reference sequence
    split_seqs = refseq.split("NNN")

    # Loop over the splits and make checks
    variable_positions = []
    variable_counter = 0
    for frag_ind, fragment in enumerate(split_seqs[:-1]):

        # Make sure fragments other than the first are divisible by 3, if not, then
        # the sections are not in frame
        frag_len = len(fragment)
        if (frag_ind > 0) and (not (frag_len % 3 == 0)):
            log_error(f"Error for {plate}-{well} refseq:"
                      "Specified variable positions are not in the same reading frame.")

        # Extend the variable positions. It will be extended by 
        # the fragment length
        variable_counter += frag_len
        for _ in range(3):
            variable_positions.append(variable_counter)
            variable_counter += 1

    # Convert to an array. Make sure everything is sorted.
    variable_positions = np.array(variable_positions)
    variable_positions.sort()

    # Make sure that all variable positions are "N" and that we found the appropriate
    # number of codons
    assert len(variable_positions) == n_positions, "Missing variable positions"
    assert all(refseq[pos] == "N" for pos in variable_positions), "Error in variable position calculation"
        
    # Validate that the provided in frame bp is correct
    if (variable_positions[0] % 3) != (inframe_bp):
        log_error(f"Error for {plate}-{well} refseq:"
                  "The provided in frame base does not match the frame of the provided variable codons.")
    
    return variable_positions

def construct_bcs_to_refseq(refseq_df, index_df):
    """Constructs global objects from the index_df and ref_seq_df."""
    
    # Join on plate and well
    merged_dfs = refseq_df.merge(index_df, on=("IndexPlate", "Well"))

    # Make sure the length of the refseq df is the same as the length of the 
    # merged dfs. We do not want to have lost anything in the merge.
    if len(merged_dfs) != len(refseq_df):
        log_error("Merging refseq_df and index_df failed. Values lost in merge.")

    # For each reference sequence, link barcode to plate information. Also 
    # calculate any expected mutagenesis positions
    bc_to_ref_plate_well = {}
    for row in merged_dfs.itertuples(index = False):

        # Find variable positions
        variable_positions = find_codons_variable_positions(
            row.ReferenceSequence,
            row.FrameDistance,
            row.IndexPlate,
            row.Well,
        )

        # Map barcode to reference sequence, plate, and well
        bc_to_ref_plate_well[(row.FBC, row.RBC)] = {
            "IndexPlate": row.IndexPlate,
            "PlateName": row.PlateName,
            "Well": row.Well,
            "FPrimer": row.FPrimer,
            "RPrimer": row.RPrimer,
            "VariableRegion": row.VariableRegion,
            "ReferenceSequence": row.ReferenceSequence,
            "FrameDistance": row.FrameDistance,
            "ExpectedVariablePositions": variable_positions,
            "BpIndStart": row.BpIndStart,
            "AAIndStart": row.AaIndStart,
        }
        
    # Return the mapping of barcode pair to well information
    return bc_to_ref_plate_well

def load_all(cl_args):
    """Loads everything. This wraps all functions defined above."""
    
    # Find matching forward and reverse reads in a folder if the files are not
    # explicitly passed in together
    if cl_args["fastq_r"] == "":
        forward_file, reverse_file, unmatched_files = find_matches(cl_args["folder"])
    
    # Otherwise, assign files
    else:
        forward_file = cl_args["folder"]
        reverse_file = cl_args["fastq_r"]
        unmatched_files = []
    
    # Log the identified input files
    log_input_file(forward_file, reverse_file, unmatched_files)
    
    # Load the index dataframe and reference sequence files
    index_df = load_dual_inds()
    refseq_df = load_ref_seq(cl_args)
    
    # Merge the two dataframes and generate the dataframe info needed to be 
    # passed in to all wells
    bc_to_ref_plate_well = construct_bcs_to_refseq(refseq_df, index_df)
    
    return forward_file, reverse_file, bc_to_ref_plate_well
