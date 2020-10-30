
# Write a function that matches forward and reverse reads in a passed in folder
def _FindMatches(seqfiles):
    
    # Define a generic regex to identify forward and reverse reads in the folder
    f_regex = re.compile("(.*)_R1_.*")
    r_regex = re.compile("(.*)_R2_.*")
    
    # Extract all forward filenames from the seqfiles
    f_searches = [f_regex.search(filename) for filename in seqfiles]
    f_filenames = [filename for filename, test in zip(seqfiles, f_searches) if test]
    f_filenames_set = set(f_filenames)
    
    # Identify remaining files
    remaining_filenames = [filename for filename in seqfiles 
                           if filename not in f_filenames_set]
    
    # Create a dictionary storing all forward files. Make the values None for now
    filepairs = {filename: None for filename in f_filenames}
    
    # Create a dictionary linking a root name to f_filename
    root_to_f = {f_regex.search(filename).group(1): filename for 
                 filename in f_filenames}
    
    # Loop over the remaining filenames and match to the forward files. Also record
    # those filenames with no match
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
    if len(final_filepairs) == 0:
        LogError("Could not match forward and reverse files. Make sure they have appropriate names.")
    
    # Return the dictionary matching files as well as the unmatched files
    return final_filepairs, unmatched_files

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

# Write a function that identifies the positions of "NNN" in a reference sequence
def FindNNN(reference_sequence):

    # Define a list which will record where the variable positions start in the 
    # reference sequence
    var_sites = []
    
    # Define a list which will record where any "N" is found in the reference
    # sequence
    N_sites = []
    
    # Define a variable to track the number of variable sites found
    N_found = 0

    # Create a variable to record whether or not we are in series
    in_series = True

    # Loop over each base in the reference
    for i, char in enumerate(reference_sequence):

        # Find each N. If it is the first in a series of 3, report it as the start
        # site for a variable position
        if char=="N" and N_found % 3==0:
            
            # Record that this is an 'N' that we found, and that it is the start
            # of a codon
            var_sites.append(i)
            N_sites.append(i)
            N_found += 1
            
            # Record that this is the latest value of N found
            latest_N = i

        elif char=="N":
            
            # Record that this is an 'N' that we found
            N_sites.append(i)
            N_found += 1
            
            # Check to make sure we are in series with the previous N found
            if i - latest_N != 1:
                in_series = False
                
            # Update latest_N
            latest_N += 1

    # Confirm that N was found in multiples of 3
    found_in_3 = True if N_found % 3 == 0 else False
    
    # Confirm that N was found in codon format
    codon_format = in_series and found_in_3
    
    # If there are no variable sites, then throw an error
    n_var_sites = len(var_sites)
    if n_var_sites == 0:
        LogError("No variable sites detected in one of the forward or reverse reference sequences.")
    
    # Check to be sure that the variable sites all exist in the same reference
    # frame. This is only important when we have more than one codon present.
    if n_var_sites > 1:
        
        # What is the distance between the location of each codon?
        diffs = np.diff(var_sites)
        
        # Are all differences divisible by 3? If they are all in the same frame
        # then they should be
        in_frame = np.all(diffs%3 == 0)
        
        # Throw an error if not all are in frame
        if not in_frame:
            LogError("Specified variable positons are not in the same reading frame. Aborting run.")
    
    # Return the variable sites, the number of variable sites, and whether or 
    # not N was included in codon-format (multiples of 3 in series) 
    return var_sites, n_var_sites, codon_format

# Associate reference sequences with barcodes
def load_refseq(ref_seq_loc):
    
    # Load the index map and reference sequence
    index_map = pd.read_csv("/home/brucejwittmann/GitRepos/ssSeq/ssSeqSupport/IndexMap.csv")
    ref_seq_crude = pd.read_csv("/home/brucejwittmann/GitRepos/ssSeq/AlignmentDev/TestData/20200205_ssSeq/RefSeqs.csv")

    # Expand each reference sequence
    updated_ref_array = []
    for row in ref_seq_crude.itertuples(index = False):
        updated_ref_array.extend([[row.PlateName, row.IndexPlate, well, row.ReferenceSequence, row.InFrameBase]
                                 for well in ALLOWED_WELLS])

    # Define the complete reference sequence dataframe
    complete_ref_seq = pd.DataFrame(updated_ref_array, columns = ("PlateName", "IndexPlate", "Well", "ReferenceSequence", "InFrameBase"))

    # Join on plate and well
    merged_dfs = complete_ref_seq.merge(index_map, on = ("IndexPlate", "Well"))

    # Map barcode to reference sequence, plate, and well
    warnings.warn("Did not yet implement calculation of variable bases")
    bc_to_ref_plate_well = {(row.FBC, row.RBC): {"IndexPlate": row.IndexPlate,
                                                 "PlateNickname": row.PlateName,
                                                 "Well": row.Well,
                                                 "ReferenceSequence": row.ReferenceSequence,
                                                "InFrameBase": row.InFrameBase,
                                                "ExpectedVariablePositions": np.array([], dtype = int),
                                                "BpIndStart": 1,
                                                "AAIndStart": 1}
                           for row in merged_dfs.itertuples(index = False)}
    
    return bc_to_ref_plate_well

# Write a function for loading and pairing fastq files
def load_fastq(f_loc, r_loc):

    # First check if file paths are gzipped or not
    gzipped_f = True if "fastq.gz" in ForwardReads_Filepath else False
    gzipped_r = True if "fastq.gz" in ReverseReads_Filepath else False

    # Assign proper parser
    if gzipped_f:
        forward_parser = gzip.open(ForwardReads_Filepath, "rt")
    else:
        forward_parser = open(ForwardReads_Filepath, "r")

    if gzipped_r:
        reverse_parser = gzip.open(ReverseReads_Filepath, "rt")
    else:
        reverse_parser = open(ReverseReads_Filepath, "r")


    # Create a dictionary that links id to sequence object
    id_to_reads = {}
    print("Loading forward reads...")
    all_f_recs = list(SeqIO.parse(f_loc, "fastq"))
    for f_record in all_f_recs:
        temp_record = SeqPair()
        temp_record.assign_f(f_record)
        id_to_reads[f_record.id] = temp_record
    
    # Associate reverse reads with the forward
    print("Loading reverse reads...")
    all_r_recs = list(SeqIO.parse(r_loc, "fastq"))
    for r_record in all_r_recs:

        # If there is no partern in id_to_reads, create a new object 
        # and continue
        if r_record.id not in id_to_reads:
            temp_record = SeqPair()
            temp_record.assign_r(r_record)
            id_to_reads[r_record.id] = temp_record

        # Otherwise, attach the reverse record
        else:
            id_to_reads[r_record.id].assign_r(r_record)
            
    # Only keep records that have a partner
    return tuple(id_to_reads.values())