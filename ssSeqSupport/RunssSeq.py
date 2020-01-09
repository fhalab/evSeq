# Import third party modules
from tqdm import tqdm
from scipy.stats import mode
import numpy as np
import os
import os.path
import re
import gzip

# # Import ssSeqSupport modules
from . import LoadRefSeq, LoadDualInds, ConstructBCsToRefSeq
from . import CheckArgs
from . import LogInputFiles, LogError, LogInfo
from . import SeqPair
from . import CreateID
from . import GenerateReadQualChart
from . import PlateObjects

# Write a function that performs all of ssSeq
# This will be called from the command line interface to run the full operation.
def RunssSeq(args):
    
    """
    Performs the full ssSeq analysis based on the arguments passed into 
    the command line.
    
    Parameters
    ----------
    args: dict: Dictionary resulting from calling 'parser.parse_args()' in the 
        initial setup of the run.
    """

    # Build the output directories
    _BuildOutputDirs(args)
            
    # Check the argument inputs. 
    all_files, folder = CheckArgs(args)
    
    # If we are working with a folder, identify matching files
    if folder:
        filepairs, unmatched_files = _FindMatches(all_files)
    
    # Otherwise, just package matching files
    else:
        filepairs = {args["folder"]: args["fastq_r"]}
        unmatched_files = []
    
    # Log the matched and unmatched files
    LogInputFiles(filepairs, unmatched_files)
        
    # Loop over the different filepairs
    for i, filepair in tqdm(enumerate(filepairs.items()), desc = "Run #", 
                         total = len(filepairs), leave = False):
        
        # Build sequencing pairs
        seq_pairs = _BuildSeqPairs(*filepair)
        
        # Set the appropriate read length in args if the default is used
        if args["read_length"] == None:
            args["read_length"] = _SetReadLength(seq_pairs)
        
        # Analyze the seq_pairs. 
        _AnalyzeSeqPairs(seq_pairs, filepair, args)
        
        # Process the seq_pairs only if the --analysis_only flag is not thrown
        if not args["analysis_only"]:

            # Map barcodes to reference sequences
            BcsToRefSeq = ConstructBCsToRefSeq(LoadRefSeq(args),
                                               LoadDualInds())

            # Process sequence pairs
            _ProcessSeqPairs(seq_pairs, BcsToRefSeq, args, filepair, i)
    
# Write a function for setting the read length argument when it is not passed in
def _SetReadLength(seq_pairs):
    
    # Get the length of every sequence pair
    seq_lengths = np.array([[pair.f_len, pair.r_len] for pair in seq_pairs])
    seq_lengths = seq_lengths.flatten()
    
    # Get the mode of the sequence lengths. This is what we use as the sequence
    # length
    read_length = mode(seq_lengths)[0][0]
    
    # Record the identified read length
    LogInfo("The identified read length for this run was: {}".format(read_length))
    
    # Return the read length
    return read_length

# Write a function that builds the output directory structure
def _BuildOutputDirs(args):
    
    # Build the folder structure if it does not exist
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
        
    # Build the summaries folder only if we are not in ts mode
    summary_dir = os.path.join(args["output"], "Summaries/")
    os.mkdir(summary_dir)
    
    # Build the read qualities folder
    qual_dir = os.path.join(args["output"], "Qualities/")
    os.mkdir(qual_dir)
    
    # Build the heatmaps folder
    heatmap_dir = os.path.join(args["output"], "Platemaps/")
    os.mkdir(heatmap_dir)
    
    # If we are in ts mode, build additional directories
    if args["troubleshoot"]:
        extra_dirs = [os.path.join(args["output"], loc) for loc in
                      ["Alignments", "AACountsFrequencies",
                       "BPCountsFrequencies", "ConsensusSequences"]]
        for directory in extra_dirs:
            os.mkdir(directory)
    
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

# Write a function that builds sequence pairs from the forward and reverse
# read files
def _BuildSeqPairs(ForwardReads_Filepath, ReverseReads_Filepath):
    
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

    # Loop over the forward file and start building objects
    with forward_parser as f:

        # Load the file to memory and count how many lines are in it
        file_list = list(f)

        # Create a dictionary so we can find seq_pairs later
        id_to_pair = {}

        # Create an index for use in constructing individual sequence objects
        block_ind = 0

        # Loop over all lines, feeding blocks of 4 into a seq_pair object
        info_block = [None, None, None, None]
        for i, line in enumerate(tqdm(file_list, position = 1, leave = False,
                                      desc = "Preparing forward reads"), 1):

            # Strip whitespace from the line
            line = line.strip()

            # Add to the info_block
            info_block[block_ind] = line
            
            # Update the block_ind
            block_ind += 1

            # If this is a line divisble by 4, create a new seq_pair object and append to the growing list
            if i % 4 == 0:

                # Add information to a dictionary
                seq_pair = SeqPair(info_block)
                id_to_pair[seq_pair.id] = seq_pair

                # Reset the info block list and block_ind
                info_block = [None, None, None, None]
                block_ind = 0

    # Loop over the reverse file and append partners to each file
    with reverse_parser as f:

        # Load the file to memory
        file_list = list(f)

        # Create an index for use in constructing individual sequence objects
        block_ind = 0

        # Loop over each line in the file.
        info_block = [None, None, None, None]
        for i, line in enumerate(tqdm(file_list, position = 1, leave = False,
                                      desc = "Preparing reverse reads"), 1):

            # Strip whitespace from line
            line = line.strip()

            # Determine what line we're on
            line_check = i % 4

            # Add to the info block
            info_block[block_ind] = line
            
            # Update the block_ind
            block_ind += 1

            # Identify if this is a line divisible by 4. If so, append to the appropriate seq_pair
            if line_check == 0:

                # Identify the uid
                uid = CreateID(id_line)

                # If there is a match for the uid, append it to the appropriate partner
                if uid in id_to_pair:
                    id_to_pair[uid].attach_partner(info_block)

                # Reset the info block
                info_block = [None, None, None, None]
                block_ind = 0

            # Identify if this is the first line. If so, this is the id line.
            elif line_check == 1:

                # Set the line as the id line
                id_line = line
            
    # Return seq_pairs as a list
    return list(id_to_pair.values())

# Write a function analyzes sequencing pairs only
def _AnalyzeSeqPairs(seq_pairs, filepair, args):
    
    # Get the path to the Q-score save location
    qual_dir = os.path.join(args["output"], "Qualities/")
    
    # Get the basenames of the two files
    f_basename = os.path.splitext(os.path.basename(filepair[0]))[0]
    r_basename = os.path.splitext(os.path.basename(filepair[1]))[0]
    
    # Get the average read qualities of the forward reads
    mean_f_qual_scores = [int(np.mean(seq_pair.f_qual_scores)) for seq_pair in seq_pairs]
    f_qual_counts = np.unique(mean_f_qual_scores, return_counts=True)
    
    # save the reverse read qualities
    mean_r_qual_scores = [int(np.mean(seq_pair.r_qual_scores)) for seq_pair in seq_pairs]
    r_qual_counts = np.unique(mean_r_qual_scores, return_counts=True)
    
    # generate quality score histogram
    counts = (f_qual_counts, r_qual_counts)
    GenerateReadQualChart(counts, "{}{}-{}-ReadQualPlot.html".format(qual_dir,
                                                                     f_basename,
                                                                     r_basename))

# Write a function that processes sequencing pairs
def _ProcessSeqPairs(seq_pairs, BcsToRefSeq, args, filepair, combo_ind):

    # Parse each pair and separate into plates and wells. Discard orphan pairs.
    well_pairs = {}
    for pair in seq_pairs:
            
        # Determine if the combination of forward and reverse barcode has been determined.
        # If the combo has not been seen yet, create a list in the dictionary
        well_ID = (pair.f_barcode, pair.r_barcode)
        if well_ID not in well_pairs:
            well_pairs[well_ID] = [pair]

        # If the combo has been seen, add to the list that's growing
        else:
            well_pairs[well_ID].append(pair)

    # Create wells for each list in the well_pairs dictionary
    all_wells = [PlateObjects.Well(pair_list, args["read_length"], BcsToRefSeq) 
                 for pair_list in well_pairs.values()]
    true_wells = [well for well in all_wells if well.real_well is True]

    # Parse each well and separate wells by plate
    plate_wells = {}
    for well in true_wells:

        # Determine if the plate has been seen yet. If not, create a list
        # in the dictionary
        if well.plate not in plate_wells:
            plate_wells[well.plate] = [well]

        # If the plate has been seen, just append to the list.
        else:
            plate_wells[well.plate].append(well)

    # Generate plates
    plates = [PlateObjects.Plate(well_list) for well_list in plate_wells.values()]
    n_plates = len(plates)

    # Loop over each plate
    for i, plate in tqdm(enumerate(plates), total = n_plates, 
                         desc = "Processing plates", position = 1,
                         leave = False):
        desc = "{} of {}".format(i + 1, n_plates)
        plate.process(args, desc, filepair, combo_ind)