# Import third party modules
import numpy as np
import os

# Import ssSeqSupport modules
from . import (IdParser, ReverseCompDict, CodonTable, AdapterLengthF,
               AdapterLengthR)
from . import LogError

# Write a function for pulling information from the id line
def GetBlockInfo(id_line):

    # Parse the block with a regex
    return IdParser.search(id_line).groups()

# Write a function to generate ids from a single id line
def CreateID(id_line):

        # Get the id information for the block
        (instrument_name, _, _, lane, tile, x_coord, y_coord, _, _, _,
         sample_number) = GetBlockInfo(id_line)

        # Create a unique id for the pair that can pair ends
        return (instrument_name, lane, tile, x_coord, y_coord, sample_number)

# Write a function that returns the reverse complement of a sequence
def ReverseComplement(seq):

    # Loop through the sequence in reverse and translate
    return "".join(ReverseCompDict[char] for char in reversed(seq))

# Write a function for translating sequences
def Translate(seq, start_ind):

    # Get the number of codons in the sequence
    n_codons = np.floor((len(seq) - start_ind)/3)

    # Identify all codons
    codons = [seq[int(start_ind + i*3): int(start_ind + (i+1)*3)] for i in range(int(n_codons))]

    # Translate all codons
    translation = []
    for codon in codons:

        # If this is "NNN" we give a question mark
        if "N" in codon:
            translation.append("?")
        else:
            translation.append(CodonTable[codon])

    # Return the translation
    return "".join(translation)

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
    if len(var_sites) == 0:
        LogError("No variable sites detected in one of the forward or reverse reference sequences.")
    
    # Return the variable sites, the number of variable sites, and whether or 
    # not N was included in codon-format (multiples of 3 in series) 
    return var_sites, len(var_sites), codon_format

# Write a function that builds the output directory structure
def BuildOutputDirs(args):
    
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