# Import required modules
import pandas as pd
import numpy as np
import re

###############################################################################
########################## Functions for use in classes #######################

# Write a function for pulling information from the id line
def get_block_info(id_line):

    # Parse the block with a regex
    return id_parser.search(id_line).groups()

# Write a function to generate ids from a single id line
def create_id(id_line):

        # Get the id information for the block
        (instrument_name, runid, flowcellid, lane, tile,
         x_coord, y_coord, pair_number, filtered_bin,
         control_bin, sample_number) = get_block_info(id_line)

        # Create a unique id for the pair that can pair ends
        return (instrument_name, lane, tile, x_coord, y_coord, sample_number)

# Write a function that returns the reverse complement of a sequence
def reverse_complement(seq):

    # Loop through the sequence in reverse and translate
    return "".join(reverse_comp_dict[char] for char in reversed(seq))

# Write a function for translating sequences
def translate(seq, start_ind, codon_table):

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
            translation.append(codon_table[codon])

    # Return the translation
    return "".join(translation)

# Define a class that represents a reference sequence
class ref_seq():

    # Define the initialization. This takes the reference sequence, identifies
    # variable sites, and translates where possible.
    def __init__(self, reference_sequence, codon_table):

        # Assign the reference sequence as a class attribute as well as its length
        reference_sequence = reference_sequence.upper()
        self.seq = reference_sequence
        self.seq_length = len(reference_sequence)

        # Define lists for collecting information
        self.var_sites = []

        # Define variables that will be used
        N_found = 0

        # Loop over each base in the reference
        for i, char in enumerate(reference_sequence):

            # Find the each N. This tells us the frame.
            if char=="N" and N_found==0:
                self.first_N = i
                N_found += 1
                self.var_sites.append(i)

            elif char=="N" and N_found % 3==0:
                self.var_sites.append(i)
                N_found += 1

            elif char=="N":
                N_found += 1

        # Identify the number of variable sites
        self.n_var_sites = len(self.var_sites)

        # Identify the translation start site
        self.trans_start = self.var_sites[0] % 3

        # Translate the reference sequence
        self.translated_seq = translate(reference_sequence, self.trans_start,
                                        codon_table)

        # Record the length of the translation
        self.trans_length = len(self.translated_seq)

        # Record where we can find the amino acids of interest
        self.var_aa_sites = [int(np.floor(site / 3)) for site in self.var_sites]

        # Write the bp and aa sequences as lists
        self.aas_as_list = [char for char in self.translated_seq]
        self.bps_as_list = [char for char in self.seq]

###############################################################################
##################### Load Relevant Files and Parameters ######################

# Create a regex to parse global parameters
global_parser = re.compile(".+: (.+)")

# Parse the global parameters file
with open("./GlobalParams.txt", "r") as f:

    # Convert the file to a list
    param_list = [line.strip() for line in f if "-" in line]

    # Extract variables from the file
    (RefSeq_Filepath,
     ForwardReads_Filepath,
     ReverseReads_Filepath,
     ts_mode,
     n_jobs,
     q_cutoff,
     alignment_cutoff,
     output_location,
     id_regex) = [global_parser.search(param_list[i]).group(1)
                               for i in range(len(param_list))]

# Convert the q-score cutoff to an integer
q_cutoff = int(q_cutoff)

# Convert the alignment cutoff to a float
alignment_cutoff = float(alignment_cutoff)

# Convert n_jobs to an integer
if n_jobs == "None":
    n_jobs = None
else:
    n_jobs = int(n_jobs)

# Convert ts_mode to a boolean
if ts_mode == "True":
    ts_mode = True
elif ts_mode == "False":
    ts_mode = False
else:
    print("ts_mode must be 'True' or 'False'. Aborting...")
    quit()

# Load in the barcode files and the allowed codons
DualInds = pd.read_excel("./ssSeq/DnaParams.xlsx", sheet_name = "DualIndices")
CodonTable_df = pd.read_excel("./ssSeq/DnaParams.xlsx", sheet_name = "CodonTable")
BPtoInd_df = pd.read_excel("./ssSeq/DnaParams.xlsx", sheet_name = "BPtoInd")
AAtoInd_df = pd.read_excel("./ssSeq/DnaParams.xlsx", sheet_name = "AAtoInd")

# Load the reference sequences
ref_seqs = pd.read_csv(RefSeq_Filepath)

###############################################################################
############## #Build Globals from the Loaded Parameters ######################
# Create a regex for parsing the id line
id_parser = re.compile(id_regex)

# Define a codon table
codon_table = {codon: aa for codon, aa in CodonTable_df.itertuples(index=False)}

# Get lists of unique dual barcodes and unique reference sequences
unique_well_BCs = list({(fbc, rbc) for fbc, rbc in DualInds.loc[:, ["F-BC", "R-BC"]].itertuples(index=False)})
unique_f_ref_seqs = ref_seqs["F-RefSeq"].unique().tolist()
unique_r_ref_seqs = ref_seqs["R-RefSeq"].unique().tolist()

# Construct a table to reference sequences to barcodes
lookup = ref_seqs.merge(DualInds, how = "left", on=["IndexPlate", "Well"])

# Take the barcode files and make dictionaries that link barcodes to plates and wells
BC_to_well = {(f_BC, r_BC): (plate, well) for plate, well, f_BC, r_BC in DualInds.itertuples(index=False)}

# Take the lookup table and make dictionaries that convert forward and reverse
# barcodes to reference sequences
BCs_to_refseq = {(f_BC, r_BC): (ref_seq(f_ref_seq, codon_table),
                                ref_seq(r_ref_seq, codon_table),
                                plate_name, n_variable_sites) for
                 f_BC, r_BC, f_ref_seq, r_ref_seq, plate_name, n_variable_sites in
                 lookup.loc[:, ["F-BC", "R-BC", "F-RefSeq",
                 "R-RefSeq", "PlateName", "NumberOfVariablePositions"]].itertuples(index=False)}

# Define barcode length
f_bc_length = len(unique_well_BCs[0][0])
r_bc_length = len(unique_well_BCs[0][1])

# Define a dictionary which links basepairs to indices and the reverse.
bp_to_ind = {bp: ind for bp, ind in BPtoInd_df.itertuples(index=False)}
ind_to_bp = {v: k for k, v in bp_to_ind.items()}
bp_opts = [ind_to_bp[ind] for ind in range(len(ind_to_bp))]

# Define a dictionary which links amino acids to indices and the reverse
aa_to_ind = {aa: ind for aa, ind in AAtoInd_df.itertuples(index=False)}
ind_to_aa = {v: k for k, v in aa_to_ind.items()}
aa_opts = [ind_to_aa[ind] for ind in range(len(ind_to_aa))]

# Write a dictionary for reverse complements
reverse_comp_dict = {"A": "T",
                     "T": "A",
                     "C": "G",
                     "G": "C",
                     "N": "N"}
