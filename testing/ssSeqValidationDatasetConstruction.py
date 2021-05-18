#!/usr/bin/env python
# coding: utf-8

# Import needed modules
import pandas as pd
import numpy as np
import csv
import pickle
from itertools import chain
from multiprocessing import Pool

# Seed
np.random.seed(27)

# First load the barcode sequences.
barcode_df = pd.read_csv("../ssSeqSupport/IndexMap.csv")

# Define the adapter sequences
f_adapter = "CACCCAAGACCACTCTCCGG"
r_adapter = "GGTAGACGGAGACAGGCGG"

# Define a codon table
CodonTable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
               'GGG': 'G', 'TAG': '*', 'TAA': '*', 'TGA': '*', "TT": "Del",
             "TTTT": "Ins"} 

LimitedCodonTable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 
               'TGT': 'C', 'TGC': 'C'} 

# Define the available codons
available_codons = list(CodonTable.keys())
limited_available_codons = list(LimitedCodonTable.keys())

# Define a reverse complement dictionary
ReverseCompDict = {"A": "T",
                     "T": "A",
                     "C": "G",
                     "G": "C",
                     "N": "N"}

# Define all basepairs
all_bp = ("A", "T", "C", "G")

# Define a length of codon associated with each amino acid
aa_to_length = {aa: len(codon) for codon, aa in CodonTable.items()}
aa_to_length.update({"Ins": 4, "Del": 2})

# Write a function that returns the reverse complement of a sequence
def ReverseComplement(seq):

    # Loop through the sequence in reverse and translate
    return "".join(ReverseCompDict[char] for char in reversed(seq))

# Quality options
quality_opts = [">", "?", "A", "B", "C", "D", "E"]
qual_to_score = {">": 29, "?": 30, "A": 31,
                 "B": 32, "C": 33, "D": 34, "E": 35}

# Write a function for adding junk DNA
def add_junk(dna_string, junk_level, i, uid_ind):
    
    # Make a set of random DNA molecules of length 150
    junk_choices = np.random.choice(available_codons, size = (junk_level, 50), replace = True)
    
    # Add to the dna string
    for junk_choice in junk_choices:
        
        # Add to the uid ind
        uid_ind += 1
        
        # Make a uid
        uid = "@M06418:33:000000000-CRL6Y:1:1101:{}:{} 1:N:0:162".format(i, uid_ind)
        
        # Make a quality score at random
        Q_choices = "".join(["".join(np.random.choice(quality_opts, size = len(codon_choice), replace = True))
                             for codon_choice in junk_choice])
        
        # Join the chosen dna
        complete_junk = "".join(junk_choice)
        
        # Add to the dna string
        dna_string += "\n" + uid
        dna_string += "\n" + complete_junk
        dna_string += "\n+"
        dna_string += "\n" + Q_choices 
        
    # Return the dna_string
    return dna_string

# Write a function that builds reference sequence files
def BuildRefSeqFile(ref_seq, pos_list, savename, n_plates = 8):
    
    # Build the reference sequence
    ref_seq_temp = (ref_seq[:pos_list[0]] + "NNN" + ref_seq[pos_list[0] + 3: pos_list[1]] +
                     "NNN" + ref_seq[pos_list[1] + 3 : pos_list[2]] +
                     "NNN" + ref_seq[pos_list[2] + 3 : pos_list[3]] + 
                     "NNN" + ref_seq[pos_list[3] + 3:])
    
    # Build the reference sequence file
    ref_seq_headers = [["PlateName", "IndexPlate", "ReferenceSequence"]]
    ref_seq_list = [["Test-DI0{}".format(ind+1), "DI0{}".format(ind+1), ref_seq_temp] 
                    for ind in range(n_plates)]
    complete_ref_seq_file = ref_seq_headers + ref_seq_list
    with open(savename.format(i+1), "w") as f:
        writer = csv.writer(f)
        writer.writerows(complete_ref_seq_file)
        
    # Return the reference sequence
    return ref_seq_temp
        
# Write a function that builds a mutated sequence
def BuildMutSeq(ref_seq, pos_list, n_seq_max, available_codons, CodonTable, fbc, rbc):
    
    # Select 4 codons at random
    codon_choices = np.random.choice(available_codons, size = 4, replace = False)
    aa_choices = [CodonTable[codon] for codon in codon_choices]

    # Construct all gene sequences
    new_seq = (ref_seq[:pos_list[0]] + codon_choices[0] + ref_seq[pos_list[0] + 3: pos_list[1]] +
                 codon_choices[1] + ref_seq[pos_list[1] + 3 : pos_list[2]] +
                 codon_choices[2] + ref_seq[pos_list[2] + 3 : pos_list[3]] + 
                 codon_choices[3] + ref_seq[pos_list[3] + 3:])

    # Generate actual fake read sequences
    f_seq = fbc + f_adapter + new_seq[:123]
    r_seq = rbc + r_adapter + ReverseComplement(new_seq[-124:])

    # Generate a number at random which dictates the number of sequences we will see for the row
    n_seqs = np.random.randint(0, n_seq_max)
    
    # Return relevant info
    return f_seq, r_seq, n_seqs, aa_choices

# Write a function that generates quality scores
def GenerateQualScores(aa1, aa2, aa3, aa4, baseline_quality, pos_list):
    
    # Select 4 sets of quality scores at random
    qual_choices = [np.random.choice(quality_opts, size = aa_to_length[aa], replace = True)
                    for aa in [aa1, aa2, aa3, aa4]]
    qual_strings = ["".join(block) for block in qual_choices]

    # Calculate whether all characters give a quality above 30/equal to 30 or not
    greater30 = [None, None, None, None]
    for l, string in enumerate(qual_strings):
        greater30[l] = all([qual_to_score[character]>=30 for character in string])

    # Generate quality scores 
    quality_scores = (baseline_quality[:pos_list[0]] + qual_strings[0] + baseline_quality[pos_list[0] + 3: pos_list[1]] +
                 qual_strings[1] + baseline_quality[pos_list[1] + 3 : pos_list[2]] +
                 qual_strings[2] + baseline_quality[pos_list[2] + 3 : pos_list[3]] + 
                 qual_strings[3] + baseline_quality[pos_list[3] + 3:])

    # Generate fake quality scores
    f_qual = "@"*27 + quality_scores[:123]
    r_qual = "@"*26 + quality_scores[-124:][::-1]
    
    # Return the forward and reverse qualities as well as greater30
    return f_qual, r_qual, greater30

# Write a function that adds to fastq files
def AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter,
               rbc, r_adapter, i, j, f_qual, r_qual):
    
    # Make a copy of our sequences
    temp_fseq = fseq
    temp_rseq = rseq

    # Set noise_f and noise_r to false
    noise_f = False
    noise_r = False

    # At random, replace fseq and rseq with noise. Do this at a low frequency (1%)
    if np.random.random() < 0.01:
        temp_fseq = fbc + f_adapter + "".join(np.random.choice(all_bp, size = 123))
        noise_f = True
        
    if np.random.random() < 0.01:
        temp_rseq = rbc + r_adapter + "".join(np.random.choice(all_bp, size = 124))
        noise_r = True

    # Build the fastq entries
    uid = "@M06418:33:000000000-CRL6Y:1:1101:{}:{} 1:N:0:162".format(i, j)

    # Add to the forward entry
    if i==0 and j==0:
        fstring += uid
    else:
        fstring += "\n" + uid
    fstring += "\n" + temp_fseq
    fstring += "\n+"
    fstring += "\n" + f_qual

    # Add to the reverse entry
    if i==0 and j==0:
        rstring += uid
    else:
        rstring += "\n" + uid
    rstring += "\n" + temp_rseq
    rstring += "\n+"
    rstring += "\n" + r_qual

    # Return fstring and rstring
    return fstring, rstring, noise_f, noise_r


################################## Function for Different Positions ####################################

# Define the generic reference sequence
ref_seq = "ATGGCGCCGACCCTGTCGGAACAGACCCGTCAGTTGGTACGTGCGTCTGTGCCTGCACTGCAGAAACACTCAGTCGCTATTAGCGCCACGATGTATCGGCTGCTTTTCGAACGGTATCCCGAAACGCGGAGCTTGTTTGAACTTCCTGAGAGACAGATACACAAGCTTGCGTCGGCCCTGTTGGCCTACGCCCGTAGTATCGACAACCCATCGGCGTTACAGGCGGCCATCCGCCGCATGGTGCTTTCCC"

# Define the positions of "NNN" in each of the desired reference sequences
pos1 = [9, 57, 111, 168]
pos2 = [56, 110, 167, 200]
pos3 = [109, 166, 199, 232]
position_lists = [pos1, pos2, pos3]
def BuildDifPosSeqs():

    # Build the additional reference sequences
    for i, pos_list in enumerate(position_lists):

        # Build the reference sequence for each position
        _ = BuildRefSeqFile(ref_seq, pos_list, "TestDatasets/DifferentPositions/Pos{}/RefSeqs.csv".format(i+1))

    # Get the length of the quality refseq
    refseq_len = len(ref_seq)
    
    # Loop over each position set
    bc_info_by_pos = [None, None, None,]
    for j, pos_list in enumerate(position_lists):

        # Now loop over each barcode row
        updated_bc_info = [None for _ in range(len(barcode_df))]
        for k, (_, row) in enumerate(barcode_df.iterrows()):

            # Get fbc and rbc
            fbc = row["F-BC"]
            rbc = row["R-BC"]
            
            # Build the mutated sequences
            f_seq, r_seq, n_seqs, aa_choices = BuildMutSeq(ref_seq, pos_list, 101, available_codons, CodonTable,
                                                          fbc, rbc)

            # Record all information
            updated_bc_info[k] = (row.IndexPlate, row.Well, fbc, rbc, n_seqs, 
                                  f_seq, r_seq, *aa_choices)

        # Store updated_bc_info in the prepared list
        bc_info_by_pos[j] = updated_bc_info

    # Loop over the different position info files
    for k, updated_bc_info in enumerate(bc_info_by_pos):

        # Define the pos_list
        pos_list = position_lists[k]

        # Define the forward and reverse fastq file strings
        fstring = ""
        rstring = ""

        # Define a list for holding information about each read
        read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "AA2", "AA3", "AA4", "Q1", "Q2", "Q3", "Q4"]]

        # Build actual fastq files
        for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1, aa2, aa3, aa4) in enumerate(updated_bc_info):

            # Loop over the number of sequences
            for j in range(nseqs):
                
                # Make a baseline quality score
                baseline_quality = "".join(np.random.choice(quality_opts, size = refseq_len))
                
                # Make Qscores
                f_qual, r_qual, greater30 = GenerateQualScores(aa1, aa2, aa3, aa4, baseline_quality, pos_list)

                # Add to the forward and reverse fastq file strings
                fstring, rstring, noise_f, noise_r = AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter, 
                                                                rbc, r_adapter, i, j, f_qual, r_qual)

                # Append to the read list
                read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, aa2, aa3, aa4, *greater30])

            # Add junk to forward and reverse
            fstring = add_junk(fstring, 10, i, j)
            rstring = add_junk(rstring, 10, i, j)

        # Save everything
        with open("TestDatasets/DifferentPositions/Pos{}/TestData_R1_test.fastq".format(k+1), "w") as f:
            f.write(fstring)
        with open("TestDatasets/DifferentPositions/Pos{}/TestData_R2_test.fastq".format(k+1), "w") as f:
            f.write(rstring)
        with open("TestDatasets/DifferentPositions/Pos{}/SequenceInfo.csv".format(k+1), "w") as f:
            writer = csv.writer(f)
            writer.writerows(read_list)


###########################################################################################################################
# Create a dataset where there is overlap between forward and reverse reads. Vary the degree of overlap and see how we do.
############################################################################################################################

# Define the generic reference sequence
overlap_ref_seq = "ATGGCGCCGACCCTGTCGGAACAGACCCGTCAGTTGGTACGTGCGTCTGTGCCTGCACTGCAGAAACACTCAGTCGCTATTAGCGCCACGATGTATCGGCTGCTTTTCGAACGGTATCCCGAAACGCGGAGCTTGTTTGAACTTCCTGAGAGACAGATACACAAGCTTGCGTCGGCCCTGTTGGCCTACGCCCGTAGTAT"

# Define the positions of "NNN" in each of the desired reference sequences
overlap_set1 = [3, 9, 111, 168] # One overlapping position
overlap_set2 = [9, 99, 111, 168] # Two overlapping positions
overlap_set3 = [81, 90, 111, 168] # Three overlapping positions
overlap_set4 = [81, 90, 111, 117] # Four overlapping positions

overlap_position_lists = [overlap_set1, overlap_set2, overlap_set3, overlap_set4]

# Get the length of the overlap_ref_seq
overlap_ref_len = len(overlap_ref_seq)

# Build the additional reference sequences
overlap_refs = [None, None, None, None]
for i, pos_list in enumerate(overlap_position_lists):
    overlap_refs[i] = BuildRefSeqFile(overlap_ref_seq, pos_list, "TestDatasets/DifferentOverlapDegrees/Set{}/RefSeqs.csv".format(i+1))

def BuildDifOverlapSeqs():

    # Loop over each position set
    bc_info_by_pos = [None, None, None, None]
    for j, pos_list in enumerate(overlap_position_lists):

        # Now loop over each barcode row
        updated_bc_info = [None for _ in range(len(barcode_df))]
        for k, (_, row) in enumerate(barcode_df.iterrows()):
            
            # Get barcodes
            fbc = row["F-BC"]
            rbc = row["R-BC"]
            
            # Build the mutated sequences
            f_seq, r_seq, n_seqs, aa_choices = BuildMutSeq(overlap_ref_seq, pos_list, 101, available_codons, CodonTable,
                                                          fbc, rbc)

            # Record all information
            updated_bc_info[k] = (row.IndexPlate, row.Well, fbc, rbc, n_seqs, 
                                  f_seq, r_seq, *aa_choices)

        # Store updated_bc_info in the prepared list
        bc_info_by_pos[j] = updated_bc_info

    # Loop over the different position info files
    for k, updated_bc_info in enumerate(bc_info_by_pos):

        # Define the pos_list
        pos_list = overlap_position_lists[k]

        # Define the forward and reverse fastq file strings
        fstring = ""
        rstring = ""

        # Define a list for holding information about each read
        read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "AA2", "AA3", "AA4", "Q1", "Q2", "Q3", "Q4", "Q1_2", "Q2_2", "Q3_2", "Q4_2"]]

        # Build actual fastq files
        for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1, aa2, aa3, aa4) in enumerate(updated_bc_info):

            # Loop over the number of sequences
            for j in range(nseqs):

                # Make a baseline quality score
                overlap_baseline_quality = "".join(np.random.choice(quality_opts, size = overlap_ref_len))
                
                # Make Qscores
                f_qual, _, greater30_1 = GenerateQualScores(aa1, aa2, aa3, aa4, overlap_baseline_quality, pos_list)
                _, r_qual, greater30_2 = GenerateQualScores(aa1, aa2, aa3, aa4, overlap_baseline_quality, pos_list)

                # Add to the forward and reverse fastq file strings
                fstring, rstring, noise_f, noise_r = AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter, 
                                                                rbc, r_adapter, i, j, f_qual, r_qual)

                # Append to the read list
                read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, aa2, aa3, aa4, *greater30_1, *greater30_2])

            # Add junk to forward and reverse
            fstring = add_junk(fstring, 100, i, j)
            rstring = add_junk(rstring, 100, i, j)

        # Save everything
        with open("TestDatasets/DifferentOverlapDegrees/Set{}/TestData_R1_test.fastq".format(k+1), "w") as f:
            f.write(fstring)
        with open("TestDatasets/DifferentOverlapDegrees/Set{}/TestData_R2_test.fastq".format(k+1), "w") as f:
            f.write(rstring)
        with open("TestDatasets/DifferentOverlapDegrees/Set{}/SequenceInfo.csv".format(k+1), "w") as f:
            writer = csv.writer(f)
            writer.writerows(read_list)


########################################## Different RefSeq by Well ##################################

def BuildDifRefSeq():

    # Define a dictionary to keep the information on each reference sequence chosen
    ref_seq_info = (
        (
            {"FAAs": ("AA1", "AA2", "AA3"),
             "FQs": ("Q1", "Q2", "Q3"),
             "RAAs": ("AA3", "AA4"),
             "RQs": ("Q3", "Q4"),
             "OverlapAAs": ("AA3",)}
        ),
        (
            {"FAAs": ("AA1", "AA2", "AA3"),
             "FQs": ("Q1", "Q2", "Q3"),
             "RAAs": ("AA2", "AA3", "AA4"),
             "RQs": ("Q2", "Q3", "Q4"),
             "OverlapAAs": ("AA2", "AA3")}
        ),
        (
            {"FAAs": ("AA1", "AA2", "AA3"),
             "FQs": ("Q1", "Q2", "Q3"),
             "RAAs": ("AA1", "AA2", "AA3", "AA4"),
             "RQs": ("Q1", "Q2", "Q3", "Q4"),
             "OverlapAAs": ("AA1", "AA2", "AA3")}
        ),
        (
            {"FAAs": ("AA1", "AA2", "AA3", "AA4"),
             "FQs": ("Q1", "Q2", "Q3", "Q4"),
             "RAAs": ("AA1", "AA2", "AA3", "AA4"),
             "RQs": ("Q1", "Q2", "Q3", "Q4"),
             "OverlapAAs": ("AA1", "AA2", "AA3", "AA4")}
        )
    )

    # Create a dictionary to store the reference sequence information for each 
    # well
    platewell_refseq_info = {}

    # Choose a random position set
    # Now loop over each barcode row
    updated_bc_info = [None for _ in range(len(barcode_df))]
    pos_index = []
    for k, (_, row) in enumerate(barcode_df.iterrows()):

        # Choose a random position list
        rand_pos_ind = np.random.choice(4)
        pos_list = overlap_position_lists[rand_pos_ind]
        pos_index.append(rand_pos_ind)
        platewell_refseq_info[row["IndexPlate"] + row["Well"]] = ref_seq_info[rand_pos_ind]

        # Get barcodes
        fbc = row["F-BC"]
        rbc = row["R-BC"]
        
        # Build the mutated sequences
        f_seq, r_seq, n_seqs, aa_choices = BuildMutSeq(overlap_ref_seq, pos_list, 101, available_codons, 
                                                       CodonTable, fbc, rbc)

        # Record all information
        updated_bc_info[k] = (row.IndexPlate, row.Well, fbc, rbc, n_seqs, 
                              f_seq, r_seq, *aa_choices)

    # Save platewell_refseq_info
    with open("./TestDatasets/DifferentRefSeqByWell/RefSeqInfoDict.pkl", "wb") as f:
        pickle.dump(platewell_refseq_info, f)

    # Define the forward and reverse fastq file strings
    fstring = ""
    rstring = ""

    # Define a list for holding information about each read
    read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "AA2", "AA3", "AA4", "Q1", "Q2", "Q3", "Q4"]]
    ref_seq_table = [["PlateName", "IndexPlate", "Well", "ReferenceSequence"]]

    # Build actual fastq files
    for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1, aa2, aa3, aa4) in enumerate(updated_bc_info):

        # Define the position list
        chosen_pos_ind = pos_index[i]
        pos_list = overlap_position_lists[chosen_pos_ind]

        # Grow the ref_seq_table
        ref_seq_table.append(("Test-{}".format(index_plate), index_plate, well, overlap_refs[chosen_pos_ind]))

        # Loop over the number of sequences
        for j in range(nseqs):

            # Make a baseline quality score
            overlap_baseline_quality = "".join(np.random.choice(quality_opts, size = overlap_ref_len))
            
            # Make Qscores
            f_qual, r_qual, greater30 = GenerateQualScores(aa1, aa2, aa3, aa4, overlap_baseline_quality, pos_list)

            # Add to the forward and reverse fastq file strings
            fstring, rstring, noise_f, noise_r = AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter, 
                                                            rbc, r_adapter, i, j, f_qual, r_qual)

            # Append to the read list
            read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, aa2, aa3, aa4, *greater30])

        # Add junk to forward and reverse
        fstring = add_junk(fstring, 10, i, j)
        rstring = add_junk(rstring, 10, i, j)

    # Save everything
    with open("TestDatasets/DifferentRefSeqByWell/TestData_R1_test.fastq", "w") as f:
        f.write(fstring)
    with open("TestDatasets/DifferentRefSeqByWell/TestData_R2_test.fastq", "w") as f:
        f.write(rstring)
    with open("TestDatasets/DifferentRefSeqByWell/SequenceInfo.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(read_list)
    with open("TestDatasets/DifferentRefSeqByWell/RefSeqs.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(ref_seq_table)


#################################### Mixed Wells #####################################

def BuildMixedWellsSeq():

    # First load the barcode sequences.
    barcode_df = pd.read_csv("../ssSeqSupport/IndexMap.csv")

    # Define the adapter sequences
    f_adapter = "CACCCAAGACCACTCTCCGG"
    r_adapter = "GGTAGACGGAGACAGGCGG"

    # Choose a position list
    pos_list = overlap_position_lists[2]

    # Build the refseq file
    _ = BuildRefSeqFile(overlap_ref_seq, pos_list, "TestDatasets/MixedWells/RefSeqs.csv")

    # Choose a random position set
    # Now loop over each barcode row
    updated_bc_info = [None for _ in range(len(barcode_df)*2)]
    for k, (_, row) in enumerate(chain(barcode_df.iterrows(), barcode_df.iterrows())):

        # Get barcodes
        fbc = row["F-BC"]
        rbc = row["R-BC"]
        
         # Build the mutated sequences
        f_seq, r_seq, n_seqs, aa_choices = BuildMutSeq(overlap_ref_seq, pos_list, 101, available_codons, CodonTable,
                                                      fbc, rbc)

        # Record all information
        updated_bc_info[k] = (row.IndexPlate, row.Well, fbc, rbc, n_seqs, 
                              f_seq, r_seq, *aa_choices)

    #Seed
    np.random.seed(5)

    # Define the forward and reverse fastq file strings
    fstring = ""
    rstring = ""

    # Define a list for holding information about each read
    read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "AA2", "AA3", "AA4", "Q1", "Q2", "Q3", "Q4"]]

    # Build actual fastq files
    for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1, aa2, aa3, aa4) in enumerate(updated_bc_info):

        # Loop over the number of sequences
        for j in range(nseqs):

            # Make a baseline quality score
            overlap_baseline_quality = "".join(np.random.choice(quality_opts, size = overlap_ref_len))
            
            # Make Qscores
            f_qual, r_qual, greater30 = GenerateQualScores(aa1, aa2, aa3, aa4, overlap_baseline_quality, pos_list)

            # Add to the forward and reverse fastq file strings
            fstring, rstring, noise_f, noise_r = AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter, 
                                                            rbc, r_adapter, i, j, f_qual, r_qual)

            # Append to the read list
            read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, aa2, aa3, aa4, *greater30])

        # Add junk to forward and reverse
        fstring = add_junk(fstring, 10, i, j)
        rstring = add_junk(rstring, 10, i, j)

    # Save everything
    with open("TestDatasets/MixedWells/TestData_R1_test.fastq", "w") as f:
        f.write(fstring)
    with open("TestDatasets/MixedWells/TestData_R2_test.fastq", "w") as f:
        f.write(rstring)
    with open("TestDatasets/MixedWells/SequenceInfo.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(read_list)


######################### Forward and Reverse Don't Agree #######################

def BuildFRDisagreeSeqs():

    # First load the barcode sequences.
    barcode_df = pd.read_csv("../ssSeqSupport/IndexMap.csv")

    # Define the adapter sequences
    f_adapter = "CACCCAAGACCACTCTCCGG"
    r_adapter = "GGTAGACGGAGACAGGCGG"

    # Loop over each position set
    bc_info_by_pos = [None, None]
    for j, pos_list in enumerate(overlap_position_lists[:2]):

        # Build refseq
        _ = BuildRefSeqFile(overlap_ref_seq, pos_list, "TestDatasets/ForwardAndReverseDisagree/Set{}/RefSeqs.csv".format(j+1))

        # Now loop over each barcode row
        updated_bc_info = [None for _ in range(len(barcode_df))]
        for k, (_, row) in enumerate(barcode_df.iterrows()):

            # Get barcodes
            fbc = row["F-BC"]
            rbc = row["R-BC"]
            
            # Build the mutated sequences
            f_seq, _, n_seqs, aa_choices1 = BuildMutSeq(overlap_ref_seq, pos_list, 101, 
                                                           limited_available_codons, LimitedCodonTable,
                                                       fbc, rbc)
            _, r_seq, _, aa_choices2 = BuildMutSeq(overlap_ref_seq, pos_list, 101, 
                                                           limited_available_codons, LimitedCodonTable,
                                                  fbc, rbc)

            # Record all information
            updated_bc_info[k] = (row.IndexPlate, row.Well, fbc, rbc, n_seqs, 
                                  f_seq, r_seq, *aa_choices1, *aa_choices2)

        # Store updated_bc_info in the prepared list
        bc_info_by_pos[j] = updated_bc_info

    # Loop over the different position info files
    for k, updated_bc_info in enumerate(bc_info_by_pos):

        # Define the pos_list
        pos_list = overlap_position_lists[k]

        # Define the forward and reverse fastq file strings
        fstring = ""
        rstring = ""

        # Define a list for holding information about each read
        read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "AA2", "AA3", "AA4", "AA1_2", "AA2_2", "AA3_2", "AA4_2", "Q1", "Q2", "Q3", "Q4"]]

        # Build actual fastq files
        for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1, aa2, aa3, aa4, aa1_2, aa2_2, aa3_2, aa4_2) in enumerate(updated_bc_info):

            # Loop over the number of sequences
            for j in range(nseqs):

                # Make a baseline quality score
                overlap_baseline_quality = "".join(np.random.choice(quality_opts, size = overlap_ref_len))
                
                 # Make Qscores
                f_qual, r_qual, greater30 = GenerateQualScores(aa1, aa2, aa3, aa4, overlap_baseline_quality, pos_list)

                # Add to the forward and reverse fastq file strings
                fstring, rstring, noise_f, noise_r = AddToFastQ(fseq, rseq, fstring, rstring, fbc, f_adapter, 
                                                                rbc, r_adapter, i, j, f_qual, r_qual)

                # Append to the read list
                read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, aa2, aa3, aa4, aa1_2, aa2_2, aa3_2, aa4_2, *greater30])

            # Add junk to forward and reverse
            fstring = add_junk(fstring, 10, i, j)
            rstring = add_junk(rstring, 10, i, j)

        # Save everything
        with open("TestDatasets/ForwardAndReverseDisagree/Set{}/TestData_R1_test.fastq".format(k+1), "w") as f:
            f.write(fstring)
        with open("TestDatasets/ForwardAndReverseDisagree/Set{}/TestData_R2_test.fastq".format(k+1), "w") as f:
            f.write(rstring)
        with open("TestDatasets/ForwardAndReverseDisagree/Set{}/SequenceInfo.csv".format(k+1), "w") as f:
            writer = csv.writer(f)
            writer.writerows(read_list)


############################################ Single Site ################################

def BuildSingleSiteSeqs():
    # Define the generic reference sequence
    single_ref_seq = "ATGGCGCCGACCCTGTCGGAACAGACCCGTCAGTTGGTACGTGCGTCTGTGCCTGCACTGCAGAAACACTCAGTCGCTATTAGCGCCACGATGTATCGGCTGCTTTTCGAACGGTATCCCGAAACGCGGAGCTTGTTTGAACTTCCTGAG"
    
    # Get the length of the reference sequence
    len_single_ref = len(single_ref_seq)
    
    # Define the positions of "NNN" in each of the desired reference sequences
    mut_pos = 75

    # Build the additional reference sequences
    single_ref_NNN = single_ref_seq[:mut_pos] + "NNN" + single_ref_seq[mut_pos + 3:]

    # Build the reference sequence file
    ref_seq_headers = [["PlateName", "IndexPlate", "ReferenceSequence"]]
    ref_seq_list = [["Test-DI0{}".format(ind+1), "DI0{}".format(ind+1), single_ref_NNN] 
                    for ind in range(8)]
    complete_ref_seq_file = ref_seq_headers + ref_seq_list
    with open("TestDatasets/SingleSite/RefSeqs.csv".format(), "w") as f:
        writer = csv.writer(f)
        writer.writerows(complete_ref_seq_file)

    # First load the barcode sequences.
    barcode_df = pd.read_csv("../ssSeqSupport/IndexMap.csv")

    # Define the adapter sequences
    f_adapter = "CACCCAAGACCACTCTCCGG"
    r_adapter = "GGTAGACGGAGACAGGCGG"

    # Loop over all barcodes
    updated_bc_info = [None for _ in range(len(barcode_df))]
    for k, (_, row) in enumerate(barcode_df.iterrows()):

        # Select 1 codon at random
        codon_choices = np.random.choice(available_codons, size = 1, replace = False)
        aa_choices = [CodonTable[codon] for codon in codon_choices]

        # Construct all gene sequences
        new_seq = (single_ref_seq[:mut_pos] + codon_choices[0] + single_ref_seq[mut_pos + 3:])

        # Generate actual fake read sequences
        f_seq = row["F-BC"] + f_adapter + new_seq[:123]
        r_seq = row["R-BC"] + r_adapter + ReverseComplement(new_seq[-124:])

        # Generate a number at random which dictates the number of sequences we will see for the row
        n_seqs = np.random.randint(0, 101)

        # Record all information
        updated_bc_info[k] = (row.IndexPlate, row.Well, row["F-BC"], row["R-BC"], n_seqs, 
                              f_seq, r_seq, *aa_choices)

    # Define the forward and reverse fastq file strings
    fstring = ""
    rstring = ""

    # Define a list for holding information about each read
    read_list = [["IndexPlate", "Well", "F-BC", "R-BC", "NoiseF", "NoiseR", "AA1", "Q1"]]

    # Build actual fastq files
    for i, (index_plate, well, fbc, rbc, nseqs, fseq, rseq, aa1) in enumerate(updated_bc_info):

        # Loop over the number of sequences
        for j in range(nseqs):

            # Make a copy of our sequences
            temp_fseq = fseq
            temp_rseq = rseq

            # Set noise_f and noise_r to false
            noise_f = False
            noise_r = False

            # At random, replace fseq and rseq with noise. Do this at a low frequency (1%)
            if np.random.random() < 0.01:
                temp_fseq = fbc + f_adapter + "".join(np.random.choice(all_bp, size = 123))
                noise_f = True

            if np.random.random() < 0.01:
                temp_rseq = rbc + r_adapter + "".join(np.random.choice(all_bp, size = 124))
                noise_r = True
            
            # Select 8 sets of quality scores at random
            qual_choices = np.random.choice(quality_opts, size = (1, aa_to_length[aa1]), replace = True)
            qual_strings = ["".join(block) for block in qual_choices]

            # Calculate whether all characters give a quality above 30/equal to 30 or not
            greater30 = [None]
            for l, string in enumerate(qual_strings):
                greater30[l] = all([qual_to_score[character]>=30 for character in string])

            # Make a baseline quality score
            single_baseline_quality = "".join(np.random.choice(quality_opts, size = len_single_ref))
                
            # Generate quality scores 
            single_quality_scores1 = (single_baseline_quality[:mut_pos] + qual_strings[0] + single_baseline_quality[mut_pos + 3:])

            # Generate fake quality scores
            f_qual = "@"*27 + single_quality_scores1[:123]
            r_qual = "@"*26 + single_quality_scores1[-124:][::-1]

            # Build the fastq entries
            uid = "@M06418:33:000000000-CRL6Y:1:1101:{}:{} 1:N:0:162".format(i, j)

            # Add to the forward entry
            if i==0 and j==0:
                fstring += uid
            else:
                fstring += "\n" + uid
            fstring += "\n" + temp_fseq
            fstring += "\n+"
            fstring += "\n" + f_qual

            # Add to the reverse entry
            if i==0 and j==0:
                rstring += uid
            else:
                rstring += "\n" + uid
            rstring += "\n" + temp_rseq
            rstring += "\n+"
            rstring += "\n" + r_qual

            # Append to the read list
            read_list.append([index_plate, well, fbc, rbc, noise_f, noise_r, aa1, *greater30,])

    # Save everything
    with open("TestDatasets/SingleSite/TestData_R1_test.fastq".format(k+1), "w") as f:
        f.write(fstring)
    with open("TestDatasets/SingleSite/TestData_R2_test.fastq".format(k+1), "w") as f:
        f.write(rstring)
    with open("TestDatasets/SingleSite/SequenceInfo.csv".format(k+1), "w") as f:
        writer = csv.writer(f)
        writer.writerows(read_list)

########################### Package and Run Everything #####################################

# Define a dictionary containing datasets to check
task_dict = {"DifPos": BuildDifPosSeqs,
            "DifOverlap": BuildDifOverlapSeqs,
            "DifRefSeq": BuildDifRefSeq,
            "MixedWells": BuildMixedWellsSeq,
            "FRDisagree": BuildFRDisagreeSeqs,
            "SingleSite": BuildSingleSiteSeqs}

# Write a helper function
def SeqBuildHelper(task):
    
    # Run the task
    task_dict[task]()
    
# Write a function to multiprocess generation of datasets
def MultiprocessBuildSeqs(njobs=6):
    
    # Run each task
    tasks = list(task_dict.keys())
    
    # Multiprocess
    with Pool(njobs) as p:
        _ = list(p.imap_unordered(SeqBuildHelper, tasks))
        
# Run buidseqs
if __name__ == "__main__":
    MultiprocessBuildSeqs()