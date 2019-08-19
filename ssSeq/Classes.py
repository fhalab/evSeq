# Import all of the required 3rd party modules
import re
import csv
import os
import os.path
import numpy as np
import pandas as pd
from Bio import pairwise2
from tqdm import tqdm
from multiprocessing import Pool

# Import ssSeq modules
from ssSeq.GlobalSetup import *

# Define a function that analyzes a well in the plate when
# we are in troubleshooting mode
def multiprocess_plate_analyzer_ts(well):

    # Make and process the alignments
    well.align()
    well.analyze_alignments()
    well.generate_consensus()
    well.generate_consensus(forward=False)

    # Generate axis labels for matrices
    f_ref_seq_aas = well.f_ref_seq.aas_as_list
    r_ref_seq_aas = well.r_ref_seq.aas_as_list
    f_ref_seq_bps = well.f_ref_seq.bps_as_list
    r_ref_seq_bps = well.r_ref_seq.bps_as_list

    # Get the variant counts
    variant_counts = well.alignment_results[-1]

    # Get the depth (the sum of all counts in the variant)
    var_depth = sum(variant_counts.values())

    # Create the well-specific variant counts table
    variant_info = [[well.plate, well.well, well.f_barcode, well.r_barcode,
                     combo, count/var_depth, var_depth]
                    for combo, count in variant_counts.items()]

    # Generate the header
    header = "{}_{}_{}_{}".format(well.plate, well.well,
                                  well.f_barcode, well.r_barcode)

    # Create the well-specific summary tables
    summary_info = []
    bp_count_freq_info = []
    aa_count_freq_info = []

    # Analyze both forward and reverse results
    for (label, aa_labels, bp_labels, insertions, deletions, bp_counts,
         aa_counts, bp_freq, aa_freq, var_aa_info) in zip(("Forward", "Reverse"),
                                                       (f_ref_seq_aas, r_ref_seq_aas),
                                                       (f_ref_seq_bps, r_ref_seq_bps),
                                                       *well.alignment_results[:-1]):

        ############### Build frequency and count forms ##############
         # Add an empty well to the amino acid and basepair labels
        aa_labels = [None] + aa_labels
        bp_labels = [None] + bp_labels

        # Convert matrices to lists
        aa_count_list = aa_counts.T.tolist()
        aa_freq_list = aa_freq.T.tolist()
        bp_count_list = bp_counts.T.tolist()
        bp_freq_list = bp_freq.T.tolist()

        # Build the amino acid table
        aa_table = [[header + "_Counts_" + label],
                    aa_labels]
        aa_table.extend([[aa] + sublist for aa, sublist in zip(aa_opts, aa_count_list)])
        aa_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         aa_labels])
        aa_table.extend([[aa] + sublist for aa, sublist in zip(aa_opts, aa_freq_list)])
        aa_table.extend([[None], [None]])

        # Build the basepair table
        bp_table = [[header + "_Counts_" + label],
                    bp_labels]
        bp_table.extend([[bp] + sublist for bp, sublist in zip(bp_opts, bp_count_list)])
        bp_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         bp_labels])
        bp_table.extend([[bp] + sublist for bp, sublist in zip(bp_opts, bp_freq_list)])
        bp_table.extend([[None], [None]])

        # Append the completed tables to the master lists
        aa_count_freq_info.extend(aa_table)
        bp_count_freq_info.extend(bp_table)

        # Append to the summary table
        for site, freq_dict, depth in var_aa_info:
            summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
                                  site + 1, aa, freq, depth] for aa, freq in freq_dict.items()])


    ##################### Format Alignment Results #######################
    # Add forward alignments
    section_break = "######################################################################################\n"
    alignment_text = section_break
    alignment_text = alignment_text + header + "_ForwardAlignments\n\n"
    for i, alignment in enumerate(well.f_alignments, 1):

        # Add alignment
        alignment_text += "Sequence {}\n".format(i)
        alignment_text += pairwise2.format_alignment(*alignment)
        alignment_text += "\n"

    # Add reverse alignments
    alignment_text = alignment_text + "\n\n" + section_break + header + "_ReverseAlignments\n\n"
    for i, alignment in enumerate(well.r_alignments, 1):

        # Add alignment
        alignment_text += "Sequence {}\n".format(i)
        alignment_text += pairwise2.format_alignment(*alignment)
        alignment_text += "\n"

    alignment_text += "\n\n\n"

    # Add consensus text
    consensus_text = header + "_ForwardConsensus\n"
    consensus_text = consensus_text + well.f_consensus + "\n\n"
    consensus_text = consensus_text + header + "_ReverseConsensus\n"
    consensus_text = consensus_text + well.r_consensus + "\n"
    consensus_text += "\n\n\n"

    # Return all information
    return (aa_count_freq_info, bp_count_freq_info, alignment_text,
            consensus_text, summary_info, variant_info)

# Define a function that analyzes a well in the plate when
# we are not in troubleshooting mode
def multiprocess_plate_analyzer(well):

    # Make and process the alignments
    well.align()
    well.analyze_alignments()
    well.generate_consensus()
    well.generate_consensus(forward=False)

    # Generate axis labels for matrices
    f_ref_seq_aas = well.f_ref_seq.aas_as_list
    r_ref_seq_aas = well.r_ref_seq.aas_as_list
    f_ref_seq_bps = well.f_ref_seq.bps_as_list
    r_ref_seq_bps = well.r_ref_seq.bps_as_list

    # Get the variant counts
    variant_counts = well.alignment_results[-1]

    # Get the depth (the sum of all counts in the variant)
    var_depth = sum(variant_counts.values())

    # Create the well-specific variant counts table
    variant_info = [[well.plate, well.well, well.f_barcode, well.r_barcode,
                     combo, count/var_depth, var_depth]
                    for combo, count in variant_counts.items()]

    # Generate the header
    header = "{}_{}_{}_{}".format(well.plate, well.well,
                                  well.f_barcode, well.r_barcode)

    # Create the well-specific summary tables
    summary_info = []

    # Analyze both forward and reverse results
    for (label, aa_labels, bp_labels, insertions, deletions, bp_counts,
         aa_counts, bp_freq, aa_freq, var_aa_info) in zip(("Forward", "Reverse"),
                                                       (f_ref_seq_aas, r_ref_seq_aas),
                                                       (f_ref_seq_bps, r_ref_seq_bps),
                                                       *well.alignment_results[:-1]):
        # Append to the summary table
        for site, freq_dict, depth in var_aa_info:
            summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
                                  site + 1, aa, freq, depth] for aa, freq in freq_dict.items()])

    # Add consensus text
    consensus_text = header + "_ForwardConsensus\n"
    consensus_text = consensus_text + well.f_consensus + "\n\n"
    consensus_text = consensus_text + header + "_ReverseConsensus\n"
    consensus_text = consensus_text + well.r_consensus + "\n"
    consensus_text += "\n\n\n"

    # Return all information
    return consensus_text, summary_info, variant_info

###############################################################################
################################## Classes ####################################
class translation():

    # Initialize
    def __init__(self, seq, start_ind, codon_table, low_quality_chars = []):

        # Get the number of codons in the sequence
        n_codons = np.floor((len(seq) - start_ind)/3)

        # Identify all codons
        codons = [seq[int(start_ind + i*3): int(start_ind + (i+1)*3)] for i in range(int(n_codons))]

        # Link bp to codon
        bp_ranges = [list(range(int(start_ind + i*3), int(start_ind + (i+1)*3)))
                     for i in range(int(n_codons))]
        bp_to_codon = {bp: codon for codon, bp_list in enumerate(bp_ranges)
                       for bp in bp_list}

        # Identify low quality codons
        self.low_quality_codons = {bp_to_codon[bp] for bp in low_quality_chars
                                   if bp in bp_to_codon.keys()}

        # Translate all codons
        translation = []
        for codon in codons:

            # If this is "NNN" we give a question mark
            if "N" in codon:
                translation.append("?")
            else:
                translation.append(codon_table[codon])


        # Store the translation
        self.translation  = "".join(translation)

class seq_pair():

    # Define the initialization. This takes the sequence information and
    # organizes it
    def __init__(self, info_block, start_f = True):

        # Parse the first-seen block to build the object. Set the forward
        # or reverse sequence as appropriate
        # Q = -10log10(e)
        # P = 10^(-Q/10)
        if start_f:
            raw_id, self.f_seq, qual_mark, f_ASCII_qual_scores = info_block
            self.f_qual_scores = [ord(char)-33 for char in f_ASCII_qual_scores]

            # 10^(-Q/10) * 10^(-Q/10) = 10^((Q1 + Q2)/-10)
        else:
            raw_id, self.r_seq, qual_mark, r_ASCII_qual_scores = info_block
            self.r_qual_scores = [ord(char)-33 for char in r_ASCII_qual_scores]

        # Get the id information for the block
        (instrument_name, runid, flowcellid, lane,
         tile, x_coord, y_coord,
         pair_number, filtered_bin,
          control_bin, sample_number) = get_block_info(raw_id)

        # Create a unique id for the pair that can pair ends
        self.id = (instrument_name, lane, tile, x_coord,
                   y_coord, sample_number)

    # Write a function to attach the partner information
    def attach_partner(self, info_block, r_partner = True):

        # Depending on the partner that's being appended, attach different values
        if r_partner:
            raw_id, self.r_seq, qual_mark, r_ASCII_qual_scores = info_block
            self.r_qual_scores = [ord(char)-33 for char in r_ASCII_qual_scores]
        else:
            raw_id, self.f_seq, qual_mark, f_ASCII_qual_scores = info_block
            self.f_qual_scores = [ord(char)-33 for char in f_ASCII_qual_scores]

        # Force everything to be uppercase
        self.f_seq = self.f_seq.upper()
        self.r_seq = self.r_seq.upper()

        # Identify the forward barcode and its prob scores
        self.f_barcode = self.f_seq[:f_bc_length]

        # Identify the reverse barcode and its prob scores
        self.r_barcode = self.r_seq[:r_bc_length]

        # Define sequences without the barcodes
        self.barcodeless_f = self.f_seq[f_bc_length:]
        self.barcodeless_f_q = self.f_qual_scores[f_bc_length:]
        self.reversed_barcodeless_r = reverse_complement(self.r_seq[r_bc_length:])
        self.reversed_barcodeless_r_q = list(reversed(self.r_qual_scores[r_bc_length:]))

    # Write a function to determine if this is a sequence without a partner
    def is_orphan(self):

        # Test to see if we have forward and reverse sequences
        if hasattr(self, "f_seq") and hasattr(self, "r_seq"):
            return False
        else:
            return True

# Write a class which holds all sequence blocks of a given set of barcodes
class well():

    """
    Desired functionality:
    1. Identify the appropriate plate and well for the sequence grouping
    2. Align all sequences to the appropriate reference sequence(s)
    3. Determine probabilities of each base in the alignment
    4. Identify the amino acid(s) at the variable position and assign a probability
        - One option for this is just to pass in reference sequences where the
          variable region(s) is denoted as "NNN"
    5. Identify any other mutations in the alignment. Return a probability score.
    """

    # Define the initialization function
    def __init__(self, seq_pairs):

        # Set the group of sequences as a class attribute
        self.seq_pairs = seq_pairs

        # Identify the forward and reverse barcodes for the object
        self.f_barcode = seq_pairs[0].f_barcode
        self.r_barcode = seq_pairs[0].r_barcode

        # Confirm that the barcode pair exists in the barcode dictionary
        bc_pair = (self.f_barcode, self.r_barcode)
        if bc_pair in BCs_to_refseq.keys():

            # Identify the reference sequence for the barcode pair
            (self.f_ref_seq, self.r_ref_seq,
             self.plate_name, self.n_variable_sites) = BCs_to_refseq[bc_pair]

            # Identify the plate and well for the barcode pair
            self.plate, self.well = BC_to_well[(self.f_barcode, self.r_barcode)]

            # Confirm that this is a real well
            self.real_well = True

        # If the barcode pair does not exist, then this is not a real well
        # It likely results from a sequencing error
        else:
            self.real_well = False

        # Get the number of sequences paired to the well
        self.n_seqs = len(seq_pairs)

    # Define a function for aligning each sequence to the appropriate reference
    def align(self):

        # Loop over the sequence pairs and make alignments to the reference sequences
        self.f_alignments = [pairwise2.align.globalxs(self.f_ref_seq.seq,
                                                             pair.barcodeless_f,
                                                             -2, -1, one_alignment_only=True,
                                                             penalize_end_gaps = False)[0]
                             for pair in self.seq_pairs]
        self.r_alignments = [pairwise2.align.globalxs(self.r_ref_seq.seq,
                                                             pair.reversed_barcodeless_r,
                                                             -2, -1, one_alignment_only=True,
                                                             penalize_end_gaps = False)[0]
                             for pair in self.seq_pairs]

    # Write a helper function that analyzes a single alignment
    def analyze_alignment(self, alignment, j, f_r_ind, bp_counts):

        # Unpack the alignment
        (reference_alignment, read_alignment, score, begin, end) = alignment

        # If the alignment score is below the threshold, return this information
        # and move on
        if f_r_ind==0:
            if score/self.f_ref_seq.seq_length < alignment_cutoff:
                return None
        else:
            if score/self.r_ref_seq.seq_length < alignment_cutoff:
                return None

        # Define the deletions and , low_quality_charsinsertions list
        insertions = []
        deletions = []

        # Define lists for collecting information
        holding_dels = []
        var_sites = []

        # Define a list for holding low quality characters
        low_quality_chars = []

        # Set the reference sequence
        if f_r_ind == 0:
            ref_seq = self.f_ref_seq
        else:
            ref_seq = self.r_ref_seq

        # Define variables that will be used
        ref_counter = -1
        read_counter = -1
        del_possible = False
        bad_alignment = False
        ins_found = False
        del_found = False
        in_ref = False

        # Loop over each base in the alignments
        for i, (ref_char, read_char) in enumerate(zip(reference_alignment, read_alignment)):

            # If the read counter is not a space, add 1.
            if read_char != "-":
                read_counter += 1

            # Until we have a non-space character, we do not record information.
            # We only increment the reference counter when we have non-space characters
            if ref_char != "-":
                in_ref = True
                ref_counter += 1

            # Our read should always start before our reference (as we need to
            # read the adapter region). If we are in reference before we are at
            # a read character, this is a bad alignment. Return None.
            if in_ref and read_counter==-1:
                return None

            # If we are not in reference, just keep going until we are
            if not in_ref:
                continue

            # If the reference counter is equal to the translation start point,
            # then record the aligned start point translating the read
            if ref_counter == ref_seq.trans_start:
                read_trans_start = read_counter

            # Once we are past the length of the reference we can terminate
            if ref_counter==ref_seq.seq_length - 1:
                bp_counts[ref_counter, bp_to_ind[read_char]] += 1
                break

            # Determine if the ref char is a space. If it is, this indicates an
            # insertion in the read sequence.
            if ref_char=="-" and in_ref and ref_counter < ref_seq.seq_length - 1:
                insertions.append([j, ref_counter, read_char])

                # Report that an insertion was found
                ins_found = True

                # Continue to the next character
                continue

            # Determine if the read character surpasses our quality score filter.
            # If it does not, record and continue
            if f_r_ind == 0:
                if self.seq_pairs[j].barcodeless_f_q[read_counter] < q_cutoff:
                    low_quality_chars.append(read_counter)
                    continue
            else:
                if self.seq_pairs[j].reversed_barcodeless_r_q[read_counter] < q_cutoff:
                    low_quality_chars.append(read_counter)
                    continue

            # Determine if the read char is a space. If it is, this indicates a deletion.
            # It is okay if we have a space if we are at the end of the read. Sometimes
            # reads are shorter than the reference.
            if read_char=="-":

                # Mark that we might have a deletion. We need to wait for a non-space
                # base to show up to be sure that it is a deletion
                del_possible = True
                holding_dels.append([j, ref_counter, ref_char])

            # Determine if the read character is a non-space and we have deletions
            # in holding. If so, append deletions to the deletions list.
            if read_char!="-" and del_possible:

                # Append deletions in holding
                deletions.extend(holding_dels)

                # Report that a deletion was found
                del_found = True

                # Reset relevant values
                del_possible = False
                holding_dels = []

            # Add to the alignment matrix to determine what base is called at
            # what position
            bp_counts[ref_counter, bp_to_ind[read_char]] += 1

        # Return relevant values
        return insertions, deletions, read_trans_start, ins_found, del_found, low_quality_chars

    # Define a function for performing data analysis on the alignments
    def analyze_alignments(self):

        """
        Notes:
        1. Assume that NNN denotes the frame.
        2. Look for frameshifts upstream and downstream of NNN.
        3. How do you report frameshifts to the user?
        4. How do you report insertions/deletions to the user?
        """

        ######################## Upfront Processing #############################
        # Determine the number of alignments in the set
        n_f_alignments = len(self.f_alignments)
        n_r_alignments = len(self.r_alignments)
        assert n_f_alignments == n_r_alignments, "Unequal number of alignments."

        # Build a matrix to count the occurences of each base type at each position
        # in the reference sequence
        bp_counts = [np.zeros([self.f_ref_seq.seq_length, 6]),
                     np.zeros([self.r_ref_seq.seq_length, 6])]
        aa_counts = [np.zeros([self.f_ref_seq.trans_length, 22]),
                     np.zeros([self.r_ref_seq.trans_length, 22])]

        # Package redundant arguments
        ref_seqs = [self.f_ref_seq, self.r_ref_seq]
        n_alignments = [n_f_alignments, n_r_alignments]

        # Identify the number of variable positions in each reference
        f_read_var_sites = len(self.f_ref_seq.var_aa_sites)
        r_read_var_sites = len(self.r_ref_seq.var_aa_sites)
        both_read_var_sites = [f_read_var_sites, r_read_var_sites]

        # Get the maximum number of sites expected to vary in each read
        highest_refs = [max(self.f_ref_seq.var_aa_sites),
                        max(self.r_ref_seq.var_aa_sites)]

        # Get the number of overlapping sites between the two
        n_overlapping_sites = f_read_var_sites + r_read_var_sites - self.n_variable_sites

        # Loop over each alignment in the set and analyze
        insertions = [[],[]]
        deletions = [[],[]]
        variant_counts = {}
        variants = []
        for j, (f_alignment, r_alignment) in enumerate(zip(self.f_alignments, self.r_alignments)):

            # Pull the unaligned read sequences
            unaligned_read_seqs = [self.seq_pairs[j].barcodeless_f,
                                   self.seq_pairs[j].reversed_barcodeless_r]

            # Create a series of tests
            poor_alignment = False
            check_combo = True

            # Analyze the forward and reverse alignments
            for i, alignment in enumerate((f_alignment, r_alignment)):

                # Analyze the alignment
                alignment_output = self.analyze_alignment(alignment,
                                                            j, i, bp_counts[i])

                # If the alignment did not meet the alignment score cutoff, continue
                if alignment_output is None:
                    poor_alignment = True
                    continue

                # Otherwise, unpack the alignment results
                (temp_insertions,
                 temp_deletions,
                 read_trans_start,
                 ins_found,
                 del_found,
                 low_quality_chars) = alignment_output

                # Record insertions and deletions
                insertions[i].extend(temp_insertions)
                deletions[i].extend(temp_deletions)

                # Don't bother moving on if there was an insertion or deletion
                # Don't bother moving on if the alignment quality was poor.
                # Don't check the combination if this is the case for either
                # alignment in the forward-reverse pair.
                if ins_found or del_found or poor_alignment:
                    check_combo = False
                    poor_alignment = False
                    continue

                # Translate each codon in the read alignment
                alignment_translation = translation(unaligned_read_seqs[i],
                                                    read_trans_start,
                                                    codon_table,
                                                    low_quality_chars)

                # Point to the appropriate reference sequence
                ref_seq = ref_seqs[i]

                # Add to the amino acid count matrix and compare the translation to the
                # reference
                bad_codon_found = False
                for k, aa_read in enumerate(alignment_translation.translation[:ref_seq.trans_length]):

                    # Add to the amino acid count matrix if it is not a low quality
                    # codon. If this codon is a variable position, continue adding
                    # to the alignment matrix, but don't try and determine the
                    # variant. Just continue in the parent loop from here
                    if k in alignment_translation.low_quality_codons:
                        if k in ref_seq.var_aa_sites:
                            bad_codon_found = True
                        continue
                    aa_counts[i][k, aa_to_ind[aa_read]] += 1

                # If either of the reads failed for any reason, break the loop.
                # We can't make combination data without both reads.
                if bad_codon_found or not check_combo:
                    break

                # Identify the highest possible reference and the number of variable
                # sites in the read
                highest_ref = highest_refs[i]
                read_var_sites = both_read_var_sites[i]

                # Go this route if this is the first iteration (we are looking at the
                # forward read)
                if i==0:

                    # Create a string to store the variant combination
                    variant_combo = ""

                    # Identify the combination of mutations at the variable position
                    # Skip if the read isn't long enough to reach all of them
                    counter = 0
                    site = ref_seq.var_aa_sites[counter]
                    while site < len(alignment_translation.translation):
                        variant_combo += alignment_translation.translation[site]
                        counter += 1

                        # Break the loop if we have reached the max
                        if counter == read_var_sites:
                            break

                        # Set the next site
                        site = ref_seq.var_aa_sites[counter]

                    # Add "-" on to the end of the variant for all missing sites
                    missing_sites = read_var_sites - counter
                    for _ in range(missing_sites):
                        variant_combo += "-"

                # If we are working with the reverse sequence, make sure our overlaps
                # match
                else:

                    # If we had any overlapping sites, identify what we had in the
                    # forward read
                    if n_overlapping_sites > 0:

                        # Pull from the end of the variant combo the number of overlaps
                        overlaps = variant_combo[-n_overlapping_sites:]

                    # Identify the combination of mutations at the variable position
                    # Skip if the read isn't long enough to reach all of them
                    counter = 0
                    site = ref_seq.var_aa_sites[counter]
                    mismatch = False
                    while site < len(alignment_translation.translation):

                        # Identify the variable position
                        var_pos = alignment_translation.translation[site]

                        # If the counter is less than the number of overlapping sites,
                        # then check to make sure the forward and reverse reads agree.
                        # if they do not agree, check to see if the mismatch is because
                        # the forward read was too short. If the forward read was too
                        # short, rescue by using the reverse call. If the forward read was long
                        # enough and we have a mismatch, set the mismatch flag to True
                        # and break
                        if (counter < n_overlapping_sites and
                            var_pos != overlaps[counter]):

                            # Check to see if this is because the forward read was
                            # too short
                            if overlaps[counter] == "-":

                                # Rescue the variant combo, up the counter,
                                # set the next site, continue
                                variant_combo = variant_combo[:-n_overlapping_sites] + var_pos
                                counter += 1

                                # Break the loop if we have reached the max
                                if counter == read_var_sites:
                                    break

                                # Set the next site
                                site = ref_seq.var_aa_sites[counter]
                                continue

                            # Otherwise, set the mismatch flag to True and break
                            else:
                                mismatch = True
                                break

                        # If the counter is less than the number of overlapping
                        # sites, add to the counter, set the next site and continue
                        elif counter < n_overlapping_sites:
                            counter += 1

                            # Break the loop if we have reached the max
                            if counter == read_var_sites:
                                break

                            # Set the next site
                            site = ref_seq.var_aa_sites[counter]
                            continue

                        # Otherwise, proceed as we did for the forward read
                        variant_combo += var_pos
                        counter += 1

                        # Break the loop if we have reached the max
                        if counter == read_var_sites:
                            break

                        # Set the next site
                        site = ref_seq.var_aa_sites[counter]

                    # Add "-" on to the end of the variant for all missing sites
                    missing_sites = read_var_sites - counter
                    for _ in range(missing_sites):
                        variant_combo += "-"

                    # If mismatch is True, continue
                    if mismatch:
                        continue

                    # If this is the first time we've seen the combination add to the dictionary
                    if variant_combo not in variant_counts:
                        variant_counts[variant_combo] = 1

                    # Otherwise, add to the count
                    else:
                        variant_counts[variant_combo] += 1

                    # Add to the list of variants
                    variants.append([j, variant_combo])

        # Convert the count matrices to frequency matrices
        bp_sums = [None, None]
        aa_sums = [None, None]
        bp_freqs = [None, None]
        aa_freqs = [None, None]
        var_aa_info = [[], []]

        for i, (bp_count_matrix,
                aa_count_matrix,
                ref_seq) in enumerate(zip(bp_counts, aa_counts, ref_seqs)):

            bp_sums[i] = bp_count_matrix.sum(axis=1)
            aa_sums[i] = aa_count_matrix.sum(axis=1)
            bp_freqs[i] = np.nan_to_num(bp_count_matrix/bp_sums[i][:, None])
            aa_freqs[i] = np.nan_to_num(aa_count_matrix/aa_sums[i][:, None])

            # Identify the frequency of each amino acid at the variable positions
            for k, site in enumerate(ref_seq.var_aa_sites):

                # Identify all indices of the amino acid where we have more than 0 total counts
                non_zero_inds = np.argwhere(aa_count_matrix[site, :] != 0).flatten()

                # Report the frequency of amino acids at the site
                aa_info_dict = {ind_to_aa[ind]: aa_freqs[i][site, ind] for ind in non_zero_inds}

                # Store information about the site
                var_aa_info[i].append([k, aa_info_dict, aa_sums[i][site]])

        # Save all relevant information
        self.alignment_results = (insertions, deletions, bp_counts, aa_counts,
                                   bp_freqs, aa_freqs, var_aa_info, variant_counts)

    # Define a function for generating a consensus sequences
    def generate_consensus(self, freq_thresh = 0.9, forward = True):

        # Decide which set of alignment results we're working off of
        if forward:
            (insertions, deletions, bp_counts, aa_counts, bp_freq, aa_freq,
             var_aa_info) = [entry[0] for entry in self.alignment_results[:-1]]
        else:
            (insertions, deletions, bp_counts, aa_counts, bp_freq, aa_freq,
             var_aa_info) = [entry[1] for entry in self.alignment_results[:-1]]

        # Parse insertion information
        insertion_count_dict = {}
        insertion_info_dict = {}
        for j, ref_ind, ins_char in insertions:

            # Count the number of insertions of each type
            # If we haven't seen this reference index yet, add it to the dictionary
            if ref_ind not in insertion_info_dict:

                # Add to the count of the specific character at the position
                insertion_info_dict[ref_ind] = {ins_char: 1}

            # If we have seen the reference index but have not seen the character,
            # add the character to the dictionary
            elif ins_char not in insertion_info_dict[ref_ind]:

                # Add the specific character at the position
                insertion_info_dict[ref_ind][ins_char] = 1

            # Otherwise, just add to the existing instance
            else:
                insertion_info_dict[ref_ind][ins_char] += 1

            # If we haven't seen the reference yet, add to the count dictionary
            if ref_ind not in insertion_count_dict:
                insertion_count_dict[ref_ind] = 1

            # Otherwise, just add to the count
            else:
                insertion_count_dict[ref_ind] += 1

        # Get the frequency that an insertion is called at each site. Record the
        # insertion position if it is above the threshold
        true_insertions = []
        for ref_ind, count in insertion_count_dict.items():

            # Determine if the frequency of occurence is above the threshold
            if count/self.n_seqs >= freq_thresh:

                # Set a binary for recording whether an identified base meets the
                # required threshold
                ins_thresh_met = False

                # Take a look at the frequency of bases in the insertion.
                for bp, count in insertion_info_dict[ref_ind].items():

                    # If the character meets the threshold, record the character
                    if count/self.n_seqs >= freq_thresh:
                        true_insertions.append([ref_ind, bp])
                        ins_thresh_met = True

                # If the insertion threshold was not met for any characters,
                # insert an ambiguous bases
                if not ins_thresh_met:
                    true_insertions.append([ref_ind, "N"])

        # Loop over the bp_freq matrix and record the most common base at each
        # position assuming it occurs with greater frequency than the threshold
        thresh_met_inds = np.argwhere(bp_freq >= freq_thresh)
        max_inds = bp_freq.argmax(axis=1)

        # Identify indices where things are not ambiguous
        non_ambiguous_inds = set(thresh_met_inds[:, 0])

        # Construct the sequence
        consensus = np.array(["N" for ind in range(len(bp_freq))])
        consensus = [ind_to_bp[max_inds[ind]] if ind in non_ambiguous_inds else
                          "N" for ind in range(len(bp_freq))]

        # Add insertions
        for ind, bp in true_insertions:
            consensus.insert(ind, bp)

        # Decide whether to assign this as the forward or reverse consensus sequence
        if forward:
            self.f_consensus = "".join(consensus)
        else:
            self.r_consensus = "".join(consensus)

# Define a class that can handle plate information
class plate():

    """
    In report, print out all of the alignments to a text file
    1. Frequency of each base at each reference position.
    2. Any amino acid changes from the expected (including silent mutations) making
        sure to report the frequency of the changes.
    3. The amino acid at the variable position making sure to report the frequency
        of the amino acid at the site.
    4. Any frameshifts (deletions or insertions). All of these values can be appended
        as attributes of the well that we can report later in plate format.
    """

    # Define the initialization function. This defines the groups of wells as a single plate
    def __init__(self, wells):

        # Determine if we are running in troubleshoot mode or not

        # Set the list of wells as an attribute of the plate
        self.wells = wells

        # Sort the wells by their physical well location
        self.wells.sort(key = lambda x: x.well)

        # Set the plate name
        self.name = self.wells[0].plate_name

    # Write a function for processing contents of the plate. This will perform alignments
    # for each well, then generate a report
    def process(self, position, ts_mode, n_jobs, desc):

        # Create variables to store information needed for writing the output
        alignment_text = ""
        consensus_text = ""
        summary_info = [["Plate", "Well", "ReadDirection", "F-BC", "R-BC", "Site",
                         "AA", "AlignmentFrequency", "WellSeqDepth"]]
        bp_count_freq_info = []
        aa_count_freq_info = []
        variant_info = [["Plate", "Well", "F-BC", "R-BC",
                         "VariantCombo", "AlignmentFrequency", "WellSeqDepth"]]

        # Process all wells in the plate
        with Pool(n_jobs) as p:

            # Change approach depending on whether or not we are in TS mode
            if ts_mode:
                results = list(tqdm(p.imap(multiprocess_plate_analyzer_ts,
                                                     self.wells),
                                    position = position, desc = desc,
                                    total = len(self.wells)))

            else:
                results = list(tqdm(p.imap(multiprocess_plate_analyzer,
                                                     self.wells),
                                    position = position, desc = desc,
                                    total = len(self.wells)))

        # Process the results
        for result in results:

            # Process differently depending on whether or not we are running
            # in troubleshoot mode
            if ts_mode:

                # Unpack results
                (aa_count_freq_info_part, bp_count_freq_info_part,
                 alignment_text_part, consensus_text_part,
                 summary_info_part, variant_info_part) = result

                # Append the completed tables to the master lists
                aa_count_freq_info.extend(aa_count_freq_info_part)
                bp_count_freq_info.extend(bp_count_freq_info_part)

                # Extend the alignment text
                alignment_text += alignment_text_part

            else:

                # Unpack differently if we are not in troubleshoot mode
                consensus_text_part, summary_info_part, variant_info_part = result

            # Extend the summary and variant info tables
            summary_info.extend(summary_info_part)
            variant_info.extend(variant_info_part)

            # Extend the consensus sequences table
            consensus_text += consensus_text_part


        ###################### Save all of the generated data #################
        # Check to be sure we have all of the required directories. Make those
        # that we don't have

        default_dirs = [os.path.join(output_location, loc) for loc in ["Summaries",
                                                                "ConsensusSequences"]]
        for dir in default_dirs:

            # Make the directory if it doesn't exist
            if not os.path.isdir(dir):
                os.makedirs(dir)


        # If we are in troubleshoot mode, save the extra data
        if ts_mode:

            # Make additional directories
            extra_dirs = [os.path.join(output_location, loc) for loc in ["Alignments",
                                                                   "AACountsFrequencies",
                                                                     "BPCountsFrequencies",
                                                                     "Pickles"]]
            for dir in extra_dirs:

                # Make the directory if it doesn't exist
                if not os.path.isdir(dir):
                    os.makedirs(dir)

            # Save alignments
            with open(extra_dirs[0]+"/{}_Alignments.txt".format(self.name), "w") as f:
                f.writelines(alignment_text)

            # Save frequency and count tables
            with open(extra_dirs[1]+"/{}_AACountsFrequencies.csv".format(self.name), "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(aa_count_freq_info)

            with open(extra_dirs[2]+"/{}_BPCountsFrequencies.csv".format(self.name), "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(bp_count_freq_info)

        # Save summary info
        with open(default_dirs[0]+"/{}_SummaryInfo.csv".format(self.name), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(summary_info)

        # Save combination count information
        with open(default_dirs[0]+"/{}_VariantInfo.csv".format(self.name), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(variant_info)

        # Save consensus sequences
        with open(default_dirs[1]+"/{}_ConsensusSequences.txt".format(self.name), "w") as f:
            f.writelines(consensus_text)
