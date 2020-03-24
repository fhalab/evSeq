# Import third party modules
from Bio import pairwise2
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import os.path
import csv

# Import necessary modules from ssSeqSupport
from . import BuildRefSeqs
from . import BpToInd, AaToInd, IndToAa, IndToBp
from . import LogError, LogWarning
from . import Translation
from . import MultiprocessPlateAnalyzer, MultiprocessPlateAnalyzerTS
from . import GenerateSequencingHeatmap

# Write a class which holds all sequence blocks of a given set of barcodes
class Well():

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
    def __init__(self, seq_pairs, read_length, BcsToRefSeq):

        # Set the group of sequences as a class attribute
        self._seq_pairs = seq_pairs

        # Identify the forward and reverse barcodes for the object
        self._f_barcode = seq_pairs[0].f_barcode
        self._r_barcode = seq_pairs[0].r_barcode

        # Confirm that the barcode pair exists in the barcode dictionary
        bc_pair = (seq_pairs[0].f_barcode, seq_pairs[0].r_barcode)
        if bc_pair in BcsToRefSeq:

            # Identify the barcode plate, plate nickname, well string, and 
            # raw reference sequence
            self._plate, self._well, self._plate_name, raw_ref = BcsToRefSeq[bc_pair]
            
            # Construct the forward and reverse reference sequences
            self._f_ref_seq, self._r_ref_seq, self._n_variable_sites = BuildRefSeqs(raw_ref,
                                                                                    read_length)
            
            # Confirm that this is a real well
            self._real_well = True

        # If the barcode pair does not exist, then this is not a real well
        # It likely results from a sequencing error
        else:
            self._real_well = False

        # Get the number of sequences paired to the well
        self._n_seqs = len(seq_pairs)

    # Define a function for aligning each sequence to the appropriate reference
    def align(self):

        # Loop over the sequence pairs and make alignments to the reference sequences
        self._f_alignments = [pairwise2.align.globalxs(self._f_ref_seq.seq,
                                                      pair.barcodeless_f,
                                                      -2, -1, one_alignment_only=True,
                                                      penalize_end_gaps = False)[0]
                             for pair in self._seq_pairs]
        self._r_alignments = [pairwise2.align.globalxs(self._r_ref_seq.seq,
                                                             pair.reversed_barcodeless_r,
                                                             -2, -1, one_alignment_only=True,
                                                             penalize_end_gaps = False)[0]
                             for pair in self._seq_pairs]

    # Write a helper function that analyzes a single alignment
    def analyze_alignment(self, alignment, j, f_r_ind, bp_counts, 
                          alignment_cutoff, q_cutoff):

        # Unpack the alignment
        (reference_alignment, read_alignment, score, _, _) = alignment

        # If the alignment score is below the threshold, return this information
        # and move on
        if f_r_ind==0:
            if score/self._f_ref_seq.seq_length < alignment_cutoff:
                return None
        else:
            if score/self._r_ref_seq.seq_length < alignment_cutoff:
                return None

        # Define lists to store identified insertions and deletions
        insertions = []
        deletions = []

        # Define lists for holding temporary information
        holding_dels = []

        # Define a list for holding low quality characters
        low_quality_chars = []

        # Set the reference sequence
        if f_r_ind == 0:
            ref_seq = self._f_ref_seq
        else:
            ref_seq = self._r_ref_seq

        # Define variables that will be used for keeping track of our position
        # along the alignment
        ref_counter = -1
        read_counter = -1
        del_possible = False
        ins_found = False
        del_found = False
        in_ref = False
        stop_point = ref_seq.seq_length - 1

        # Loop over each base in the alignments
        for ref_char, read_char in zip(reference_alignment, read_alignment):

            # If the read counter is not a space, add 1.
            if read_char != "-":
                read_counter += 1

            # Until we have a non-space character in the reference, we do not 
            # record information. We only increment the reference counter when
            # we have non-space characters
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
            if ref_counter == stop_point:
                bp_counts[ref_counter, BpToInd[read_char]] += 1
                break

            # Determine if the ref char is a space and we are in the reference 
            # sequence. If it is, this indicates an insertion in the read sequence.
            if ref_char=="-" and in_ref:
                insertions.append([j, ref_counter, read_char])

                # Report that an insertion was found
                ins_found = True

                # Continue to the next character
                continue

            # Determine if the read character surpasses our quality score filter.
            # If it does not, record and continue
            if f_r_ind == 0:
                if self._seq_pairs[j].barcodeless_f_q[read_counter] < q_cutoff:
                    low_quality_chars.append(read_counter)
                    continue
            else:
                if self._seq_pairs[j].reversed_barcodeless_r_q[read_counter] < q_cutoff:
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
            bp_counts[ref_counter, BpToInd[read_char]] += 1

        # Return relevant values
        return insertions, deletions, read_trans_start, ins_found, del_found, low_quality_chars

    # Define a function for performing data analysis on the alignments
    def analyze_alignments(self, alignment_cutoff, q_cutoff):

        ######################## Upfront Processing #############################
        # Determine the number of alignments in the set
        n_f_alignments = len(self._f_alignments)
        n_r_alignments = len(self._r_alignments)
        
        # Log an error if we don't have an equal number of alignments
        if n_f_alignments != n_r_alignments:
            LogError("There are an unequal number of forward and reverse alignments.")
    
        # Build a matrix to count the occurences of each base type at each position
        # in the reference sequence
        bp_counts = [np.zeros([self._f_ref_seq.seq_length, 6]),
                     np.zeros([self._r_ref_seq.seq_length, 6])]
        aa_counts = [np.zeros([self._f_ref_seq.trans_length, 22]),
                     np.zeros([self._r_ref_seq.trans_length, 22])]

        # Package redundant arguments
        ref_seqs = [self._f_ref_seq, self._r_ref_seq]

        # Identify the number of variable positions in each reference
        f_read_var_sites = len(self._f_ref_seq.var_aa_sites)
        r_read_var_sites = len(self._r_ref_seq.var_aa_sites)
        both_read_var_sites = [f_read_var_sites, r_read_var_sites]

        # Get the number of overlapping sites between the two
        n_overlapping_sites = f_read_var_sites + r_read_var_sites - self._n_variable_sites

        # Loop over each alignment in the set and analyze
        insertions = [[],[]]
        deletions = [[],[]]
        variant_counts = {}
        for j, (f_alignment, r_alignment) in enumerate(zip(self._f_alignments, self._r_alignments)):

            # Pull the unaligned read sequences
            unaligned_read_seqs = [self._seq_pairs[j].barcodeless_f,
                                   self._seq_pairs[j].reversed_barcodeless_r]

            # Create a series of tests
            check_combo = True
            bad_codon_found = False

            # Analyze the forward and reverse alignments
            for i, alignment in enumerate((f_alignment, r_alignment)):

                # Analyze the alignment
                alignment_output = self.analyze_alignment(alignment, j, i, 
                                                          bp_counts[i], 
                                                          alignment_cutoff,
                                                          q_cutoff)

                # If the alignment did not meet the alignment score cutoff, 
                # record it as a bad alignment
                if alignment_output is None:
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
                if ins_found or del_found:
                    check_combo = False
                    continue

                # Translate each codon in the read alignment
                alignment_translation = Translation(unaligned_read_seqs[i],
                                                    read_trans_start,
                                                    low_quality_chars)

                # Point to the appropriate reference sequence
                ref_seq = ref_seqs[i]

                # Add to the amino acid count matrix and compare the translation to the
                # reference
                for k, aa_read in enumerate(alignment_translation.translation[:ref_seq.trans_length]):

                    # Add to the amino acid count matrix if it is not a low quality
                    # codon. If this codon is a variable position, continue adding
                    # to the alignment matrix, but don't try and determine the
                    # variant. Just continue in the parent loop from here
                    if k in alignment_translation.low_quality_codons:
                        if k in ref_seq.var_aa_sites:
                            bad_codon_found = True
                        continue
                    aa_counts[i][k, AaToInd[aa_read]] += 1

                # If the forward reads failed for any reason, continue. If the 
                # # reverse reads failed for any reason, break the loop.
                # We can't make combination data without both reads.
                if (bad_codon_found or not check_combo) and i == 0:
                    continue
                elif (bad_codon_found or not check_combo) and i == 1:
                    break

                # Identify the highest possible reference and the number of variable
                # sites in the read
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
                aa_info_dict = {IndToAa[ind]: aa_freqs[i][site, ind] for ind in non_zero_inds}

                # Store information about the site
                var_aa_info[i].append([k, aa_info_dict, aa_sums[i][site]])

        # Save all relevant information
        self._alignment_results = (insertions, deletions, bp_counts, aa_counts,
                                   bp_freqs, aa_freqs, var_aa_info, variant_counts)

    # Define a function for generating a consensus sequences
    def generate_consensus(self, freq_thresh = 0.9, forward = True):

        # Decide which set of alignment results we're working off of
        if forward:
            (insertions, _, _, _, bp_freq, _, _) = [entry[0] for entry in 
                                                    self._alignment_results[:-1]]
        else:
            (insertions, _, _, _, bp_freq, _, _) = [entry[1] for entry in
                                                     self._alignment_results[:-1]]

        # Parse insertion information
        insertion_count_dict = {}
        insertion_info_dict = {}
        for _, ref_ind, ins_char in insertions:

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
            if count/self._n_seqs >= freq_thresh:

                # Set a binary for recording whether an identified base meets the
                # required threshold
                ins_thresh_met = False

                # Take a look at the frequency of bases in the insertion.
                for bp, count in insertion_info_dict[ref_ind].items():

                    # If the character meets the threshold, record the character
                    if count/self._n_seqs >= freq_thresh:
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
        consensus = [IndToBp[max_inds[ind]] if ind in non_ambiguous_inds else
                          "N" for ind in range(len(bp_freq))]

        # Add insertions
        for ind, bp in true_insertions:
            consensus.insert(ind, bp)

        # Decide whether to assign this as the forward or reverse consensus sequence
        if forward:
            self._f_consensus = "".join(consensus)
        else:
            self._r_consensus = "".join(consensus)
            
    # Set up properties
    @property
    def well(self):
        return self._well
    
    @property
    def real_well(self):
        return self._real_well
    
    @property
    def plate_name(self):
        return self._plate_name
    
    @property
    def f_barcode(self):
        return self._f_barcode
    
    @property
    def r_barcode(self):
        return self._r_barcode
    
    @property
    def plate(self):
        return self._plate

# Define a class that can handle plate information
class Plate():

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
        self._wells = wells

        # Sort the wells by their physical well location
        self._wells.sort(key = lambda x: x.well)

        # Set the plate name
        self._name = self._wells[0].plate_name

    # Write a function for processing contents of the plate. This will perform alignments
    # for each well, then generate a report
    def process(self, args, desc, filepair, combo_ind):

        # Create variables to store information needed for writing the output
        alignment_text = ""
        consensus_text = ""
        summary_info = [["Plate", "Well", "ReadDirection", "F-BC",
                         "R-BC", "Site", "AA", "AlignmentFrequency", "WellSeqDepth"]]
        bp_count_freq_info = []
        aa_count_freq_info = []
        variant_info = [["Plate", "Well", "F-BC", "R-BC",
                         "VariantCombo", "AlignmentFrequency", "WellSeqDepth"]]

        # Process all wells in the plate
        with Pool(args["jobs"]) as p:

            # Generate multiprocessing args
            multiprocessess_args = [[well, args["alignment_filter"], args["q_cutoff"]]
                                    for well in self._wells]

            # Change approach depending on whether or not we are in TS mode
            if args["troubleshoot"]:
                results = list(tqdm(p.imap(MultiprocessPlateAnalyzerTS,
                                           multiprocessess_args),
                                    position = 2, desc = desc,
                                    total = len(self._wells), leave = False))

            else:
                results = list(tqdm(p.imap(MultiprocessPlateAnalyzer,
                                           multiprocessess_args),
                                    position = 2, desc = desc,
                                    total = len(self._wells), leave = False))

        # Process the results
        for result in results:

            # Process differently depending on whether or not we are running
            # in troubleshoot mode
            if args["troubleshoot"]:

                # Unpack results
                (aa_count_freq_info_part, bp_count_freq_info_part,
                 alignment_text_part, consensus_text_part,
                 summary_info_part, variant_info_part) = result

                # Append the completed tables to the master lists
                aa_count_freq_info.extend(aa_count_freq_info_part)
                bp_count_freq_info.extend(bp_count_freq_info_part)

                # Extend the alignment text
                alignment_text += alignment_text_part
                
                # Extend the consensus sequences table
                consensus_text += consensus_text_part


            else:

                # Unpack differently if we are not in troubleshoot mode
                summary_info_part, variant_info_part = result
 
            # Extend the summary and variant info tables
            summary_info.extend(summary_info_part)
            variant_info.extend(variant_info_part)

        # Convert summary info and variant infor to dataframes
        summary_info = pd.DataFrame(summary_info[1:], columns = summary_info[0])
        variant_info_df = pd.DataFrame(variant_info[1:], columns = variant_info[0])
        
        # Add the source file names on
        summary_info["R1"] = filepair[0]
        summary_info["R2"] = filepair[1]
        variant_info_df["R1"] = filepair[0]
        variant_info_df["R2"] = filepair[1]

        # Merge with full DataFrame unless we had no alignment. If there is no
        # alignment, then skip the merge and assign df_full to be a copy of
        # variant_info_df
        if len(variant_info_df) == 0:
            
            # Just make a copy of df_full
            df_full = variant_info_df.copy()
            
            # Log a warning that we didn't find any alignments for this plate
            LogWarning("\nNo sequences passed alignment QC for '{}'. Consider lowering filtering thresholds.".format(self._name))
            
        else:
            
            # Find the max value for Alignment Frequency for each well
            df_max = variant_info_df.groupby('Well')[['Well', 'AlignmentFrequency']].max().reset_index(drop=True)
            df_full = variant_info_df.merge(df_max)

            # Round for easier reading
            df_full['AlignmentFrequency'] = np.round(df_full['AlignmentFrequency'].values, 3)
            
            # Generate heatmap and save it
            # Build the heatmaps folder
            hm_output_file = os.path.join(args["output"],
                                            "Platemaps/{}-{}_SequencingHeatmap".format(self._name,
                                                                                        combo_ind))
            GenerateSequencingHeatmap(df_full, self._name, hm_output_file)


        ###################### Save all of the generated data #################
        # Get the paths to all save locations
        summary_dir = os.path.join(args["output"], "Summaries")
        extra_dirs = [os.path.join(args["output"], loc) for loc in
                      ["Alignments", "AACountsFrequencies",
                       "BPCountsFrequencies", "ConsensusSequences"]]
        
        # If we are in troubleshoot mode, save the extra data
        if args["troubleshoot"]:

            # Save alignments
            with open(extra_dirs[0]+"/{}-{}_Alignments.txt".format(self._name,
                                                                   combo_ind), "w") as f:
                f.writelines(alignment_text)

            # Save frequency and count tables
            with open(extra_dirs[1]+"/{}-{}_AACountsFrequencies.csv".format(self._name,
                                                                            combo_ind), "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(aa_count_freq_info)

            with open(extra_dirs[2]+"/{}-{}_BPCountsFrequencies.csv".format(self._name,
                                                                            combo_ind), "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(bp_count_freq_info)
                                
            # Save consensus sequences
            with open(extra_dirs[3]+"/{}-{}_ConsensusSequences.txt".format(self._name,
                                                                           combo_ind), "w") as f:
                f.writelines(consensus_text)

        # Save summary info
        summary_info.to_csv(summary_dir+"/{}-{}_SummaryInfo.csv".format(self._name,
                                                                           combo_ind),
                               index = False)

        
        # Save combination count information
        variant_info_df.to_csv(summary_dir+"/{}-{}_VariantInfo.csv".format(self._name,
                                                                           combo_ind),
                               index = False)
           
        # Save the limited table with which captures the maximum frequency alignment    
        df_full.to_csv(summary_dir+"/{}-{}_MaxInfo.csv".format(self._name, combo_ind),
                       index=False)
        