# Import third party modules
from Bio import pairwise2

# Import ssSeqSupport modules
from . import AaOpts, BpOpts

#Define a function to analyze wells in plate for bp frequencies only
def MultiProcessPlateAnalyzerTiles(args):
    """
    In ssSeqMultiprocessing.py
    
    AMK
    
    Analyze wells in the plate for their basepair frequencies.
    Used for determining the positions mutated in each well.

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #print("start of multiprocessplateanalyzertiles")
    
    # Unpack args
    well, alignment_cutoff, q_cutoff = args

    # Make and process the alignments
    well.align()
    
    #print("Before analyze_alignments_tiles")
    #Replace analyze alignments with generate frequency matrices
    well.analyze_alignments_tiles(alignment_cutoff, q_cutoff)
    #print("After analyze_alignments_tiles")
    
    
    # Generate consensus bp sequences for consensus outputs
    well.generate_consensus()
    well.generate_consensus(forward=False)
    
    # AMK function, determine mutation positions by comparing top frequencies
    # to the reference sequence
    
    #print("Before compare_freq_to_ref")
    well.compare_freq_to_ref()
    #print("After 1 compare_freq_to_ref")
    well.compare_freq_to_ref(forward=False)
    #print("After both compare_freq_to_ref")
    
    # Generate axis labels for matrices
    f_ref_seq_aas = well._f_ref_seq.aas_as_list
    r_ref_seq_aas = well._r_ref_seq.aas_as_list
    f_ref_seq_bps = well._f_ref_seq.bps_as_list
    r_ref_seq_bps = well._r_ref_seq.bps_as_list

    # # Get the variant counts
    # variant_counts = well._alignment_results[-1]

    # # Get the depth (the sum of all counts in the variant)
    # var_depth = sum(variant_counts.values())
    


    # # Create the well-specific variant counts table
    # variant_info = [[well.plate, well.well, well.f_barcode, well.r_barcode,
    #                   combo, count/var_depth, var_depth]
    #                 for combo, count in variant_counts.items()]

    # Generate the header
    header = "{}_{}_{}_{}".format(well.plate, well.well,
                                  well.f_barcode, well.r_barcode)

    # Create the well-specific summary tables
    summary_info = []
    
    #AMK added posiiton summary info
    pos_summary_info = []
    
    
    bp_count_freq_info = []
    aa_count_freq_info = []

    # Analyze both forward and reverse results
    for (label, aa_labels, bp_labels, _, _, bp_counts,
         aa_counts, bp_freq, aa_freq, var_aa_info, var_pos_info) in zip(
                                             ("Forward", "Reverse"),
                                             (f_ref_seq_aas, r_ref_seq_aas),
                                             (f_ref_seq_bps, r_ref_seq_bps),
                                              *well._alignment_results[:-1],
                                              (well.f_var_pos_info,
                                              well.r_var_pos_info)):

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
        aa_table.extend([[aa] + sublist for aa, sublist in zip(AaOpts, aa_count_list)])
        aa_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         aa_labels])
        aa_table.extend([[aa] + sublist for aa, sublist in zip(AaOpts, aa_freq_list)])
        aa_table.extend([[None], [None]])

        # Build the basepair table
        bp_table = [[header + "_Counts_" + label],
                    bp_labels]
        bp_table.extend([[bp] + sublist for bp, sublist in zip(BpOpts, bp_count_list)])
        bp_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         bp_labels])
        bp_table.extend([[bp] + sublist for bp, sublist in zip(BpOpts, bp_freq_list)])
        bp_table.extend([[None], [None]])

        # Append the completed tables to the master lists
        aa_count_freq_info.extend(aa_table)
        bp_count_freq_info.extend(bp_table)

        # # Append to the summary table
        # for site, freq_dict, depth in var_aa_info:
        #     summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
        #                           site + 1, aa, freq, depth] for aa, freq in freq_dict.items()])
        
        #print(var_pos_info)
        # # Append to the summary table
        # for refcodon, bp_index, codon, refaa, aa_position, aa in var_pos_info:
        #     pos_summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
        #                           refcodon, bp_index, codon, refaa, aa_position, aa]])
        
        # debug
        # Append to the summary table
        for i in var_pos_info:
            pos_summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
                                  i]])


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
            consensus_text, summary_info, pos_summary_info)
    
    
    

# Define a function that analyzes a well in the plate when
# we are in troubleshooting mode
def MultiprocessPlateAnalyzerTS(args):
    """
    In SupportFuncs.py
    
    Anlyzes wells in the plate when in troubleshooting mode.

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    aa_count_freq_info : TYPE
        DESCRIPTION.
    bp_count_freq_info : TYPE
        DESCRIPTION.
    alignment_text : TYPE
        DESCRIPTION.
    consensus_text : TYPE
        DESCRIPTION.
    summary_info : TYPE
        DESCRIPTION.
    variant_info : TYPE
        DESCRIPTION.

    """
    
    # Unpack args
    well, alignment_cutoff, q_cutoff = args

    # Make and process the alignments
    well.align()
    well.analyze_alignments(alignment_cutoff, q_cutoff)
    well.generate_consensus()
    well.generate_consensus(forward=False)

    # Generate axis labels for matrices
    f_ref_seq_aas = well._f_ref_seq.aas_as_list
    r_ref_seq_aas = well._r_ref_seq.aas_as_list
    f_ref_seq_bps = well._f_ref_seq.bps_as_list
    r_ref_seq_bps = well._r_ref_seq.bps_as_list

    # Get the variant counts
    variant_counts = well._alignment_results[-1]

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
    for (label, aa_labels, bp_labels, _, _, bp_counts,
         aa_counts, bp_freq, aa_freq, var_aa_info) in zip(("Forward", "Reverse"),
                                                       (f_ref_seq_aas, r_ref_seq_aas),
                                                       (f_ref_seq_bps, r_ref_seq_bps),
                                                       *well._alignment_results[:-1]):

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
        aa_table.extend([[aa] + sublist for aa, sublist in zip(AaOpts, aa_count_list)])
        aa_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         aa_labels])
        aa_table.extend([[aa] + sublist for aa, sublist in zip(AaOpts, aa_freq_list)])
        aa_table.extend([[None], [None]])

        # Build the basepair table
        bp_table = [[header + "_Counts_" + label],
                    bp_labels]
        bp_table.extend([[bp] + sublist for bp, sublist in zip(BpOpts, bp_count_list)])
        bp_table.extend([[None],
                         [header + "_Frequencies_" + label],
                         bp_labels])
        bp_table.extend([[bp] + sublist for bp, sublist in zip(BpOpts, bp_freq_list)])
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
def MultiprocessPlateAnalyzer(args):
    """
    In SupportFuncs.py
    
    Analyzes a well in the plate when not in troubleshooting mode

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    summary_info : TYPE
        DESCRIPTION.
    variant_info : TYPE
        DESCRIPTION.

    """

    # Unpack args
    well, alignment_cutoff, q_cutoff = args

    # Make and process the alignments
    well.align()
    well.analyze_alignments(alignment_cutoff, q_cutoff)

    # Generate axis labels for matrices
    f_ref_seq_aas = well._f_ref_seq.aas_as_list
    r_ref_seq_aas = well._r_ref_seq.aas_as_list
    f_ref_seq_bps = well._f_ref_seq.bps_as_list
    r_ref_seq_bps = well._r_ref_seq.bps_as_list

    # Get the variant counts
    variant_counts = well._alignment_results[-1]

    # Get the depth (the sum of all counts in the variant)
    var_depth = sum(variant_counts.values())

    # Create the well-specific variant counts table
    variant_info = [[well.plate, well.well, well.f_barcode, well.r_barcode,
                     combo, count/var_depth, var_depth]
                    for combo, count in variant_counts.items()]

    # Create the well-specific summary tables
    summary_info = []

    # Analyze both forward and reverse results
    for (label, _, _, _, _, _, _, _, _, var_aa_info) in zip(("Forward", "Reverse"),
                                                            (f_ref_seq_aas, r_ref_seq_aas),
                                                            (f_ref_seq_bps, r_ref_seq_bps),
                                                            *well._alignment_results[:-1]):
        # Append to the summary table
        for site, freq_dict, depth in var_aa_info:
            summary_info.extend([[well.plate, well.well, label, well.f_barcode, well.r_barcode,
                                  site + 1, aa, freq, depth] for aa, freq in freq_dict.items()])

    # Return all information
    return summary_info, variant_info
