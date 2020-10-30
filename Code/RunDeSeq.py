warnings.warn("RunDeSeq must still be integrated with RunssSeq")

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
            


# Write a function for filtering out bad seqpairs
def qc_seqpairs(all_seqpairs, read_length = None, length_cutoff = 0.9, 
                average_q_cutoff = 25):
    
    print("Running read qc...")
    
    # If we don't have the read length determine it
    if read_length is None:

        # Get the most common read length. We will assign this as our read length
        all_readlengths = np.array([seqpair.read_lengths() for seqpair in all_seqpairs])
        read_length = ss.mode(all_readlengths, axis = None, nan_policy = "omit").mode[0]
        
    # Calculate the read filter
    read_filter = read_length * length_cutoff
        
    # Run QC on every read
    for seqpair in all_seqpairs:
        seqpair.qc_reads(read_filter, average_q_cutoff)
    
    # Eliminate any duds, which are those seqpairs with both a forward and a reverse that failed qc
    no_duds = tuple(filter(lambda x: not x.is_dud(), all_seqpairs))
    
    return no_duds

# Write a function for assigning seqpairs to a well
def assign_seqpairs_to_well(filtered_seqpairs, bc_to_ref_plate_well, savedir):

    # Loop over all seqpairs and assign to wells
    print("Assigning sequences to wells...")
    well_pairs = {}
    for pair in filtered_seqpairs:

        # Grab the well ID and see if it is a real well. Continue
        # if it is not. "Fake" wells are those that result from 
        # sequencing errors
        well_id = (pair.f_barcode, pair.r_barcode)
        if well_id not in bc_to_ref_plate_well:
            continue
        
        # Check to see if we have seen this well already.
        # If we have seen it, append to growing list. If we have not,
        # start a new list
        if well_id in well_pairs:
            well_pairs[well_id].append(pair)
        else:
            well_pairs[well_id] = [pair]
            
    # Now build and return the well objects
    return [Well(pair, bc_to_ref_plate_well[well_id], savedir) 
            for well_id, pair in well_pairs.items()] 

# Write a function that can process a single well
def process_well(well, return_alignments = False, 
                 bp_q_cutoff = 30,
                 variable_thresh = 0.1,
                 variable_count = 1):

    # Align
    well.align()

    # Analyze alignments. 
    has_reads = well.analyze_alignments(bp_q_cutoff, variable_count)

    # If we don't have any reads that passed
    # QC, skip straight to analyzing counts. 
    if has_reads:
        # Build count matrices
        well.build_unit_count_matrices()

        # Identify variable positions
        well.identify_variable_positions(variable_thresh)

    # Analyze reads with decoupled counts
    well.analyze_unpaired_counts(variable_thresh)

    # Analyze reads with coupled counts
    well.analyze_paired_counts(variable_thresh, variable_count)

    # If we are returning alignments, generate them
    if return_alignments:
        formatted_alignments = well.format_alignments()
    else:
        formatted_alignments = None

    # Return relevant information for downstream processing
    return (well.unpaired_bp_output, well.unpaired_bp_output_max,
            well.unpaired_aa_output, well.unpaired_aa_output_max,
            well.paired_bp_output, well.paired_aa_output, 
            formatted_alignments) 
    
def format_and_save_outputs(well_results, saveloc, return_alignments):

    # Write a function that processes the output of analyzing all wells
    unpacked_output = tuple(zip(*well_results))

    # Concatenate all dataframes
    full_dfs = tuple(pd.concat(df_list, ignore_index = True) for df_list in unpacked_output[:-1])

    # Get just the max of each of the paired outputs
    max_outs = [None, None]
    for i, paired_df in enumerate(full_dfs[4:]):

        # Get the columns of interest
        limited_df = paired_df.loc[:, ["IndexPlate", "Well", "AlignmentFrequency"]]

        # Group by plate and well
        grouped_df = limited_df.groupby(by = ["IndexPlate", "Well"])
        max_outs[i] = paired_df.loc[grouped_df.idxmax().AlignmentFrequency.values].copy()

    # Loop over all dataframes, sort by plate and well, and save
    savenames = ("Bases_Decoupled_All.csv", "Bases_Decoupled_Max.csv",
                "AminoAcids_Decoupled_All.csv", "AminoAcids_Decoupled_Max.csv",
                 "Bases_Coupled_All.csv", "Combos_Coupled_All.csv",
                 "Bases_Coupled_Max.csv", "Combos_Coupled_Max.csv")
    for savename, output_df in zip(savenames, chain(full_dfs, max_outs)):

        # Sort by plate and well
        output_df.sort_values(by = ["IndexPlate", "Well"], inplace = True)

        # Save the dataframe
        output_df.to_csv(os.path.join(saveloc, "OutputCounts", savename), index = False)

    # Loop over and save all alignments if asked to do so
    if return_alignments:
        for savename, savestr in unpacked_output[-1]:
            with open(savename, "w") as f:
                f.write(savestr)    

# Write a function that runs deSeq
def run_deseq(ref_seq_loc, forward_read_loc, reverse_read_loc,
              saveloc, n_cpus, stop_after_qualities = False,
              stop_after_fastq = False, return_alignments = False,
              read_length = None, length_cutoff = 0.9,
              average_q_cutoff = 25, bp_q_cutoff = 30,
              variable_thresh = 0.1, variable_count = 1):
    
    # Load the reference sequence file and associate reference sequence
    # information with wells
    bc_to_ref_plate_well = load_refseq(ref_seq_loc)
    
    # Load fastq files
    all_seqpairs = load_fastq(forward_read_loc, reverse_read_loc)
    
    # Filter the seqpairs
    filtered_seqpairs = qc_seqpairs(all_seqpairs, 
                                    read_length = read_length,
                                    length_cutoff = length_cutoff, 
                                    average_q_cutoff = average_q_cutoff)
    
    # Plot qualities
    warnings.warn("Need to implement quality plotting")
#     plot_qualities()
    
    # Return if we stop after plot qualities
    if stop_after_qualities:
        return
    
    # Assign seqpairs to a well
    all_wells = assign_seqpairs_to_well(filtered_seqpairs, bc_to_ref_plate_well, saveloc)
    
    # Save the fastq files
    for well in all_wells:
        well.write_fastqs()

    # Return if we stop after fastq
    if stop_after_fastq:
        return
        
    # Complete the multiprocessing function, then process all wells
    complete_multiprocessor = partial(process_well, 
                                      bp_q_cutoff = bp_q_cutoff,
                                      return_alignments = return_alignments,
                                     variable_thresh = variable_thresh,
                                     variable_count = variable_count)
        
    # Multiprocess to handle wells
    with Pool(n_cpus) as p:
        processed_well_results = list(tqdm(p.imap_unordered(complete_multiprocessor, all_wells),
                                      desc = "Processing wells:", total = len(all_wells)))
        
    # Handle processed output
    format_and_save_outputs(processed_well_results, saveloc, return_alignments)