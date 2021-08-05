# Import evSeq objects
from .util.logging import log_info, log_warning
from .util.input_validation import check_args
from .util.input_processing import load_all
from .seq_pair import SeqPair
from .well import Well
from .data_visualization import (generate_read_qual_chart, 
                                 generate_sequencing_heatmaps,
                                 save_heatmap_to_file)

# Import other required modules
import os
import numpy as np
import pandas as pd
import scipy.stats as ss
from Bio import SeqIO
from itertools import chain
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import gzip


def build_seqpairs(f_loc, r_loc):
    """Load and pair fastq file entries."""

    # Create a dictionary that links id to sequence object
    id_to_reads = {}
    
    # Load and parse forward reads
    print("Loading forward reads...")

    # Determine how to open file if zipped or not
    f_opener, f_encoding = (gzip.open, 'rt') if "fastq.gz" in f_loc \
                      else (open, 'r')

    # Store as list in memory
    with f_opener(f_loc, f_encoding) as handle:
        all_f_recs = list(SeqIO.parse(handle, "fastq"))

    for f_record in tqdm(all_f_recs, desc = "Parsing forward reads..."):
        temp_record = SeqPair()
        temp_record.assign_f(f_record)
        id_to_reads[f_record.id] = temp_record
    
    # Associate reverse reads with the forward
    print("Loading reverse reads...")

    # Determine how to open
    r_opener, r_encoding = (gzip.open, 'rt') if "fastq.gz" in r_loc \
        else (open, 'r')

    # Store as list
    with r_opener(r_loc, r_encoding) as handle:
        all_r_recs = list(SeqIO.parse(handle, "fastq"))
    
    for r_record in tqdm(all_r_recs, "Pairing reverse reads..."):

        # If there is no partern in id_to_reads, create a new object 
        # and continue
        if r_record.id not in id_to_reads:
            temp_record = SeqPair()
            temp_record.assign_r(r_record)
            id_to_reads[r_record.id] = temp_record

        # Otherwise, attach the reverse record
        else:
            id_to_reads[r_record.id].assign_r(r_record)
            
    # Return all records
    return tuple(id_to_reads.values())


def qc_seqpairs(
    all_seqpairs,
    read_length,
    length_cutoff,
    average_q_cutoff
):
    """Filters out bad seqpairs."""

    print("Running read qc...")
    
    # If we don't have the read length determine it
    if read_length is None:

        # Get the most common read length. We will assign this as our read length
        all_readlengths = np.array([seqpair.read_lengths() for seqpair in all_seqpairs])
        read_length = ss.mode(all_readlengths, axis = None, nan_policy = "omit").mode[0]
        
        log_info(f"A read length of {read_length} was calculated for this run")
        
    # Calculate the read filter
    read_filter = read_length * length_cutoff
        
    # Run QC on every read
    for seqpair in all_seqpairs:
        seqpair.qc_reads(read_filter, average_q_cutoff)
    
    # Eliminate any duds, which are those seqpairs with both a forward and a reverse that failed qc
    no_duds = tuple(filter(lambda x: not x.is_dud(), all_seqpairs))
    
    return no_duds


def assign_seqpairs_to_well(filtered_seqpairs, bc_to_ref_plate_well, savedir):
    """Assigning seqpairs to a well."""

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


def process_well(
    well,
    return_alignments=False,
    bp_q_cutoff=30,
    variable_thresh=0.1,
    variable_count=1,
):
    """Processes a single well."""
    
    # Define a frequency warning variable
    freq_warning = ""
    
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
        freq_warning = well.identify_variable_positions(variable_thresh)

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
            formatted_alignments, freq_warning) 
    
def format_and_save_outputs(well_results, saveloc, return_alignments):
    """Processes the output of analyzing all wells and saves results
    to the disk."""
    unpacked_output = tuple(zip(*well_results))

    # Concatenate all dataframes
    full_dfs = tuple(pd.concat(df_list, ignore_index = True) for df_list in unpacked_output[:-2])

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
                 "Bases_Coupled_All.csv", "AminoAcids_Coupled_All.csv",
                 "Bases_Coupled_Max.csv", "AminoAcids_Coupled_Max.csv")
    for savename, output_df in zip(savenames, chain(full_dfs, max_outs)):

        # Sort by plate and well
        output_df = output_df.sort_values(by=["IndexPlate", "Well"])

        # Save the dataframe
        output_df.to_csv(os.path.join(saveloc, "OutputCounts", savename), index=False)
        
    # Generate heatmaps from the Combos_Coupled_Max dataframe
    heatmaps = generate_sequencing_heatmaps(max_outs[-1])
    save_heatmap_to_file(heatmaps, saveloc)

    # Loop over and save all alignments if asked to do so
    if return_alignments:
        for savename, savestr in unpacked_output[-2]:
            with open(savename, "w") as f:
                f.write(savestr)   
                
    # Report all warnings
    for warning in unpacked_output[-1]:
        if warning != "":
            log_warning(f"High mutational frequency in {warning}. You may "
                        "want to check alignments for accuracy.") 


def run_evSeq(cl_args):
    """Runs the evSeq protocol. Should be run from command line."""
    # Check the input arguments
    check_args(cl_args)
    
    # Identify sequencing files and load reference sequence information
    forward_file, reverse_file, bc_to_ref_plate_well = load_all(cl_args)
    
    # Pair all sequences
    all_seqpairs = build_seqpairs(forward_file, reverse_file)
    
    # Generate quality plots
    generate_read_qual_chart(all_seqpairs, cl_args["output"])
    
    # Return if we stop after plot qualities
    if cl_args["analysis_only"]:
        return
    
    # Run QC on the seqpairs
    filtered_seqpairs = qc_seqpairs(all_seqpairs, cl_args["read_length"],
                                    cl_args["length_cutoff"], 
                                    cl_args["average_q_cutoff"])
    
    # Assign seqpairs to a well
    all_wells = assign_seqpairs_to_well(filtered_seqpairs, 
                                        bc_to_ref_plate_well,
                                        cl_args["output"])
    
    # Save the fastq files
    if cl_args["keep_parsed_fastqs"] or cl_args["only_parse_fastqs"]:
        print("Saving parsed well-filtered fastq files...")
        for well in all_wells:
            well.write_fastqs()

        # Return if we stop after fastq
        if cl_args["only_parse_fastqs"]:
            return
        
    # Complete the multiprocessing function, then process all wells
    complete_multiprocessor = partial(process_well, 
                                      bp_q_cutoff = cl_args["bp_q_cutoff"],
                                      return_alignments = cl_args["return_alignments"],
                                      variable_thresh = cl_args["variable_thresh"],
                                      variable_count = cl_args["variable_count"])
        
    # Multiprocess to handle wells
    with Pool(cl_args["jobs"]) as p:
        processed_well_results = list(tqdm(
            p.imap_unordered(complete_multiprocessor, all_wells),
            desc="Processing wells:", total=len(all_wells))
        )  
    
    # Handle processed output. This saves the summary dataframes, generates 
    # platemaps, and saves alignments (if requested)
    print("Saving outputs to disk...")
    format_and_save_outputs(processed_well_results,
                            cl_args["output"],
                            cl_args["return_alignments"])
