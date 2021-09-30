# Import evSeq objects
from .util.logging import log_info, log_warning
from .util.input_validation import check_args
from .util.input_processing import load_all
from .seq_pair import SeqPair
from .well import Well
from .data_visualization import (generate_qualplot, 
                                 generate_platemaps,
                                 save_platemap_to_file)

# Import other required modules
import os
import numpy as np
import pandas as pd
import scipy.stats as ss
from Bio import SeqIO
from itertools import chain
from functools import partial
from multiprocessing import Pool
import tqdm
import gzip

def build_seqpairs(f_loc, r_loc, tqdm_fn=tqdm.tqdm):
    """Load and pair fastq file entries."""

    # Create a dictionary that links id to sequence object
    id_to_reads = {}

    # Determine how to open files if zipped or not
    f_opener, f_encoding = (gzip.open, 'rt') if "fastq.gz" in f_loc \
                      else (open, 'r')
    r_opener, r_encoding = (gzip.open, 'rt') if "fastq.gz" in r_loc \
                      else (open, 'r')

    # Store reads as list in memory
    print("Loading forward reads...")
    with f_opener(f_loc, f_encoding) as handle:
        all_f_recs = list(SeqIO.parse(handle, "fastq"))

    # For Gooey (non-tqdm) progress
    simple = True if tqdm_fn.__name__ == 'blank' else False
    
    # Parse the forward reads
    if simple:
        print('Parsing forward reads...')
        percents = list(range(0, 101))
    for i, f_record in enumerate(tqdm_fn(all_f_recs,
                                         desc="Parsing forward reads..."), 1):
        temp_record = SeqPair()
        temp_record.assign_f(f_record)
        assert f_record.id not in id_to_reads, "Duplicate IDs detected"
        id_to_reads[f_record.id] = temp_record

        # For simple Gooey progress
        if simple:
            percent = 100*(i / len(all_f_recs))
            if int(percent) == percents[0]:
                progress = percents.pop(0)
                print(f"Progress: {progress}%")
    
    # Associate reverse reads with the forward
    print("Loading reverse reads...")
    with r_opener(r_loc, r_encoding) as handle:
        all_r_recs = list(SeqIO.parse(handle, "fastq"))
    
    if simple:
        print('Pairing reverse reads...')
        percents = list(range(0, 101))
    observed_reads = set()
    for i, r_record in enumerate(tqdm_fn(all_r_recs,
                                         desc="Pairing reverse reads..."), 1):

        # If there is no partern in id_to_reads, create a new object 
        # and continue
        if r_record.id not in id_to_reads:
            temp_record = SeqPair()
            temp_record.assign_r(r_record)
            id_to_reads[r_record.id] = temp_record

        # Otherwise, attach the reverse record
        else:
            assert r_record.id not in observed_reads, "Duplicate ID found"
            observed_reads.add(r_record.id)
            id_to_reads[r_record.id].assign_r(r_record)
            
        # For simple Gooey progress
        if simple:
            percent = 100*(i / len(all_r_recs))
            if int(percent) == percents[0]:
                progress = percents.pop(0)
                print(f"Progress: {progress}%")
            
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
    match = None,
    mismatch = None,
    open_penalty = None,
    extend = None
):
    """Processes a single well."""
    
    # Define a frequency warning variable
    freq_warning = ""
    
    # Align
    well.align({"match": match,
                "mismatch": mismatch,
                "open_penalty": open_penalty,
                "extend": extend})

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
    heatmaps = generate_platemaps(max_outs[-1])
    save_platemap_to_file(heatmaps, saveloc)

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


def run_evSeq(cl_args, tqdm_fn=tqdm.tqdm):
    """Runs the evSeq protocol. Should be run from command line."""
    # Check the input arguments
    check_args(cl_args)
    
    # Identify sequencing files and load reference sequence information
    forward_file, reverse_file, bc_to_ref_plate_well = load_all(cl_args)
    
    # Pair all sequences
    all_seqpairs = build_seqpairs(forward_file, reverse_file, tqdm_fn)
    
    # Generate quality plots
    generate_qualplot(all_seqpairs, cl_args["output"])
    
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
                                      variable_count = cl_args["variable_count"],
                                      match = cl_args["match_score"],
                                      mismatch = -cl_args["mismatch_penalty"],
                                      open_penalty = -cl_args["gap_open_penalty"],
                                      extend = -cl_args["gap_extension_penalty"])
        
    # Multiprocess to handle wells
    with Pool(cl_args["jobs"]) as p:
        iterator = tqdm_fn(
            p.imap_unordered(complete_multiprocessor, all_wells),
            desc="Processing wells...",
            total=len(all_wells)
        )
        # For Gooey (non-tqdm) progress
        simple = True if tqdm_fn.__name__ == 'blank' else False
        if simple:
            processed_well_results = [_ for _ in all_wells]
            for i, processed_well_result in enumerate(iterator, 1):
                processed_well_results[i] = processed_well_result
                percent = int(100*(i / len(all_wells)))
                print(f"Progress: {percent}%")
        else:
            processed_well_results = list(iterator)
    
    # Handle processed output. This saves the summary dataframes, generates 
    # platemaps, and saves alignments (if requested)
    print("Saving outputs to disk...")
    format_and_save_outputs(processed_well_results,
                            cl_args["output"],
                            cl_args["return_alignments"])
