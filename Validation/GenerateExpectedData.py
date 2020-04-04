#!/usr/bin/env python
# coding: utf-8

# Import required modules
from glob import glob
from multiprocessing import Pool
import os.path
import re
import pandas as pd
import numpy as np
import csv
import pickle

# Define sets of positions that struggle with insertions or deletions in the last position
indel_struggles = {"./TestDatasets/DifferentOverlapDegrees/Set4",
                   "./TestDatasets/DifferentRefSeqByWell/"}

################## Functions to Build Expected Summary Information #####################

# Write a function to build expected summary files for different overlap degrees
def BuildExpectedSummary(seq_info_dir, r1_loc, r2_loc, aa_info_dict, all_aas,
                         different_AAs = False, different_Qs = False, different_by_well = False):

    # Create seq_info_path
    seq_info_path = os.path.join(seq_info_dir, "SequenceInfo.csv")
    
    # Open the SequenceInfo.csv file
    sequence_info_df = pd.read_csv(seq_info_path)

    # Create a column that combines plate and well
    sequence_info_df["PlateWell"] = sequence_info_df["IndexPlate"] + sequence_info_df["Well"]

    # Get unique plates
    unique_plates = sequence_info_df.IndexPlate.unique()

    # Get unique plate-well combinations
    plate_well_combos = sequence_info_df.PlateWell.unique()

    # Create a list to store expected info for each summary data type
    summary_cols = ("Plate", "Well", "ReadDirection", "F-BC", "R-BC", "Site", "AA", "AlignmentFrequency", "WellSeqDepth", "R1", "R2")
    variant_cols = ("Plate", "Well", "F-BC", "R-BC", "VariantCombo", "AlignmentFrequency", "WellSeqDepth", "R1", "R2")
    summary_info_list = []
    variant_info_list = []

    # Get the forward and reverse AAs as well as their q-scores if we are using the same ones each time
    if not different_by_well:
        faas = aa_info_dict["FAAs"]
        fqs = aa_info_dict["FQs"]
        raas = aa_info_dict["RAAs"]
        rqs = aa_info_dict["RQs"]

        # Get the shared aa info
        shared_aas = set(aa_info_dict["OverlapAAs"])
      
    # Create a variable to store the insertion-deletion checks in. Reads that carry both an insertion
    # and a deletion will be skipped in our tests.
    skips = []
        
    # Loop over each plate-well combo
    n_total = len(plate_well_combos)
    for test_ind, platewell in enumerate(plate_well_combos):
        
        # If we have a different ref seq by each well, then we need to construct a different set of
        # faas, fqs, raas, rqs, and shared_aas for each platewell combo
        if different_by_well:
            
            # Identify the correct dictionary for the platewell index
            faas = aa_info_dict[platewell]["FAAs"]
            fqs = aa_info_dict[platewell]["FQs"]
            raas = aa_info_dict[platewell]["RAAs"]
            rqs = aa_info_dict[platewell]["RQs"]

            # Get the shared aa info
            shared_aas = set(aa_info_dict[platewell]["OverlapAAs"])
        
        # Pull the limited dataframe
        limited_df = sequence_info_df.loc[sequence_info_df.PlateWell == platewell, :]

        # Define a dictionary to store counts at each AA
        count_dict_f = {aa_ind: {} for aa_ind in faas}
        count_dict_r = {aa_ind: {} for aa_ind in raas}
        combo_dict = {}

        # Package aa_sets and count_dicts
        iterables = (("F", faas, fqs, count_dict_f),
                    ("R", raas, rqs, count_dict_r))
        
        # Loop over the limited_df row by row
        for ind, row in limited_df.iterrows():
                
            # Loop over each set of amino acid identities and qualities
            for read_ind, (direction, aa_set, q_set, count_dict) in enumerate(iterables):
                
                # If a read has both insertions and deletions in it, skip
                all_read_aas = set(row[aa] for aa in aa_set)
                if "Ins" in all_read_aas and "Del" in all_read_aas:
                    skips.append(platewell)
                                    
                # If a read has either an insertion or deletion, don't record info
                if "Ins" in all_read_aas or "Del" in all_read_aas:
                    
                    # If this is position 4 of set 4 of different overlap degrees and an insertion or deletion, skip
                    if seq_info_dir in indel_struggles and aa == "AA4":
                        skips.append(platewell)
                        
                    continue
                
                # If we failed the alignment filter, discard the run
                if direction == "F" and row["NoiseF"]:
                    continue
                if direction == "R" and row["NoiseR"]:
                    continue
                
                # Loop over each aa and its quality
                for aa, q in zip(aa_set, q_set):
                    
                    # Check to make sure we are not an insertion or deletion
                    insdel_check = (row[aa] != "Ins" and row[aa] != "Del") 
                    
                    # If we meet the Q-threshold and are not an insertion or deletion
                    if row[q] and insdel_check:

                        # If we have not recorded this AA previously, add to the dictionary
                        if row[aa] not in count_dict[aa]:
                            count_dict[aa][row[aa]] = 1

                        # Otherwise, increment
                        else:
                            count_dict[aa][row[aa]] += 1

            # Loop over all amino acids
            combo = ""
            record_combo = True
            for aa_ind, aa in enumerate(all_aas, 1):
                
                # Determine the forward aa and reverse aa and their qualities
                faa = aa
                fq = "Q" + str(aa_ind)
                raa = aa + "_2" if different_AAs else aa
                rq = fq + "_2" if different_Qs else fq
                                
                # Check to make sure we don't have a noisy sequence
                f_noise_check = row["NoiseF"]
                r_noise_check = row["NoiseR"]
                
                # If either sequence is noisy, break the loop
                if f_noise_check or r_noise_check:
                    record_combo = False
                    break
                
                # Make sure we are not an insertion or deletion. If we are, break
                f_insdel_check = (row[faa] != "Ins" and row[faa] != "Del")
                r_insdel_check = (row[faa] != "Ins" and row[faa] != "Del")
                if not (f_insdel_check and r_insdel_check):
                    record_combo = False
                    break
                
                # If this is not in shared aas but is in forward, check forward quality only
                if aa not in shared_aas and faa in faas:

                    # If we meet the Q-threshold, record
                    if row[fq]:
                        combo += row[faa]

                    # Otherwise, record that we should not record the combo and break
                    else:
                        record_combo = False
                        break

                # If this is not in shared aas but is in reverse, check reverse quality only
                elif aa not in shared_aas and raa in raas:

                    # If we meet the Q-threshold, record
                    if row[rq]:
                        combo += row[raa]

                    # Otherwise record that we should not record the combo and break
                    else:
                        record_combo = False
                        break

                # If this is in shared, check to be sure that both qualities and both aas match one another
                elif aa in shared_aas:

                    # Make sure that both AAs are the same
                    # Note that this might need to be changed in the future!
                    if row[fq] and row[rq] and row[faa] == row[raa]:
                        combo += row[aa]

                    # Otherwise, record taht we should not record the combo and break
                    else:
                        record_combo = False
                        break

            # If we are recording the combo, do this now
            if record_combo:

                # If we have not seen the combo in the dictionary, add to the dictionary
                if combo not in combo_dict:
                    combo_dict[combo] = 1

                # Otherwise, increment what is already there
                else:
                    combo_dict[combo] += 1

        # Extract plate, well, and bc info 
        index_plate = limited_df.IndexPlate.values[0]
        index_well = limited_df.Well.values[0]
        fbc = limited_df["F-BC"].values[0]
        rbc = limited_df["R-BC"].values[0]
        
        # Package aas
        iterables = (("Forward", faas, count_dict_f),
                    ("Reverse", raas, count_dict_r))
        
        # Construct summary dataframe by looping over forward and reverse aas
        for direction, aa_set, count_dict in iterables:

            # Construct the rows for the summary dataframe from forward reads
            for site, aa in enumerate(aa_set, 1):

                # Get the total count of reads for the well at the position
                total_read_count = sum(val for key, val in count_dict[aa].items())

                # Record all information for the specific amino acid identified
                summary_info = [[index_plate, index_well, direction, fbc, rbc, site, 
                                aa, count/total_read_count, total_read_count, r1_loc, r2_loc]
                                for aa, count in count_dict[aa].items()]

                # Extend summary info list
                summary_info_list.extend(summary_info)

        # Construct the rows for the combo data
        total_combo_count = sum(val for key, val in combo_dict.items())
        combo_info = [[index_plate, index_well, fbc, rbc, combo, count/total_combo_count,
                       total_combo_count, r1_loc, r2_loc ]for combo, count in combo_dict.items()]

        # Extend the variant info list
        variant_info_list.extend(combo_info)
        
    # Create a skips file for the entire directory
    skips_loc = os.path.join(seq_info_dir, "Skips.pkl")
    with open(skips_loc, "wb") as f:
        pickle.dump(set(skips), f)

    # Once everything is finished, convert each list to a dataframe
    summary_df = pd.DataFrame(summary_info_list, columns = summary_cols)
    variant_df = pd.DataFrame(variant_info_list, columns = variant_cols)

    # Loop over the unique plates and create a dataframe for each plate
    n_unique_plates = len(unique_plates)
    by_plate_summary_dfs = [None for _ in range(n_unique_plates)]
    by_plate_variant_dfs = [None for _ in range(n_unique_plates)]
    by_plate_max_dfs = [None for _ in range(n_unique_plates)]
    for i, unique_plate in enumerate(unique_plates):

        # Pull the part of the dataframe corresponding to the unique plate
        by_plate_summary_dfs[i] = summary_df.loc[summary_df.Plate == unique_plate, :].copy()
        by_plate_variant_dfs[i] = variant_df.loc[variant_df.Plate == unique_plate, :].copy()

        # Build a max plate
        max_alignment_complete = by_plate_variant_dfs[i].copy()
        group_plate = max_alignment_complete.groupby(by = "Well")[["Well", "AlignmentFrequency"]].max().reset_index(drop=True)
        
        try:
            max_alignment_complete = max_alignment_complete.merge(group_plate)
        except:
            by_plate_max_dfs[i] = max_alignment_complete
            continue

        # Construct a logseqdepth column
        max_alignment_complete["logseqdepth"] = np.log(max_alignment_complete["WellSeqDepth"])
        
        # Record the dataframe
        by_plate_max_dfs[i] = max_alignment_complete.copy()
        
    # Return all dataframes
    return by_plate_summary_dfs, by_plate_variant_dfs, by_plate_max_dfs

# Write a function that builds the expected summary information for each different overlap set
def BuildExpectedSummaryDifOverlap():
    
    # Define the setnames
    setnames = ("Set1", "Set2", "Set3", "Set4")
    
    # Define all_aas
    all_aas = ("AA1", "AA2", "AA3", "AA4")
    
    # Define the locations of all sets
    set_locations = [os.path.join("./TestDatasets/DifferentOverlapDegrees", setname)
                    for setname in setnames]
    
    # Get the locations of all test data
    test_data_locations_R1 = [os.path.join(set_loc, "TestData_R1_test.fastq")
                              for set_loc in set_locations]
    test_data_locations_R2 = [os.path.join(set_loc, "TestData_R2_test.fastq")
                              for set_loc in set_locations]
    
    # Define the aa_info_dicts
    info_by_set = {"Set1": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                            "RAAs": ("AA3", "AA4"),
                            "RQs": ("Q3_2", "Q4_2"),
                            "OverlapAAs": ("AA3",)},
                   "Set2": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                           "RAAs": ("AA2", "AA3", "AA4"),
                            "RQs": ("Q2_2", "Q3_2", "Q4_2"),
                            "OverlapAAs": ("AA2", "AA3")},
                   "Set3": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                           "RAAs": ("AA1", "AA2", "AA3", "AA4"),
                            "RQs": ("Q1_2", "Q2_2", "Q3_2", "Q4_2"),
                            "OverlapAAs": ("AA1", "AA2", "AA3")},
                   "Set4": {"FAAs": ("AA1", "AA2", "AA3", "AA4"),
                            "FQs": ("Q1", "Q2", "Q3", "Q4"),
                           "RAAs": ("AA1", "AA2", "AA3", "AA4"),
                            "RQs": ("Q1_2", "Q2_2", "Q3_2", "Q4_2"),
                            "OverlapAAs": ("AA1", "AA2", "AA3", "AA4")}
                  }
    
    # Generate all expected outputs
    expected_outputs = [BuildExpectedSummary(set_locations[i],
                                            test_data_locations_R1[i],
                                            test_data_locations_R2[i],
                                            info_by_set[setname],
                                            all_aas,
                                            different_Qs = True)
                       for i, setname in enumerate(setnames)]
    
    # Return all expected outputs
    return expected_outputs

# Write a function that builds the expected summary information for each different combinatorial sites
def BuildExpectedSummaryDifSites():
    
    # Define the setnames
    setnames = ("Pos1", "Pos2", "Pos3")
    
    # Define all_aas
    all_aas = ("AA1", "AA2", "AA3", "AA4")
    
    # Define the locations of all sets
    set_locations = [os.path.join("./TestDatasets/DifferentPositions", setname)
                    for setname in setnames]
    
    # Get the locations of all test data
    test_data_locations_R1 = [os.path.join(set_loc, "TestData_R1_test.fastq")
                              for set_loc in set_locations]
    test_data_locations_R2 = [os.path.join(set_loc, "TestData_R2_test.fastq")
                              for set_loc in set_locations]
    
    # Define the aa_info_dicts
    info_by_set = {"Pos1": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                            "RAAs": ("AA4",),
                            "RQs": ("Q4",),
                            "OverlapAAs": ()},
                   "Pos2": {"FAAs": ("AA1", "AA2"),
                            "FQs": ("Q1", "Q2"),
                           "RAAs": ("AA3", "AA4"),
                            "RQs": ("Q3", "Q4"),
                            "OverlapAAs": ()},
                   "Pos3": {"FAAs": ("AA1",),
                            "FQs": ("Q1",),
                           "RAAs": ("AA2", "AA3", "AA4"),
                            "RQs": ("Q2", "Q3", "Q4"),
                            "OverlapAAs": ()}
                  }
    
    # Generate all expected outputs
    expected_outputs = [BuildExpectedSummary(set_locations[i],
                                            test_data_locations_R1[i],
                                            test_data_locations_R2[i],
                                            info_by_set[setname],
                                            all_aas)
                       for i, setname in enumerate(setnames)]
    
    # Return all expected outputs
    return expected_outputs

# Write a function that builds the expected summary information for when forward and reverse disagree
def BuildExpectedSummaryFRDisagree():
    
    # Define the setnames
    setnames = ("Set1", "Set2")
    
    # Define all_aas
    all_aas = ("AA1", "AA2", "AA3", "AA4")
    
    # Define the locations of all sets
    set_locations = [os.path.join("./TestDatasets/ForwardAndReverseDisagree", setname)
                    for setname in setnames]
    
    # Get the locations of all test data
    test_data_locations_R1 = [os.path.join(set_loc, "TestData_R1_test.fastq")
                              for set_loc in set_locations]
    test_data_locations_R2 = [os.path.join(set_loc, "TestData_R2_test.fastq")
                              for set_loc in set_locations]
    
    # Define the aa_info_dicts
    info_by_set = {"Set1": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                            "RAAs": ("AA3_2", "AA4_2"),
                            "RQs": ("Q3", "Q4"),
                            "OverlapAAs": ("AA3",)},
                   "Set2": {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                           "RAAs": ("AA2_2", "AA3_2", "AA4_2"),
                            "RQs": ("Q2", "Q3", "Q4"),
                            "OverlapAAs": ("AA2", "AA3")}
                  }
    
    # Generate all expected outputs
    expected_outputs = [BuildExpectedSummary(set_locations[i],
                                            test_data_locations_R1[i],
                                            test_data_locations_R2[i],
                                            info_by_set[setname],
                                            all_aas,
                                            different_AAs = True)
                       for i, setname in enumerate(setnames)]
    
    # Return all expected outputs
    return expected_outputs

# Write a function that builds the expected summary information for when the well is mixed
def BuildExpectedSummaryMixedWells():
    
    # Define all_aas
    all_aas = ("AA1", "AA2", "AA3", "AA4")
    
    # Get the locations of all sequence info
    set_location = "./TestDatasets/MixedWells/"
                         
    # Get the locations of all test data
    test_data_location_R1 = "./TestDatasets/MixedWells/TestData_R1_test.fastq"
    test_data_location_R2 = "./TestDatasets/MixedWells/TestData_R2_test.fastq"      
        
    # Define the aa_info_dicts
    aa_info_dict = {"FAAs": ("AA1", "AA2", "AA3"),
                            "FQs": ("Q1", "Q2", "Q3"),
                           "RAAs": ("AA1", "AA2", "AA3", "AA4"),
                            "RQs": ("Q1", "Q2", "Q3", "Q4"),
                            "OverlapAAs": ("A1", "AA2", "AA3")}
    
    # Generate all expected outputs
    expected_outputs = BuildExpectedSummary(set_location,
                                            test_data_location_R1,
                                            test_data_location_R2,
                                            aa_info_dict,
                                            all_aas)
                       
    # Return all expected outputs
    return expected_outputs

# Write a function that builds the expected summary information for using a single site
def BuildExpectedSummarySingleSite():
    
    # Define all_aas
    all_aas = ("AA1",)
    
    # Get the locations of all sequence info
    set_location = "./TestDatasets/SingleSite/"
                         
    # Get the locations of all test data
    test_data_location_R1 = "./TestDatasets/SingleSite/TestData_R1_test.fastq"
    test_data_location_R2 = "./TestDatasets/SingleSite/TestData_R2_test.fastq"
                              
    # Define the aa_info_dicts
    aa_info_dict = {"FAAs": ("AA1",),
                            "FQs": ("Q1",),
                           "RAAs": ("AA1",),
                            "RQs": ("Q1",),
                            "OverlapAAs": ("AA1",)}
    
    # Generate all expected outputs
    expected_outputs = BuildExpectedSummary(set_location,
                                            test_data_location_R1,
                                            test_data_location_R2,
                                            aa_info_dict,
                                            all_aas)
                       
    # Return all expected outputs
    return expected_outputs

# Write a function that builds the expected summary information for when we have 
# different reference sequences in each well
def BuildExpectedSummaryDifRefSeqs():
    
    # Define all_aas
    all_aas = ("AA1", "AA2", "AA3", "AA4")
    
    # Get the locations of all sequence info
    set_location = "./TestDatasets/DifferentRefSeqByWell/"
                         
    # Get the locations of all test data
    test_data_location_R1 = "./TestDatasets/DifferentRefSeqByWell/TestData_R1_test.fastq"
    test_data_location_R2 = "./TestDatasets/DifferentRefSeqByWell/TestData_R2_test.fastq"
                              
    # Load the aa_info_dicts. These are different by plate and well.
    with open("./TestDatasets/DifferentRefSeqByWell/RefSeqInfoDict.pkl", "rb") as f:
        aa_info_dict = pickle.load(f)
    
    # Generate all expected outputs
    expected_outputs = BuildExpectedSummary(set_location,
                                            test_data_location_R1,
                                            test_data_location_R2,
                                            aa_info_dict,
                                            all_aas,
                                           different_by_well = True)
                       
    # Return all expected outputs
    return expected_outputs

######################### Execute Construction #################################

# Define a list of tasks for building the expected data
summary_tasks = {"DifOverlap": BuildExpectedSummaryDifOverlap,
                "DifSites": BuildExpectedSummaryDifSites,
                "FRDisagree": BuildExpectedSummaryFRDisagree,
                "MixedWells": BuildExpectedSummaryMixedWells,
                "SingleSite": BuildExpectedSummarySingleSite,
                "DifRefSeqs": BuildExpectedSummaryDifRefSeqs}

# Define a set of storage locations for the expected output
storage_locs = {"DifOverlap": "./TestDatasets/DifferentOverlapDegrees/",
                "DifSites": "./TestDatasets/DifferentPositions/",
                "FRDisagree": "./TestDatasets/ForwardAndReverseDisagree/",
                "MixedWells": "./TestDatasets/MixedWells/",
                "SingleSite": "./TestDatasets/SingleSite/",
                "DifRefSeqs": "./TestDatasets/DifferentRefSeqByWell/"}

# Define single-output functions
single_output_funcs = {"DifRefSeqs", "SingleSite", "MixedWells"}

# Define summary file names
sum_filenames = ("ExpectedSummaryInfo.csv", "ExpectedVariantInfo.csv", "ExpectedMaxInfo.csv")

# Define a function to save outputs
def SaveExpectedSummaries(save_directory, summaries):
    
    # Save each output summary
    for i, summary_list in enumerate(zip(*summaries), 1):
        for j, summary_df in enumerate(summary_list):
        
            # Create the savename
            savename = os.path.join(save_directory, "Expected_Output",
                                    f"Test-DI0{i}-0_{sum_filenames[j]}")
            
            # Save the dataframe
            summary_df.to_csv(savename, index=False)

# Write a function that heps multiprocessing in BuildExpectedSummaryAll
def BuildExpectedSummaryHelper(task):
    
    # Pull the function from the list of tasks and execute
    expected_outputs = summary_tasks[task]()
    
    # Save the output
    if task in single_output_funcs:
        SaveExpectedSummaries(storage_locs[task], expected_outputs)
    else:
        for i, pos_output in enumerate(expected_outputs, 1):
            
            # If we are looking at Different, use "Pos"
            if task == "DifSites":
                save_dir = os.path.join(storage_locs[task], f"Pos{i}")
            else:
                save_dir = os.path.join(storage_locs[task], f"Set{i}")
                
            # Save the output
            SaveExpectedSummaries(save_dir, pos_output)    

# Write a function that constructs all expected dataframes for all input synthetic data.
def BuildExpectedSummaryAll(njobs = 1):
        
    # Make an empty Skips.pkl file for the real data
    with open("./TestDatasets/RealData/Skips.pkl", "wb") as f:
        pickle.dump(set(), f)
        
    # Package tasks
    tasks = list(storage_locs.keys())

    # Multiprocess
    with Pool(njobs) as p:
        _ = list(p.imap(BuildExpectedSummaryHelper, tasks))
        
# Build the expected summaries
if __name__ == "__main__":
    BuildExpectedSummaryAll()