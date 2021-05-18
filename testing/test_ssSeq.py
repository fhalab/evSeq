# Import required modules
from glob import glob
import os.path
import re
import pandas as pd
import numpy as np
import csv
import pickle

# Write a function that selects base files based on a regex pattern
def SelectBaseFiles(search_dir, regex_pattern):
    
    # Compile regex
    compiled_regex = re.compile(regex_pattern)
    
    # Pull all files in the search dir and sort
    all_files = sorted(glob(os.path.join(search_dir, "*")))

    # Loop over the files and find matches to the regex pattern
    selected_files = [file for file in all_files if compiled_regex.search(file)]      
        
    # Return selected files
    return selected_files

# Write a function that loads all summary information
def LoadSummaryInfo(output_dir):
    
    # Get the summary directory
    summary_dir = os.path.join(output_dir, "Summaries")
    
    # Load MaxInfo, SummaryInfo, and VariantInfo
    max_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "MaxInfo.csv")]
    summary_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "SummaryInfo.csv")]
    variant_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "VariantInfo.csv")]
    
    # Return all three files
    return max_info_files, summary_info_files, variant_info_files

# Write a function that loads the data from the latest 2 ssSeq output folders
# The first folder will be standard data and the second folder will be troubleshoot data
def LoadLatestData(ssSeq_output):
    
    # Identify all runs in the folder and sort
    all_folders = sorted(glob(os.path.join(ssSeq_output, "*")))
    
    # Get the last 2 folders (most recent folders)
    most_recents = all_folders[-2:]
    
    # Load the different summaries for the non troubleshooting and troubleshooting data
    nonts_summaries = LoadSummaryInfo(most_recents[0])
    ts_summaries = LoadSummaryInfo(most_recents[1])
    
    # Return summary information
    return nonts_summaries, ts_summaries

# Write a function that loads expected summary info from a target directory
def LoadExpectedSummary(target_dir):
    
    # Get expected summary directory
    summary_dir = os.path.join(target_dir, "Expected_Output")
    
    # Load MaxInfo, SummaryInfo, and VariantInfo
    max_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "MaxInfo.csv")]
    summary_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "SummaryInfo.csv")]
    variant_info_files = [pd.read_csv(file) for file in SelectBaseFiles(summary_dir, "VariantInfo.csv")]
    
    # Return all files
    return max_info_files, summary_info_files, variant_info_files

# Write a function that modifies a dataframe to remove skipped plate-well combos
def SkipSkips(dataframe, skips):
    
    # Create a platewell column in the dataframe
    dataframe["PlateWell"] = dataframe["Plate"] + dataframe["Well"]
    
    # Filter the dataframe
    dataframe = dataframe.loc[~dataframe.PlateWell.isin(skips), :].copy()
    
    # Return the dataframe
    return dataframe

# Do the troubleshoot and non-troubleshoot summary files match?
def test_ts_match(summary_files, target_dir):
    
    # Loop over the summary groups
    for summary_groups in zip(*summary_files):
        
        # Loop over the pairs of dataframes
        for nonts_df, ts_df in zip(*summary_groups):
            
            # Make sure the two dataframes match
            assert nonts_df.equals(ts_df), f"Mismatch between non-ts and ts for {target_dir}"


###############################################################################
############################## Test Functions #################################

def test_max_info_match(ssSeq_generated_max_info, expected_max_info, target_dir, plate_name, skips):
    
    # Skip appropriate plate-wells in each dataframe
    ssSeq_generated_max_info = SkipSkips(ssSeq_generated_max_info, skips)
    expected_max_info = SkipSkips(expected_max_info, skips)
    
    # Sort each file
    ssSeq_generated_max_info.sort_values(by = ["Plate", "Well", "VariantCombo"], 
                                         inplace = True)
    expected_max_info.sort_values(by = ["Plate", "Well", "VariantCombo"], 
                                  inplace = True)
    
    # Reset index of each file
    ssSeq_generated_max_info.reset_index(drop = True, inplace = True)
    expected_max_info.reset_index(drop = True, inplace = True)
    
    # Round all numerics to 3 decimal points
    ssSeq_generated_max_info = ssSeq_generated_max_info.round(3)
    expected_max_info = expected_max_info.round(3)
    
    # Make sure the dataframes are equal
    assert ssSeq_generated_max_info.equals(expected_max_info), f"Mismatch between MaxInfo.csv expected and ssSeq output for {target_dir}, {plate_name}"

# Write a fucntion that confirms SummaryInfo.csv is correct
def test_summary_info_match(ssSeq_generated_summary_info, expected_summary_info, target_dir, plate_name, skips):
    
    # Skip appropriate plate-wells in each dataframe
    ssSeq_generated_summary_info = SkipSkips(ssSeq_generated_summary_info, skips)
    expected_summary_info = SkipSkips(expected_summary_info, skips)
    
    # Sort each file
    ssSeq_generated_summary_info.sort_values(by = ["Plate", "Well", "ReadDirection", "Site", "AlignmentFrequency", "AA"],
                                            inplace = True)
    expected_summary_info.sort_values(by = ["Plate", "Well", "ReadDirection", "Site", "AlignmentFrequency", "AA"],
                                      inplace = True)
    
    # Reset index of each file
    ssSeq_generated_summary_info.reset_index(drop = True, inplace = True)
    expected_summary_info.reset_index(drop = True, inplace = True)
    
    # Force appropriate columns to integer and float
    for df in (ssSeq_generated_summary_info, expected_summary_info):
        for column in ("Site", "WellSeqDepth"):
            df[column] = df[column].astype(int)
        df["AlignmentFrequency"] = df["AlignmentFrequency"].astype(float)
    
    # Round all numerics to 3 decimal places
    ssSeq_generated_summary_info = ssSeq_generated_summary_info.round(3)
    expected_summary_info = expected_summary_info.round(3)
    
    # Confirm that the dataframes are equal
    assert ssSeq_generated_summary_info.equals(expected_summary_info), f"Mismatch between SummaryInfo.csv expected and ssSeq output for {target_dir}, {plate_name}"

# Write a function that confirms that VariantInfo.csv is correct
def test_variant_info_match(ssSeq_generated_variant_info, expected_variant_info, target_dir, plate_name, skips):
    
    # Skip appropriate plate-wells in each dataframe
    ssSeq_generated_variant_info = SkipSkips(ssSeq_generated_variant_info, skips)
    expected_variant_info = SkipSkips(expected_variant_info, skips)
    
    # Sort each file
    ssSeq_generated_variant_info.sort_values(by = ["Plate", "Well", "AlignmentFrequency", "VariantCombo"],
                                            inplace = True)
    expected_variant_info.sort_values(by = ["Plate", "Well", "AlignmentFrequency", "VariantCombo"],
                                     inplace = True)
    
    # Reset index of each file
    ssSeq_generated_variant_info.reset_index(drop = True, inplace = True)
    expected_variant_info.reset_index(drop = True, inplace = True)
    
    # Round all numerics to 3 decimal places
    ssSeq_generated_variant_info = ssSeq_generated_variant_info.round(3)
    expected_variant_info = expected_variant_info.round(3)
    
    assert ssSeq_generated_variant_info.equals(expected_variant_info), f"Mismatch between VariantInfo.csv expected and ssSeq output for {target_dir}, {plate_name}"
    
# Define a function that reports passes
def ReportPass(target_dir, plate_name, test_type):
    
    print(f"""
          
===================================================================
Test on {target_dir}, {plate_name} {test_type} PASSED

          """)    
    
def ReportFail(e):
    
    print(f"""

===================================================================
FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL 
{e}

          """)

################################################################################
################################## RUN #########################################
# Define the target directories
target_dirs = ("./TestDatasets/DifferentOverlapDegrees/Set1/",
              "./TestDatasets/DifferentOverlapDegrees/Set2/",
              "./TestDatasets/DifferentOverlapDegrees/Set3/",
              "./TestDatasets/DifferentOverlapDegrees/Set4/",
               "./TestDatasets/DifferentPositions/Pos1/",
               "./TestDatasets/DifferentPositions/Pos2/",
               "./TestDatasets/DifferentPositions/Pos3/",
               "./TestDatasets/DifferentRefSeqByWell/",
               "./TestDatasets/ForwardAndReverseDisagree/Set1/",
               "./TestDatasets/ForwardAndReverseDisagree/Set2/",
               "./TestDatasets/MixedWells/",
               "./TestDatasets/SingleSite/",
              "./TestDatasets/RealData/")

# Write a script that loops over all directories containing testable information
for target_dir in target_dirs:
    
    # Get the location of ssSeq_output 
    ssseq_output_path = os.path.join(target_dir, "ssSeq_Output")
    
    # Load the latest data
    ssseq_output = LoadLatestData(ssseq_output_path)
    
    # Compare the troubleshoot and non-troubleshoot files
    test_ts_match(ssseq_output, target_dir)
    
    # Load the expected data
    expected_summaries = LoadExpectedSummary(target_dir)
    
    # Load the skips file
    with open(os.path.join(target_dir, "Skips.pkl"), "rb") as f:
        skips = pickle.load(f)

    # Now compare loop over the summary data side-by-side and compare
    for i, ((ssseq_max_dfs, ssseq_summ_dfs, ssseq_var_dfs), \
    (expect_max_dfs, expect_summ_dfs, expect_var_dfs))\
    in enumerate(zip(zip(*ssseq_output[0]), zip(*expected_summaries)), 1):

        # Get the plate name
        plate_name = f"DI0{i}"
        
        # Run the max test
        try:
            test_max_info_match(ssseq_max_dfs, expect_max_dfs, target_dir, plate_name, skips)
            ReportPass(target_dir, plate_name, "MaxInfo.csv")
        except Exception as e:
            ReportFail(e)
            input("Enter to continue")
            
        # Run the summary test
        try:
            test_summary_info_match(ssseq_summ_dfs, expect_summ_dfs, target_dir, plate_name, skips)
            ReportPass(target_dir, plate_name, "SummaryInfo.csv")
        except Exception as e:
            ReportFail(e)
            input("Enter to continue")
        
        # Run the variant test
        try:
            test_variant_info_match(ssseq_var_dfs, expect_var_dfs, target_dir, plate_name, skips)
            ReportPass(target_dir, plate_name, "VariantInfo.csv")
        except Exception as e:
            ReportFail(e)
            input("Enter to continue")