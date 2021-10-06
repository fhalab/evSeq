"""
Contains the functions needed for running and evaluating stress tests.
"""
# Import evSeq stuff
import tests.data_generation.globals as test_glob
from .globals import COUPLED_SORT_COLS, DECOUPLED_SORT_COLS, SAVELOC
from .run_generator import FakeRun

# Import 3rd party code
import os
import random
import itertools
import pickle
import shutil
import numpy as np
import pandas as pd
from glob import glob

def make_aa_checks(true_row, expected_row):
    
    # Check if the rows are equal.
    row_passes = True
    for key, true_val in true_row.items():

        # Get the expected value
        expected_val = expected_row[key]

        # Special cases 
        if key == "AlignmentFrequency":
            is_equiv = (np.isclose(true_val, expected_val))
        elif key == "Flags" and (true_val is np.nan):
            is_equiv = (true_val is expected_val)
        elif key == "AaPosition":
            is_equiv = (str(true_val) == str(expected_val))
        else:
            is_equiv = (true_val == expected_val)

        # Break the loop if the true and expected values are not equivalent
        if not is_equiv:
            row_passes = False
            break
        
    return row_passes    

def evaluate_max_position(true_out_by_pos, expected_out_by_pos):
    
    # Convert the maximum expected dataframe to a dictionary
    maximum_row_true = [row._asdict() for row in 
                        true_out_by_pos.itertuples(index = False)]
    assert len(maximum_row_true) == 1
    
    # Get the set of possible maxima from the expected dataframe
    max_expected_rows = []
    max_freq_observed = 0.0
    for expected_row in expected_out_by_pos.itertuples(index = False):

        # Convert to a dictionary
        expected_row = expected_row._asdict()

        # Record the maxima seen so far in the expected rows
        if expected_row["AlignmentFrequency"] >= max_freq_observed:
            max_freq_observed = expected_row["AlignmentFrequency"]
            max_expected_rows.append(expected_row)
            
    # Confirm that at least one of the expected maximum rows matches the
    # true maximum
    return any(make_aa_checks(maximum_row_true[0], expected_maximum_row)
               for expected_maximum_row in max_expected_rows)

def check_max_df(limited_true_out_max, limited_expected_out, is_coupled):
    
    # If this is not coupled, we loop over all positions and apply the analysis
    if not is_coupled:
        
        # Get all unique positions. Confirm that the dataframes have matching
        # positions. This well failed if the positions do not match, and we return
        # `False`
        unique_positions = limited_expected_out.AaPosition.unique().tolist()
        if set(unique_positions) != set(limited_true_out_max.AaPosition.unique().tolist()):
            return False
        
        # Loop over unique positions
        for pos in unique_positions:
            
            # Grab the expected and true for this position
            expected_out_by_pos = limited_expected_out.loc[limited_expected_out.AaPosition == pos]
            true_out_by_pos = limited_true_out_max.loc[limited_true_out_max.AaPosition == pos]
            
            # Evaluate the dataframes
            max_matches = evaluate_max_position(true_out_by_pos, expected_out_by_pos)
            
            # Break the loop if the maximum ever does not match
            if not max_matches:
                break
            
    # Otherwise, just test to be sure that the two dataframes match
    else:
        max_matches = evaluate_max_position(limited_true_out_max, limited_expected_out)
        
    return max_matches

def check_nonmax_df(limited_true_out, limited_expected_out):
    
    # Loop over the two dataframes and make checks
    for expected_row, true_row in itertools.zip_longest(limited_expected_out.itertuples(index = False),
                                                        limited_true_out.itertuples(index = False)):

        # Convert to dicts
        true_row = true_row._asdict()
        expected_row = expected_row._asdict()
                
        # Check if the rows are equals
        row_passes = make_aa_checks(true_row, expected_row)

        # Report if the rows are not equivalent
        if not row_passes:
            break
            
    # Report whether or not the row passed. If all passed, then we should see `True`
    return row_passes
    
def run_aa_stress_test(expected_out, true_out, true_out_max):
    
    # Get the unique plates and wells between the two. 
    expected_platewell = {tuple(platewell) for platewell in 
                          expected_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    true_platewell = {tuple(platewell) for platewell in 
                      true_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    
    # Confirm that the true plates and wells of the maximum match those of the non-maximum
    true_max_platewell = {tuple(platewell) for platewell in 
                          true_out_max.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    assert true_platewell == true_max_platewell
    
    # Find the differences between the true and expected. These are our first
    # bad plate-well combos.
    bad_platewells = list(expected_platewell ^ true_platewell)
    good_platewells = list(expected_platewell & true_platewell)
    
    # Create reports for the bad platewells
    error_reports = ["Missing well" for _ in range(len(bad_platewells))]
    max_bad_wells = bad_platewells.copy()
    max_error_reports = error_reports.copy()
    
    # Get the columns that we will use for sorting
    is_coupled = ("SimpleCombo" in true_out.columns)
    sort_list = (COUPLED_SORT_COLS if is_coupled else DECOUPLED_SORT_COLS)
    
    # For all potentially good platewells, test to be sure that we see the same
    # rows coming out of each
    for plate, well in good_platewells:
        
        # Get the limited dataframe
        limited_expected_out = expected_out.loc[(expected_out.IndexPlate == plate)&
                                                (expected_out.Well == well)]
        limited_true_out = true_out.loc[(true_out.IndexPlate == plate)&
                                        (true_out.Well == well)]
        limited_true_out_max = true_out_max.loc[(true_out_max.IndexPlate == plate)&
                                                (true_out_max.Well == well)]
        
        # Assert the non-max dataframes are the same size
        len_expected_out = len(limited_expected_out)
        len_true_out = len(limited_true_out)
        assert (len_expected_out > 0) and (len_true_out > 0)
        if len_expected_out != len_true_out:
            bad_dfs = (limited_expected_out.copy(), limited_true_out.copy())
            error_reports.append(bad_dfs)
            bad_platewells.append((plate, well))
            max_error_reports.append(bad_dfs)
            max_bad_wells.append((plate, well))
            continue
        
        # Sort the expected and true dataframes
        limited_expected_out = limited_expected_out.sort_values(by = sort_list, ignore_index = True)
        limited_true_out = limited_true_out.sort_values(by = sort_list, ignore_index = True)        
            
        # Evaluate to see if the maximum and non-maximum dfs are correct
        max_passed = check_max_df(limited_true_out_max, limited_expected_out, is_coupled)
        non_max_passed = check_nonmax_df(limited_true_out, limited_expected_out)
        
        # Report errors
        if not non_max_passed:
            bad_platewells.append((plate, well))
            error_reports.append([limited_true_out, limited_expected_out])
            
        if not max_passed:
            max_bad_wells.append((plate, well))
            max_error_reports.append([limited_true_out_max, limited_expected_out])
                
    # Evaluate how many bad platewells were found.
    successful_test = (len(bad_platewells) == 0) and (len(max_bad_wells) == 0)
    return successful_test, bad_platewells, error_reports, max_bad_wells, max_error_reports

def build_new_dir(new_dir):
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)

def compare_datasets(expected_out, true_out, true_out_max, flavor, counter):
    
    # Test the two dataframes from each set to make sure they agree
    (test_passed,
     bad_platewells,
     reports,
     max_bad_wells,
     max_error_reports) = run_aa_stress_test(expected_out, true_out, true_out_max)
    
    # Note success if all tests passed
    if test_passed:
        print(f"All {flavor} tests passed for seed: {counter}")

    # Save the error reports if there were any
    else:
        # Report errors
        for plate, well in bad_platewells:
            print(f"Errors found for non-max {flavor} {plate}-{well} with seed {counter}.")
        for plate, well in max_bad_wells:
            print(f"Errors found for max {flavor} {plate}-{well} with seed {counter}.")

        # Save the messed up components
        error_loc = os.path.join(SAVELOC, "ErrorReports")
        build_new_dir(error_loc)

        error_loc = os.path.join(error_loc, flavor)
        build_new_dir(error_loc)

        with open(os.path.join(error_loc, f"{counter}.pkl"), "wb") as f:
            pickle.dump([bad_platewells, reports, max_bad_wells, max_error_reports], f)
            
    return test_passed

def run_evseq_stress_test(detailed, include_nnn, 
                          keep_output = False, seed = 0):
    
    # Run until we break something
    counter = seed
    while True:
                
        # Report what we're working on
        print(f"Working on tests for seed {counter}...")
                
        # Update the global RNG to match the counter (for reproducbility)
        test_glob.RANDOM_SEED = counter
        test_glob.NP_RNG = np.random.default_rng(counter)
        test_glob.RANDOM_RNG = random.Random(counter)
    
        # Build a test run and the associated output files
        test_run = FakeRun(detailed = detailed)
        test_run.build_fastq()
        test_run.build_refseq(include_nnn)

        # Run evSeq on the generated data
        # test_run.run_evseq()

        # Get the expected outputs
        expected_decoupled, expected_coupled = test_run.build_expected_aa()

        # Get the true outputs. 
        most_recent_run_path = sorted(glob(os.path.join(SAVELOC, "evSeqOutput", "*")))[-1]
        count_path = os.path.join(most_recent_run_path, "OutputCounts")
        true_decoupled = pd.read_csv(os.path.join(count_path, "AminoAcids_Decoupled_All.csv"))
        true_coupled = pd.read_csv(os.path.join(count_path, "AminoAcids_Coupled_All.csv"))
        true_decoupled_max = pd.read_csv(os.path.join(count_path, "AminoAcids_Decoupled_Max.csv"))
        true_coupled_max = pd.read_csv(os.path.join(count_path, "AminoAcids_Coupled_Max.csv"))
        
        # Test the two dataframes from each set to make sure they agree
        uncoupled_passed = compare_datasets(expected_decoupled,
                                            true_decoupled,
                                            true_decoupled_max,
                                            "Uncoupled",
                                            counter)
        coupled_passed = compare_datasets(expected_coupled,
                                          true_coupled,
                                          true_coupled_max,
                                          "Coupled",
                                          counter)
        
        # If both tests passed, delete output
        if uncoupled_passed and coupled_passed and not keep_output:
            shutil.rmtree(most_recent_run_path)
            
        # Update the counter and print a line break for the next test
        counter += 1
        print("\n")