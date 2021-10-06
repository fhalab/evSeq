"""
Contains the functions needed for running and evaluating stress tests.
"""
# Import evSeq stuff
import tests.data_generation.globals as test_glob
from .globals import SAVELOC
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

def run_aa_stress_test(expected_out, true_out):
    
    # Get the unique plates and wells between the two. 
    expected_platewell = {tuple(platewell) for platewell in 
                          expected_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    true_platewell = {tuple(platewell) for platewell in 
                      true_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    
    # Find the differences between the true and expected. These are our first
    # bad plate-well combos.
    bad_platewells = list(expected_platewell ^ true_platewell)
    good_platewells = list(expected_platewell & true_platewell)
    
    # Create a report for the bad platewells
    error_reports = ["Missing well" for _ in range(len(bad_platewells))]
    
    # For all potentially good platewells, test to be sure that we see the same
    # rows coming out of each
    for plate, well in good_platewells:
        
        # Get the limited dataframe
        limited_expected_out = expected_out.loc[(expected_out.IndexPlate == plate)&
                                                (expected_out.Well == well)]
        limited_true_out = true_out.loc[(true_out.IndexPlate == plate)&
                                        (true_out.Well == well)]
        
        # Assert the dataframes are the same size
        len_expected_out = len(limited_expected_out)
        len_true_out = len(limited_true_out)
        assert (len_expected_out > 0) and (len_true_out > 0)
        if len_expected_out != len_true_out:
            error_reports.append((limited_expected_out.copy(), limited_true_out.copy()))
            bad_platewells.append((plate, well))
            continue
            
        # Go row by row and make sure the two dataframes align
        for expected_row, true_row in itertools.zip_longest(limited_expected_out.itertuples(index = False),
                                                            limited_true_out.itertuples(index = False)):

            # Convert to dicts
            expected_row = expected_row._asdict()
            true_row = true_row._asdict()

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

            # Report if the rows are not equivalent
            if not row_passes:
                bad_platewells.append((plate, well))
                error_reports.append([limited_expected_out.copy(),
                                      limited_true_out.copy()])
                break
                
    # Evaluate how many bad platewells were found.
    successful_test = (len(bad_platewells) == 0)
    return successful_test, bad_platewells, error_reports

def build_new_dir(new_dir):
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)

def compare_datasets(expected_out, true_out, flavor, counter):
    
    # Test the two dataframes from each set to make sure they agree
    test_passed, bad_platewells, reports = run_aa_stress_test(expected_out, true_out)
    
    # Note success if all tests passed
    if test_passed:
        print(f"All {flavor} tests passed for seed: {counter}")

    # Save the error reports if there were any
    else:
        # Report errors
        for plate, well in bad_platewells:
            print(f"Errors found for {flavor} {plate}-{well} with seed {counter}.")

        # Save the messed up components
        error_loc = os.path.join(SAVELOC, "ErrorReports")
        build_new_dir(error_loc)

        error_loc = os.path.join(error_loc, flavor)
        build_new_dir(error_loc)

        with open(os.path.join(error_loc, f"{counter}.pkl"), "wb") as f:
            pickle.dump([bad_platewells, reports], f)
            
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
        test_run.run_evseq()

        # Get the expected outputs
        expected_decoupled, expected_coupled = test_run.build_expected_aa()

        # Get the true outputs. Sort the true output in the same
        # way the expected was sorted.
        most_recent_run_path = sorted(glob(os.path.join(SAVELOC, "evSeqOutput", "*")))[-1]
        true_decoupled = pd.read_csv(os.path.join(most_recent_run_path, "OutputCounts", 
                                                  "AminoAcids_Decoupled_All.csv"))
        true_coupled = pd.read_csv(os.path.join(most_recent_run_path, "OutputCounts",
                                                "AminoAcids_Coupled_All.csv"))
        
        true_decoupled.sort_values(by = ["IndexPlate", "Well", "AaPosition", "Aa"],
                             inplace = True)
        true_coupled.sort_values(by = ["IndexPlate", "Well", "AlignmentFrequency", "SimpleCombo"],
                                 inplace = True)
        
        # Test the two dataframes from each set to make sure they agree
        uncoupled_passed = compare_datasets(expected_decoupled,
                                            true_decoupled,
                                            "Uncoupled",
                                            counter)
        coupled_passed = compare_datasets(expected_coupled,
                                          true_coupled,
                                          "Coupled",
                                          counter)
        
        # If both tests passed, delete output
        if uncoupled_passed and coupled_passed and not keep_output:
            shutil.rmtree(most_recent_run_path)
            
        # Update the counter and print a line break for the next test
        counter += 1
        print("\n")