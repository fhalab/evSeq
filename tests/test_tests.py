"""
Yeah, this one's meta. This file contains quick unit tests to make sure that the
functions used to confirm that true evSeq output matches expected evSeq output
in the stress testing.
"""
# Import required evSeq machinery
from .data_generation.globals import SAVELOC
from .data_generation.stress_tests import test_aa

# Import 3rd party code
import os
import pandas as pd

# Write a test to make sure that the stress test checker actually works
def test_aa_checker():
    
    # Get the paths to the aa files
    sample1_path = os.path.join(SAVELOC, "sample_count_files", "sample1")
    sample2_path = os.path.join(SAVELOC, "sample_count_files", "sample2")
    
    # Load in the expected aa files
    testdf_1 = pd.read_csv(os.path.join(sample1_path, "AminoAcids_Decoupled_All.csv"))
    testdf_2 = pd.read_csv(os.path.join(sample2_path, "AminoAcids_Decoupled_All.csv"))

    # Get the number of unique plate-well combos
    test1_platewells = testdf_1.loc[:, ["IndexPlate", "Well"]].drop_duplicates()
    test2_platewells = testdf_2.loc[:, ["IndexPlate", "Well"]].drop_duplicates()
    combined_platewells = test1_platewells.merge(test2_platewells, how = "outer")
    n_expected_failures = len(combined_platewells)

    # Run the comparison
    test_out = test_aa(testdf_1, testdf_2)

    # Make sure we have as many fails as expected
    assert not test_out[0]
    assert len(test_out[1]) == n_expected_failures
    assert len(test_out[2]) == n_expected_failures