"""
Yeah, this one's meta. This file contains quick unit tests to make sure that the
functions used to confirm that true evSeq output matches expected evSeq output
in the stress testing.
"""
# Import required evSeq machinery
from .data_generation.globals import SAVELOC
from .data_generation.stress_tests import run_aa_stress_test

# Import 3rd party code
import os
import pandas as pd

# Get the paths to the aa files
SAMPLE1_PATH = os.path.join(SAVELOC, "sample_count_files", "sample1")
SAMPLE2_PATH = os.path.join(SAVELOC, "sample_count_files", "sample2")

# Load in the expected aa files
DECOUPLED_1 = pd.read_csv(os.path.join(SAMPLE1_PATH, "AminoAcids_Decoupled_All.csv"))
DECOUPLED_2 = pd.read_csv(os.path.join(SAMPLE2_PATH, "AminoAcids_Decoupled_All.csv"))

COUPLED_1 = pd.read_csv(os.path.join(SAMPLE1_PATH, "AminoAcids_Coupled_All.csv"))
COUPLED_2 = pd.read_csv(os.path.join(SAMPLE2_PATH, "AminoAcids_Coupled_All.csv"))

DECOUPLED_1_MAX = pd.read_csv(os.path.join(SAMPLE1_PATH, "AminoAcids_Decoupled_Max.csv"))
DECOUPLED_2_MAX = pd.read_csv(os.path.join(SAMPLE2_PATH, "AminoAcids_Decoupled_Max.csv"))

COUPLED_1_MAX = pd.read_csv(os.path.join(SAMPLE1_PATH, "AminoAcids_Coupled_Max.csv"))
COUPLED_2_MAX = pd.read_csv(os.path.join(SAMPLE2_PATH, "AminoAcids_Coupled_Max.csv"))

# Get the number of unique plate-well combos. This is how many failures we expect
# in each test.
TEST1_PLATEWELLS = DECOUPLED_1.loc[:, ["IndexPlate", "Well"]].drop_duplicates()
TEST2_PLATEWELLS = DECOUPLED_2.loc[:, ["IndexPlate", "Well"]].drop_duplicates()
N_EXPECTED_FAILURES = len(TEST1_PLATEWELLS.merge(TEST2_PLATEWELLS, how = "outer"))

# Write a helper function for comparing outputs
def helper(*args):
    
    # Run the function
    test_out = run_aa_stress_test(*args)

    # Make sure we have as many fails as expected
    assert not test_out[0] # Make sure we failed
    assert len(test_out[1]) == N_EXPECTED_FAILURES
    assert len(test_out[2]) == N_EXPECTED_FAILURES
    assert len(test_out[3]) == N_EXPECTED_FAILURES
    assert len(test_out[4]) == N_EXPECTED_FAILURES
    
# Run a comparison between the sets of files. They should not match.
def test_decoupled_1():
    helper(DECOUPLED_1, DECOUPLED_2, DECOUPLED_2_MAX)

def test_decoupled_2():
    helper(DECOUPLED_2, DECOUPLED_1, DECOUPLED_1_MAX)

def test_coupled_1():
    helper(COUPLED_1, COUPLED_2, COUPLED_2_MAX)

def test_coupled_2():
    helper(COUPLED_2, COUPLED_1, COUPLED_1_MAX)