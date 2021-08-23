"""
Module for comparing evSeqOutputs
"""
import os
import glob
import pandas as pd


def compare_files(file1, file2, obj='files'):
    """Compares two csv files. Returns an exception if one is given."""
    file1_df = pd.read_csv(file1)
    file2_df = pd.read_csv(file2)

    try:
        pd.testing.assert_frame_equal(file1_df, file2_df, obj=obj)
        e = None
    except AssertionError as e:
        pass

    return e


def compare_to_expected(test_path, expected_path, print_errors=True):
    """Compares the OutputCounts files of two evSeq runs, given a path
    to a date-time output result of evSeq for each argument.
    """
    # Go through sorted OutputCounts
    expected_counts = sorted(
        glob.glob(os.path.join(expected_path, 'OutputCounts/*'))
    )
    test_counts = sorted(
        glob.glob(os.path.join(test_path, 'OutputCounts/*'))
    )

    # Validate each file
    bad_files = {}
    for exp, test in zip(expected_counts, test_counts):
        # Get expection (or None)
        e = compare_files(exp, test, obj='OutputCounts files')

        # Add file and exception to list of bad files
        if e is not None:
            bad_files[test] = e

    # Print errors
    if bad_files and print_errors:
        print(
            f'Found {len(bad_files)} bad files of {len(expected_counts)} total.'
        )
        print('\n')
        for file, message in bad_files.items():
            print(file)
            print(''.join(['-']*len(file)))
            print(message)
            print('\n')

    if bad_files:
        raise AssertionError(
            'Outputs are not the same. See printed errors for details.'
        )