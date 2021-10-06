"""
Runs stress tests of evSeq. This script will randomly generate input data and
input conditions whose expected outputs from evSeq are known. The generated data
is then run through evSeq and the evSeq outputs are compared to the known outputs.
Any situation where evSeq outputs do not match known outputs are printed to terminal
and saved. From the root evSeq directory, outputs for these tests will be saved
at ./tests/test_data/evSeqOutput and error reports will be saved at 
./tests/test_data/ErrorReports

This script will run indefinitely until you kill it.
"""
# Import required modules
import argparse
from tests.data_generation.stress_tests import run_evseq_stress_test

# Define the core function
def main():
    
    # Instantiate argparser
    parser = argparse.ArgumentParser()
    
    # Add arguments
    parser.add_argument("--detailed",
                        help = "Raise this flag to run tests using detailed evSeq inputs.",
                        required = False,
                        action = "store_true")
    parser.add_argument("--include_nnn",
                        help = "Raise this flag to run tests where 'NNN' is used to designate variable positions.",
                        required = False,
                        action = "store_true")
    parser.add_argument("--seed",
                        help = "Integer used to seed the run. Default = 0.",
                        default = 0,
                        required = False,
                        type = int)
    parser.add_argument("--keep_output",
                        help = "Raise this flag to keep ALL output. By default, this script deletes evSeq output for which all tests passed.",
                        required = False,
                        action = "store_true")
    
    # Parse argument
    args = parser.parse_args()
    
    # Run the stress test
    run_evseq_stress_test(args.detailed,
                          args.include_nnn,
                          keep_output = args.keep_output,
                          seed = args.seed)

if __name__ == "__main__":
    main()