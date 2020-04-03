Validation Information
======================
This folder contains multiple different synthetic datasets designed to test the functionality of ssSeq. It is designed for use by developers who want to make sure they haven't broken anything in the code after making changes to the software.

# Quick Use
Running inside of the ssSeqTest conda environment (see ssSeqTest.yml), the shell script ssSeqTest.sh should be run to test the existing instance of ssSeq. This script will...

1. Run ssSeq in both troubleshooting and standard mode for each test dataset found in TestDatasets/.
2. Run test_ssSeq.py to evaluate the output of ssSeq and make sure that it matches the expected output (which has been precalculated and uploaded).

Importantly, the expected data is generated using current alignment rules and default parameters. Obviously, if you change either the rules of the defaults, then you're going to need to develop new expected datasets.

After running ssSeqTest.sh, image outputs from ssSeq should be evaluated manually to confirm there are no errors -- tests are not written for this output type. Note also that tests are currently not run for the troubleshooting output -- it is only for the summary output generated in non troubleshoot mode.

# General notes on test datasets and expected datasets:
1. Only barcode plates 1-8 are currently deployed and tested in ssSeq. Further test data must be made to incorporate more barcode plates.
2. There are some discrepencies between the expected and ssSeq output. These all happen as a result of the alignment algorithm: Insertions and deletions toward the end of the sequence are harder to capture, so ssSeq does not always catch insertions and deletions at the end of a sequence. Errors in alignment can be remedied with overlapping forward and reverse reads! When reads do not overlap, some insertions and deletions may go unnoticed. All exceptions were manually evaluated for this problem, and code is written in the test data to ignore expected discrepencies. 
    

# test_ssSeq.py
This python script contains all tests that will be run on ssSeq output. All developers should feel free to add more tests to/improve the detail of the tests in this script as issues are encountered. 

# TestDatasets Folder Architecture
The folder TestDatasets contains a folder for each style of test dataset. Each test dataset folder contains...

1. Forward and reverse FastQ files pertaining to the test.
2. A reference sequence file pertaining to the test.
3. A folder containing expected output. Currently, only the expected output for the summary data is included.

The different test styles are...

1. DifferentOverlapDegrees: This folder contains subfolders with 4 different 4-site combinatorial libraries. The first set has one overlapping position between the forward and reverse reads, the second has two overlapping positions, and so on.
2. DifferentPositions: This folder contains subfolders with 3 different 4-site combinatorial libraries. No positions overlap in any of these sets, but the positions of mutation change in each set.
3. DifferentRefSeqByWell: This folder contains a test dataset where different wells have different reference sequences. Different reference sequences have different mutated positions and different degrees of overlap. All reference sequences are derived from the "DifferentOverlapDegrees" reference sequences. Note that this data is designed explicitly as a combination test of "DifferentOverlapDegrees" and "DifferentPositions" as well as to test the "detailed_refseq" style of running ssSeq.
4. ForwardAndReverseDisagree: This folder contains subfolders with 2 different 4-site combinatorial libraries. In these tests the forward and reverse reads do not necessarily agree sequence-wise on overlapping positions. The first set has one overlapping position while the second set has two overlapping positions. 
5. MixedWells: This folder contains a single test dataset of a 4-site combinatorial library. In this dataset, multiple sequences have been added to each well. This dataset is to test how well ssSeq distinguishes between different sequences in the same well.
6. RealData: This folder contains the output of a real-life ssSeq run (all other datasets are synthetic). The output of this dataset should always be the same. It is included as a failsafe in case any assumptions not directly tested in the synthetic tests no longer hold after a change to the ssSeq software.
7. SingleSite: This folder contains a single test dataset where only a single position is mutated (all other test datasets are 4-site). This test is included to make sure ssSeq works for both combinatorial and single-site libraries. Note that both forward and reverse reads cover the single site.

For each sequence pair in each test dataset, the following "noise" conditions are added in at random...
1. The mutated site may have a Q-score below the target threshold.
2. The read may have been replaced with a noisy sequence. This is to make sure that alignment filtering is working correctly. 

Each well also has random DNA fragments added to it. These should be filtered out by alignment filtering.

# On Generation of Expected Data
The synthetic datasets were all constructed in a Jupyter Notebook. The notebook is included in html and .ipynb form in this folder "ssSeqValidationDatasetConstruction.ipynb/html". 

The expected data was generated using the python script "GenerateExpectedData.py", which is included in this folder. If ssSeq rules or defaults change, this script will need to be modified to generate the expected test data.