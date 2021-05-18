#!/bin/bash
# Define a list of folders to parse
declare -a refseqs=(
"./TestDatasets/DifferentOverlapDegrees/Set1/RefSeqs.csv"
"./TestDatasets/DifferentOverlapDegrees/Set2/RefSeqs.csv"
"./TestDatasets/DifferentOverlapDegrees/Set3/RefSeqs.csv"
"./TestDatasets/DifferentOverlapDegrees/Set4/RefSeqs.csv"
"./TestDatasets/DifferentPositions/Pos1/RefSeqs.csv"
"./TestDatasets/DifferentPositions/Pos2/RefSeqs.csv"
"./TestDatasets/DifferentPositions/Pos3/RefSeqs.csv"
"./TestDatasets/ForwardAndReverseDisagree/Set1/RefSeqs.csv"
"./TestDatasets/ForwardAndReverseDisagree/Set2/RefSeqs.csv"
"./TestDatasets/MixedWells/RefSeqs.csv"
"./TestDatasets/SingleSite/RefSeqs.csv"
"./TestDatasets/RealData/RefSeqs.csv"
)

declare -a directory=(
"./TestDatasets/DifferentOverlapDegrees/Set1"
"./TestDatasets/DifferentOverlapDegrees/Set2"
"./TestDatasets/DifferentOverlapDegrees/Set3"
"./TestDatasets/DifferentOverlapDegrees/Set4"
"./TestDatasets/DifferentPositions/Pos1"
"./TestDatasets/DifferentPositions/Pos2"
"./TestDatasets/DifferentPositions/Pos3"
"./TestDatasets/ForwardAndReverseDisagree/Set1"
"./TestDatasets/ForwardAndReverseDisagree/Set2"
"./TestDatasets/MixedWells"
"./TestDatasets/SingleSite"
"./TestDatasets/RealData/"
)


# Loop through the array
for i in ${!refseqs[@]};
do
	ssSeq ${refseqs[$i]} ${directory[$i]} --output ${directory[$i]}
	ssSeq ${refseqs[$i]} ${directory[$i]} --output ${directory[$i]} --troubleshoot
done

# Handle the separate flag for using a detailed reference sequence file
ssSeq ./TestDatasets/DifferentRefSeqByWell/RefSeqs.csv ./TestDatasets/DifferentRefSeqByWell --detailed_refseq --output ./TestDatasets/DifferentRefSeqByWell
ssSeq ./TestDatasets/DifferentRefSeqByWell/RefSeqs.csv ./TestDatasets/DifferentRefSeqByWell --detailed_refseq --output ./TestDatasets/DifferentRefSeqByWell --troubleshoot

# Report that we're moving on to test
echo "STARTING ssSEQ TESTS"

# Run the test script
python ./test_ssSeq.py