# ssSeq_Parser: Beta Release
Package for analyzing site-saturation next-generation sequencing data. This is a beta-release, so there will be bugs. Let me know what features you would like changed/added as well as what bugs you come across. 

The code behind this is quite complicated, so I know I will have made logical errors somewhere that don't take full advantage of the dataset. If you feel up to it, I'd really appreciate feedback on the code/logic checking. I annotate quite extensively, but just come and talk to me if you want to discuss any of the decisions I made. Also, I apologize that the code isn't the neatest...once we've settled on features and the overall infrastructure is set I'll organize everything so that other developers can more easily contribute.

## Note on resource usage
Any kind of next-gen sequence analysis is computationally intensive. Analyzing a single run file will use ~6 GB of RAM in troubleshoot mode, ~4 GB in regular mode. Drop the number of processors and/or do not run in ts_mode if you can't devote this amount of RAM. If anyone has suggestions for lowering resource intensiveness, I'm all ears.

## Required Packages:
- Anaconda default distribution
- Biopython (conda install -c anaconda biopython)

## Setting up a run:
### RefSeqs.csv
This contains the reference sequences to be used in each plate and well. Fill out the columns as below:
1. "PlateName": Whatever nickname you want to give to your plate.
2. "IndexPlate": Whatever index plate you used to barcode the samples of this plate. Index plate 1 is DI01, index plate 2 is DI02, and so on.
3. "Well": The well of the sample in question.
4. "NumberOfVariablePositions": How many amino acid positions do you expect to be variable in this library? For single-site saturation mutagenesis, this is 1, for triple-site, this is 3, etc.
5. "F-RefSeq": The reference sequence to be used for your forward read. The codons that are variable should be written as "NNN". Use the sense strand as the reference sequence.
6. "R-RefSeq": The reference sequence to be used for your reverse read. The codons that are variable should be written as "NNN". Use the sense strand as the reference sequence (do NOT use the reverse complement to the reference).

### GlobalParams.txt
This file contains the parameters used by the script. These include the below:
1. "RefSeq_Filepath": The path to the RefSeqs.csv file.
2. "ForwardReads_Filepath": The path to the forward NGS reads. This is the file that contains "R1" in the name. Ignore "I1" files.
3. "ReverseReads_Filepath": The path to the reverse NGS reads. This is the file that contains "R2" in the name. Ignore "I2" files.
4. "ts_mode": Whether or not to run in troubleshooting mode. Running outside of ts_mode will return summary information about the calls made at the variable positions as well as consensus sequences. Troubleshooting mode is more resource intensive and will output additional information about the run, including frequency count matrices used for calculations, all alignments made, and pickles of the Python objects generated during the run. 
5. "n_jobs": The number of processors to use for processing. Multiprocessing is performed on the well level to maximize memory efficiency. Processing multiple wells simultaneously will still require more RAM, however.
6. "q_cutoff": The quality score cutoff to use when processing NGS reads. Each basepair read by a MiSeq is assigned a quality score which directly corresponds to the confidence in the basepair call. Low quality reads will result in miscalls by this and any other NGS data analysis program. Lowering the quality score cutoff will result in more reads per reference basepair used, but lower overall confidence in the calls made. Raising the quality score cutoff will result in less reads per reference basepair used, but better overall confidence in the calls made.
7. "alignment_cutoff": The per-base average alignment score against a reference needed for a read to be used. A score of 0.5 approximately means that half of the bases must match the reference for a read to be used in analysis. Lowering the alignment cutoff means that reads not pertaining to your reference sequence (such as leftover primer dimer) will be more likely to be used in analysis. Raising the alignment cutoff means that only reads that match the reference sequence very well will be used. This parameter is meant to filter out noise in the reads, and so should be kept low! If it is too high, then you will might miss certain types of critical errors in library preparation.
8. "OutputLocation": The name and location of the folder to which all output data will be saved.
9. "id_regex": The regular expression used for parsing the id line of each NGS read call. This parameter should not be changed. 

### ssSeq_Parser.py
Once the reference sequences and parameters are set up, move to the directory containing "ssSeq_Parser.py" in terminal. Run the program from terminal using "python ssSeq_Parser.py".



