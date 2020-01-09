# ssSeq_Parser: Beta Release V2
Package for analyzing site-saturation next-generation sequencing data. This is still in beta, as thorough validation of functionality has not yet been performed. The code has been reorganized into easier to manage packages, so it should be easier for other developers to contribute now. Docstrings are still on the to-do list.

## Required Packages:
All required packages can be found in the conda environment file, "ssSeq.yml". Run ssSeq from within this conda environment to ensure all functionality can be used.

## Setting up a run:
### Command line interface
If added to PATH, ssSeq can be run directly from command line by calling "ssSeq" followed by the appropriate below arguments. If not added to path, ssSeq can be run from its home directory by calling "./ssSeq", again followed by the appropriate below arguments.

Positional arguments:
 - refseq: csv containing reference sequences. See more information below.
 - folder/fastq_f: A folder containing a set of fastq or fastq.gz to be used in ssSeq. If pointing to a folder, ssSeq will automatically match forward and reverse files by locating files with identical names (except for "_R1_" for a forward read file and "_R2_" for a reverse read file). If pointing to a file, it should be the fastq or fastq.gz file containing the forward reads.
  
 Optional arguments:
 - --fastq_r: Optional argument which points to the fastq or fastq.gz file containing the reverse reads. Must be specified if "folder/fastq_f" points to a file. Will be ignored if "folder/fastq_f" is a folder.
 - --analysis_only: Throw this flag to only perform analysis on the fastq/fastq.gz files. Alignment steps will be skipped.
 - --detailed_refseq: Throw this flag if you are passing in a detailed reference sequence file. See below for more detail.
 - --jobs: Controls the number of processors used for computation. Default is the number available on your machine minus one.
 - --troubleshoot: Throw this flag to force more detailed output, including consensus sequences and the underlying matrices used to calculate alignment frequencies and make variant calls.
 - --read_length: Optional argument for specifying the read length of your input files. If not specified, the read length will automatically be assumed to be the most common fragment length in the fastq/fastq.gz files.
 - --q_cutoff: The quality score filter, default = 30. Individual bases with a Q-score below this threshold will not be used for determining alignment frequencies.
 - --alignment_filter: A fraction of the maximum possible alignment score below which an alignment is discarded. This flag is used for eliminating off-target sequencing reads (e.g. from primer dimer). The default is 0.5, and it must vary between 0 and 1.
    
### RefSeqs.csv
Your reference sequences should be passed in to ssSeq as a csv file with the headers "PlateName", "IndexPlate", and "ReferenceSequence". An example csv file can be found in the ssSeq home directory. "PlateName" is a user-specified nickname for the plate in question. "IndexPlate" is the dual index plate used to prepare the ssSeq sample for sequencing -- it should take the form "DI##". "ReferenceSequence" is the gene fragment amplicon sequenced as part of ssSeq with the variable codons denoted by "NNN".