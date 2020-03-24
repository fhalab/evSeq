ssSeq
=====
Analyze raw fastq and fastq.gz files returned from next-gen sequencing of libraries prepared by the ssSeq wetlab protocol.

Table of Contents
-----------------
- [Installation](#Installation)
    - [Non-Programmers](#Non-Programmers)
    - [General Instructions](#General Instructions)
- [Program Arguments](#Program Arguments)
    - [Required Inputs](#Required Inputs)
    - [Optional Arguments](#Optional Arguments)
- [Working with the GUI](#Working with the GUI)
- [Working through Command Line](#Working through Command Line)
- [Biological Protocols](#Biological Protocols)
    - [Primer Design](#Primer Design)
    - [Library Preparation](#Library Preparation)

## Installation
### Non-Programmers
This section details installation of high level dependencies: gitbash (Windows users), git, and anaconda. If you have installed and are familiar with these items, you can skip this section. Installation on Linux is not detailed here, as we just assume you know what you're doing. Once this section is completed, continue installation by moving to [General Instructions](#General Instructions).

#### Installing git
Windows users: Install Git/Git Bash by following the instructions [here](https://www.stanleyulili.com/git/how-to-install-git-bash-on-windows/). Unless you know what you're doing, sticking to the default options during install is best. 

Mac users: Install Git by following the instructions [here](https://www.atlassian.com/git/tutorials/install-git). Some commands must be performed through the terminal. You can open a terminal by typing command + spacebar, then searching "Terminal"; opening the search result will give you a terminal window.

#### Installing Anaconda
Anaconda is an open-source package managment framework for scientific computing with Python. For details, look at their website [here](https://www.anaconda.com/). All software that supports ssSeq is handled by the Anaconda package manager. See below for installation instructions on both Windows and Mac:

Windows users: Install Anaconda following the instructions [here](https://docs.anaconda.com/anaconda/install/windows/). At step 8, I recommend adding Anaconda to your PATH environment variable. Note that this is in contrast to the recommendation of Anaconda, but their concerns shouldn't apply for our use case.

Mac users: Install Anaconda following the instructions [here](https://docs.anaconda.com/anaconda/install/mac-os/). 

#### Opening a Terminal Window
Later steps in our installation will require you to write commands in terminal. To open a terminal in Windows, right-click on either the desktop background or within your file explorer, then open "Git Bash". For Mac users, type command + spacebar, search "Terminal", then open the search result. 

#### Construction of "GitRepos" folder
The next step will be to install ssSeq. If you have not worked with Git repos before, I recommend creating a folder where you can store all of them. Wherever seems reasonable to you (most likely your home directory), create a folder called "GitRepos".

### General Instructions
#### ssSeq Installation
Open a terminal window and navigate to your GitRepos folder. This is accomplished by entering the below command in the terminal

	cd PATH_TO_GIT_REPOS

For instance, in my case I would type

	cd /home/brucejwittmann/GitRepos

After hitting "Enter" on your keyboard, you will notice that the prefix of your command line has changed. It will look something like the below after successfully executing the "cd" command:

![Command line example](./GitImages/CommandLineExample.png "Command line example")

From the ssSeq GitHub page, find the green box labeled "Clone or download". A screenshot giving the box location is below. Click on this box, then copy the presented url. 

![Clone or download example](./GitImages/CloneDownloadImage.png "Clone or download example")

Now in the command line, type 

	git clone COPIED_URL

replacing "COPIED_URL" with the link you just copied from GitHub. You can paste the link into command line by right-clicking and selecting "paste". 

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

ssSeq expects to see "NNN" in both the forward and reverse reads. Even if you don't have a variable position in the reverse read, you must denote one in a region of the reference sequence covered by sequencing; ssSeq will crash if you do not. 
