ssSeq
=====
Analyze raw fastq and fastq.gz files returned from next-gen sequencing of libraries prepared by the ssSeq wetlab protocol.

Table of Contents
-----------------
- [Installation](#Installation)
    - [Non-Programmers](#Non-Programmers)
        - [Installing Git](#Installing-Git)
        - [Installing Anaconda](#Installing-Anaconda)
        - [Opening a Terminal Window](#Opening-a-Terminal-Window)
        - [Construction of "GitRepos" Folder](#Construction-of-GitRepos-Folder)
    - [General Instructions](#General-Instructions)
        - [ssSeq Installation](#ssSeq-Installation)
        - [Conda Environment Setup](#Conda-Environment-Setup)
            - [Dependencies](#Dependencies)
        - [PATH Variable Setup](#Path-Variable-Setup)
- [Program Arguments](#Program-Arguments)
    - [Required Inputs](#Required-Inputs)
        - [refseq](#refseq)
        - [folder](#folder)
    - [Optional Arguments](#Optional-Arguments)
- [Working with the GUI](#Working-with-the-GUI)
- [Working through Command Line](#Working-through-Command-Line)
- [Example Data](#Example-Data)
- [Understanding ssSeq Output](#Understanding-ssSeq-Output)
    - [Summaries](#Summaries)
        - [MaxInfo.csv](#MaxInfo.csv)
        - [SummaryInfo.csv](#SummaryInfo.csv)
        - [VariantInfo.csv](#VariantInfo.csv)
    - [Platemaps](#Platemaps) 
    - [Qualities](#Qualities)
    - [Alignments](#Alignments)
    - [AACountsFrequencies](#AACountsFrequencies)
    - [BPCountsFrequencies](#BPCountsFrequencies)
    - [ConsensusSequences](#ConsensusSequences)
    - [ssSeqLog](#LogFile)
- [Biological Protocols](#Biological-Protocols)
    - [Primer Design](#Primer-Design)
    - [Library Preparation](#Library-Preparation)


# Installation
## Non-Programmers
This section details installation of high level dependencies: gitbash (Windows users), git, and anaconda. If you have installed and are familiar with these items, you can skip this section. Installation on Linux is not detailed here, as we just assume you know what you're doing. Once this section is completed, continue installation by moving to [General Instructions](#General Instructions).

### Installing Git
Windows users: Install Git/Git Bash by following the instructions [here](https://www.stanleyulili.com/git/how-to-install-git-bash-on-windows/). Unless you know what you're doing, sticking to the default options during install is best. 

Mac users: Install Git by following the instructions [here](https://www.atlassian.com/git/tutorials/install-git). Some commands must be performed through the terminal. You can open a terminal by typing command + spacebar, then searching "Terminal"; opening the search result will give you a terminal window.

### Installing Anaconda
Anaconda is an open-source package managment framework for scientific computing with Python. For details, look at their website [here](https://www.anaconda.com/). All software that supports ssSeq is handled by the Anaconda package manager. See below for installation instructions on both Windows and Mac:

Windows users: Install Anaconda following the instructions [here](https://docs.anaconda.com/anaconda/install/windows/). At step 8, I recommend adding Anaconda to your PATH environment variable. Note that this is in contrast to the recommendation of Anaconda, but their concerns shouldn't apply for our use case.

Mac users: Install Anaconda following the instructions [here](https://docs.anaconda.com/anaconda/install/mac-os/). 

### Opening a Terminal Window
Later steps in our installation will require you to write commands in terminal. To open a terminal in Windows, right-click on either the desktop background or within your file explorer, then open "Git Bash". For Mac users, type command + spacebar, search "Terminal", then open the search result. 

### Construction of "GitRepos" Folder
The next step will be to install ssSeq. If you have not worked with Git repos before, I recommend creating a folder where you can store all of them. Wherever seems reasonable to you (most likely your home directory), create a folder called "GitRepos".

## General Instructions
### ssSeq Installation
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

replacing "COPIED_URL" with the link you just copied from GitHub. You can paste the link into command line by right-clicking and selecting "paste". If you successfully installed Git earlier, then this should begin the installation process. 

### Conda Environment Setup
We will now set up the conda environment for ssSeq. This conda environment will neatly package all of the support software needed to run ssSeq. To begin, navigate to the ssSeq folder via the command line. On my machine, I would type 

    cd /home/brucejwittmann/GitRepos/ssSeq

Note that ssSeq was installed within my GitRepos folder. Depending on where you installed ssSeq you will need to navigate to a different folder. Next, type the command 

    conda env create -f ssSeq.yml

If you successfully installed conda earlier, then this command should run without problems. 

#### Dependencies
Advanced users: If you would rather not use the ssSeq environment and run in a custom environment (or, if you're a brave soul, your base environment), below are the ssSeq dependencies, all of which are available through conda. These dependencies are, of course, explicitly listed in the ssSeq.yml environment file. The explicit version call with bokeh handles an incompatibility between the most recent bokeh and holoviews versions at the time of writing (03/24/2020).

    # ssSeq without GUI
    - biopython
    - colorcet
    - holoviews
    - bokeh=1.4.0
    - numpy
    - pandas
    - python=3.7
    - tqdm
    - scipy
    
    # With GUI, also need:
    - gooey

### PATH Variable Setup
This step is optional. However, if you want to run ssSeq from anywhere on your computer then you should do it. For the uninitiated, your PATH variable contains the locations your computer will look for a program when you ambiguously tell it to run something. By putting the location of ssSeq in my PATH variable, this means I can run ssSeq from anywhere on my computer by just calling

    ssSeq ARGS

instead of having to call

    /home/brucejwittmann/GitRepos/ssSeq/ssSeq ARGS

Instructions for adding ssSeq to PATH are [here](https://helpdeskgeek.com/windows-10/add-windows-path-environment-variable/) for Windows, [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/#.Uydjga1dXDg) for Mac, and [here](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path) for Linux (Ubuntu). You should add the directory containing the ssSeq GitHub repository to your PATH variable. For me, this means I would add the below to PATH

    /home/brucejwittmann/GitRepos/ssSeq/

You may also need to make the ssSeq run files executable. In Mac and Linux, accomplish this by first navigating to the ssSeq folder in command line. For instance, for me this is

    cd /home/brucejwittmann/GitRepos/ssSeq/

Once in the ssSeq folder, make both the GUI and command line versions executable by writing the below commands in command line

    chmod +x ssSeq
    chmod +x ssSeqGui

ssSeq should now be fully ready for your use.

# Program Arguments
## Required Arguments
The only two required arguments for ssSeq are a table giving the reference sequences for the expected amplicons ("refseq") and the folder containing the fastq files resulting from the sequencing run ("folder"). Both of these arguments are explained in detail below.

### refseq

### folder

## Optional Arguments


# Working with the GUI
The GUI is designed for use by non-programming experts. If you are comfortable with a command line interface, that is the recommended way to use ssSeq. If using the GUI, make sure you check the log file after each run to check for warnings or errors encountered. See details on the log file [here](#ssSeqLog). 

When using the GUI, begin by opening a terminal window then activating the ssSeq conda environment as below

    conda activate ssSeq

With the conda environment active, you can now launch the graphic user interface (GUI). There are two ways to launch ssSeq through the GUI:
1. If ssSeq was added to your PATH and ssSeqGui made executable (see [PATH Variable Setup](#PATH-Variable-Setup)), then the GUI can be opened by typing

    ssSeqGui

2. If ssSeq was not added to your PATH and is not executable, then you can activate the GUI by first navigating to the ssSeq git repo folder (installed above) through command line and explicitly invoking Python as below

    cd ssSEQ_LOCATION
    python ssSeqGui

Executing either of the above commands will launch an instance of the GUI. It should look like the below:

![GUI](./GitImages/GUI.png "GUI")

Note that the two required arguments are at the top of the GUI, details on these arguments can be found [here](#Required-Arguments). All other arguments can be accessed by scrolling down. To populate "refseq", browse to the location with your completed "refseq" csv file and click "ok". To populate "folder", browse to the folder with your fastq or fastq.gz files and click "Open". The typical ssSeq run can be started once these fields are populated. Click "Start" in the main GUI window to run ssSeq; the progress of the program will be printed to the GUI along with any encountered warnings and errors. 

Other arguments can also be specified in the GUI. See [OptionalArguments](#Optional-Arguments) for details.

# Working through Command Line
This is the recommended way to use ssSeq, as any warnings or errors encountered are printed directly to the terminal. Begin by opening a terminal window then activating the ssSeq conda environment as below

    conda activate ssSeq 

With the conda environment active, ssSeq can be run. There are two ways to launch ssSeq:
1. If ssSeq was added to your PATH and ssSeq made executable (see [PATH Variable Setup](#PATH-Variable-Setup)), then ssSeq can be run by typing

    ssSeq refseq folder OPTIONAL_ARGS FLAGS

2. If ssSeq was not added to your PATH and  is not executable, then you can run ssSeq by first navigating to the ssSeq git repo folder (installed above) through command line and explicitly invoking Python as below

    cd ssSEQ_LOCATION
    python ssSeq refseq folder OPTIONAL_ARGS FLAGS

In both of cases, "refseq" is the path to your reference sequence file and "folder" is the location of the folder with your fastq or fastq.gz files. Details on these required arguments can be found [here](#Required-Arguments). Optional arguments and flags are passed in after the two positional arguments. For information on the potential optional arguments and flags, type

    ssSeq -h

or

    python ssSeq -h

depending on whether or not you added ssSeq to PATH. The "-h" flag will pull up the help window detailing all possible ssSeq arguments. Note that these arguments are also detailed in [Optional Arguments](#Optional-Arguments). The help window will look like below:

![Help Window](./GitImages/HelpWindow.png "Help Window")

# Example Data
The folder InstallationConfirmationData contains an example reference sequence file and example fastq.gz files. These files can be used for confirming installation or just playing around with ssSeq. 

# Understanding ssSeq Output
The output location of ssSeq is controlled with the "output" optional argument (see [here](#Optional-Arguments)). If the "output" argument is not set, then ssSeq will save to the current working directory (command line) or the ssSeq Git repository folder (GUI). If the save location has not previously been used, then ssSeq will create a folder titled "ssSeq_Output" in the output location which contains a folder giving the date-time of the run initialization (in yyyymmdd-hhmmss format). If the save location has been previously used, then ssSeq will add another date-time folder with the previously generated ssSeq_Output folder. All ssSeq outputs of a specific run (except the log file) are contained in the associated date-time folder. The below sections detail the folders found within the date-time folder.

## Summaries
The summaries 

### MaxInfo.csv

### SummaryInfo.csv

### VariantInfo.csv

## Platemaps

## Qualities

## Alignments

## AACountsFrequencies

## BPCountsFrequencies

## ConsensusSequences

## ssSeqLog
ssSeq keeps a log of every run. This is the only output not found in the generate date-time folder. Within the ssSeq Git repository folder, the log can be found here: ssSeqSupport/ssSeqLog.log. Information captured by the log file includes:

1. The start time of the ssSeq run, given as 'yyyymmdd-hhmmss' followed by a series of underscores. This is the first line of each log block.
2. The values of all parameters input to ssSeq. Note that if parameters are unspecified, the log records the default parameters.
3. Calculated parameters, including
    1. The forward and reverse read file pairs identified in the 'folder' argument
    2. Any files within 'folder' that were not matched. 
    3. The calculated read length for the run if read length was unspecified. Otherwise, the input read length.
4. Any warnings encountered during the run. These warnings will also be printed to the console during the run. 
5. Fatal errors. If the program completed successfully, the last line in the log entry will read "Run completed. Log may contain warnings."

The amount of information stored in the log file is small (bytes per run), but will build with continued use of ssSeq. If the file gets too large (this will take a long time...) you can delete ssSeqLog.log; on the next run a fresh ssSeqLog.log file will be instantiated.







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
