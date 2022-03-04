# Troubleshooting
For any issues that are not covered below, please reach out to us for assistance. Issues with the software should be reported on GitHub [here](https://github.com/fhalab/evSeq/issues); we also have a forum on GitHub for discussing questions/problems with running the `evSeq` protocol (both computational and wet-lab) that can be found [here](https://github.com/fhalab/evSeq/discussions/).

## Poor reverse read quality
If your forward read quality is excellent but your reverse read quality is very poor, this can ruin your `Coupled` output files. While you can still get useful information from the forward reads in the `Decoupled` output files, this issue usually comes from library preparation and should be fixed as soon as possible.

This can occur due to certain primer designs. While the `evSeq` adapters have been optimized (both computationally and through trial-and-error) to work with most inner primer sequences, we can't expect this to always be perfect. If possible, re-design your inner primers in slightly new locations (re-designing both and testing all combinations is usually a safe bet).

This also occurs due to too much primer-dimer in your submitted sample. This is a subset of the above problem (better primer design means less primer-dimer) but can be resolved by running a longer agarose gel to separate the primer-dimer (~75-150 bp) from the actual barcoded `evSeq` amplicon. This is especially important when the amplicon is small.
## Poor results but good quality sequencing
If your plate seems to contain good quality reads but mostly dead wells, and your average read length in the fastqs is not unreasonably low, it is likely that your `refseq` file information is not correct for the given sequence. Often this comes from having the wrong sequence information (`FPrimer/RPrimer`s and/or `VariableRegion`, setting the alignment process off.

#### Tips to solve
- Ensure that your reference sequence information accurately reflects the sequences you sent your NGS provider. Small mistakes in the reference frame can ruin the analysis.
- Run `evSeq` with the `--keep_parsed_fastqs` or `--only_parse_fastqs` flag to confirm that they match your reference sequence information.
- Run `evSeq` with the `--return_alignments` flag to see where the alignment of each read to your reference sequence is failing.

## Progress bar not showing up in Jupyter
This is likely because you do not have the `ipywidgets` extension enabled. Run the below command:

```
jupyter nbextension enable --py widgetsnbextension
```

## Windows: `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'`
On Windows, you may receive the below error the first time you try to activate an environment:
```
CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
```
The error can be fixed by entering the command
```
conda init bash
```
and then repeating the `conda activate` command. This step should fix the error permanently.

## Windows: The GUI Will Not Open
Double clicking on the desktop application should open the GUI. If it does not and you are running Windows, the issue is likely caused by one of two things:

1. Most commonly, this is caused by Anaconda not being in your PATH environment variable. The easiest fix is to start over by reinstalling Anaconda and making sure to check the box that adds it to your PATH environment variable. Note that the default install on Windows does *not* add Anaconda to your path, and so you will need to pay special attention to the options you select during Anaconda install. **Warning: Reinstalling Anaconda will remove any environments that you have saved already (so you may want to look into alternate solutions if you've had Anaconda on your system for a while).** Once you have reinstalled Anaconda, you can follow the [evSeq installation instructions](3-installation.md#installing-from-github-with-the-conda-environment) to reinstall `evSeq`. An alternate option to reinstalling is to manually add Anaconda to your PATH variable; you should talk to someone with programming experience to accomplish this, however. 

2. We have also seen this problem when users install `evSeq` through Windows Command Prompt rather than Git Bash. Instructions for opening Git Bash are [here](2-basics.md#opening-a-terminal-window). To solve this problem, you will first need to remove the `evSeq` conda environment. This can be accomplished by entering the below into a Git Bash terminal:

```
conda env remove -n evSeq
```

Then, follow the [evSeq installation instructions](3-installation.md#installing-from-github-with-the-conda-environment) to reinstall `evSeq` using Git Bash rather than Command Prompt.

## macOS: `PermissionError: [Errno 1] Operation not permitted`
Any time after upgrading to a newer macOS, you might randomly find that `evSeq` has stopped working and gives `PermissionError: [Errno 1] Operation not permitted`.

This is a security issue for software being run on macOS, and can be solved by following [these directions](https://stackoverflow.com/questions/58479686/permissionerror-errno-1-operation-not-permitted-after-macos-catalina-update) in your System Preferences/Security and Privacy settings.

## macOS: `xcrun` Error
If `xcode` is not installed on Mac, then you will see an error like the below:

`Error: xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun`

To fix this, install `xcode` by running the below in the terminal:

`xcode-select --install`

## Linux: `CondaEnvException: Pip failed`
If you are running on a Linux distribution and see the error `CondaEnvException: Pip failed` during `conda install`, it is likely an issue with the wxPython installation. To confirm that this is a wxPython problem, look up a few lines from the bottom of the traceback -- if you see the line `No package 'gtk+-3.0' found` above a few lines beginning with `***`, then this is a wxPython installation problem. Documentation on the challenges of installing wxPython on Linux can be found [here](https://wxpython.org/blog/2017-08-17-builds-for-linux-with-pip/index.html) and [here](https://github.com/wxWidgets/Phoenix/issues/1831). As a quick fix, you may be able to run `sudo apt-get install build-essential libgtk-3-dev` (modified as appropriate for installing packages on different distributions) before re-attempting installation.

---
*Back to the [main page](index.md).*