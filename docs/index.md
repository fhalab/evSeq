## **E**very-**V**ariant **Seq**uencing: No sequence-function pair left behind.

Every-Variant Sequencing (`evSeq`) provides resources for extremely low cost massively parallel sequencing of protein/enzyme variants arrayed in 96-well plates. This library preparation technique and data processing tool enables sequencing all variants produced during a protein engineering or biochemical mutagenesis experiment at a cost of cents per variant, even for labs that do not have expertise in or access to next-generation sequencing (NGS) technology.

[Read the paper!](one_day...)

### Overview
![Overview](assets/overview.png)

## Documentation
### Biology
#### [Theoretical overview](bio/theory.md)

#### [Library preparation](bio/lib_prep.md)
- [Inner Primer Design](bio/lib_prep.md#inner-primer-design)
- 

### Computation
#### [Programming basics](comp/basics.md)
#### [Installation](comp/installation.md)
- [Installing from GitHub](comp/installation.md#installing-from-github)
- [Installing with PyPI (pip)](comp/installation.md#installing-from-pypi)
- [Using the `evSeq` environment](comp/installation.md#using-the-evseq-environment)
#### [Running `evSeq`](comp/usage.md)
- [Post Installation](comp/usage.md#post-installation)
- [The `refseq` file](comp/usage.md#the-refseq-file)
- [Using `evSeq` from the command line or GUI](comp/usage.md#using-evseq-from-the-command-line-or-gui)
- [Using `evSeq` in a Python environment](comp/usage.md#using-evseq-in-a-python-environment)
#### [Understanding the outputs](comp/outputs.md)
- [Qualities](comp/outputs.md#qualities)
- [Summaries](comp/outputs.md#summaries)
- [Platemaps](comp/outputs.md#platemaps)
#### Additional tools
- Mapping sequence to function

### Troubleshooting
- Poor reverse read quality
- Poor alignments but good quality

## Usage
```
evSeq refseq data_dir OPTIONAL_ARGS FLAGS
```
where `evSeq` has been set up in PATH, `refseq` is a csv file containing information about the experiment, and `data_dir` contains the raw fastq files for the experiment.

For information on optional arguments and flags, run `evSeq -h`.

## Output
`evSeq` returns sequencing information mapped to individual wells in up to eight 96-well plates per single sample in a multiplexed NGS experiment:

[plate map?]

## Installation
For non-programmers and those unfamiliar with [Anaconda](https://www.anaconda.com/) or [GitHub](https://www.github.com), see the [programing basics page](basics.md) for information on how to set up your computer environment to run evSeq.
