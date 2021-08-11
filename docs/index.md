## Every Variant Sequencing: No sequence-function pair left behind.

Every Variant Sequencing (`evSeq`) provides resources for extremely low cost massively parallel sequencing of protein/enzyme variants arrayed in 96-well plates. This library preparation technique and data processing tool enables sequencing all variants produced during a protein engineering or biochemical mutagenesis experiment at a cost of cents per variant, even for labs that do not have expertise in or access to next-generation sequencing (NGS) technology.

### [Read the Paper!](one_day...)

### Read the Docs!
Navigate to [individual pages below](#documentation) or [start at the overview](bio/theory.md).

### General Overview
#### Workflow
The `evSeq` workflow
![Workflow](assets/figure2.png)

#### Sequence-Function Data
![SeqFunc](assets/figure3.png)

## Documentation
### Biology
#### [Theoretical overview](bio/theory.md)

#### [Library preparation](bio/lib_prep.md)
- [Dual-Index Barcode Plates](bio/lib_prep.md#dual-index-barcode-plates)
- [Inner Primer Design](bio/lib_prep.md#inner-primer-design)
- [Inner Primer Test PCR](bio/lib_prep.md#inner-primer-test-pcr)
- [PCR Protocol](bio/lib_prep.md#pcr-protocol)
- [PCR Product Purification](bio/lib_prep.md#pcr-product-purification)

### Computation
#### [Computational basics](comp/basics.md)
#### [Installation](comp/installation.md)
- [Installing from GitHub with the `conda` environment](comp/installation.md#installing-from-github-with-the-conda-environment)
- [Standard `pip` Install and Dependencies](comp/installation.md#standard-pip-install-and-dependencies)
#### [Running `evSeq`](comp/usage.md)
- [Post Innstallation](comp/usage.md#post-installation)
- [Using `evSeq` from the command line or GUI](comp/usage.md#using-evseq-from-the-command-line-or-gui)
- [Required Arguments](comp/usage.md#required-arguments)
  - [The `refseq` file](comp/usage.md#the-refseq-file)
  - [`folder`](comp/usage.md#folder)
- [Optional Arguments](comp/usage.md#optional-arguments)
#### [Understanding the outputs](comp/outputs.html)
- [Qualities](comp/outputs.html#qualities)
- [OutputCounts](comp/outputs.html#outputcounts)
- [Platemaps](comp/outputs.html#platemaps)

<!--
#### [Using `evSeq` tools in a Python environment](comp/additional.html)
- [Running `evSeq` functions]()
- [Mapping sequence to function]()
-->
### Troubleshooting
- [Poor reverse read quality](troubleshooting.md#poor-reverse-read-quality)
- [Poor results but good quality sequencing](troubleshooting.md#poor-results-but-good-quality-sequencing)
- [Windows: `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'`](troubleshooting.md#windows-commandnotfounderror-your-shell-has-not-been-properly-configured-to-use-conda-activate)
- [macOS: `PermissionError: [Errno 1] Operation not permitted`](troubleshooting.md#macos-permissionerror-errno-1-operation-not-permitted)
