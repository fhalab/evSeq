## Every Variant Sequencing: No sequence-function pair left behind.

Every Variant Sequencing (`evSeq`) is a library preparation and analysis protocol designed to slot neatly into existing protein engineering workflows to enable extremely low cost massively parallel sequencing of heterologously expressed protein variants arrayed in 96-well plates. This workflow enables sequencing all variants produced during a protein engineering or biochemical mutagenesis experiment at a cost of cents per variant, even for labs that do not have expertise in or access to next-generation sequencing (NGS) technology.

### Read the Paper!
This repository accompanies the work ["evSeq: Cost-Effective Amplicon Sequencing of Every Variant in Protein Mutant Libraries"](LINK_TO_PAPER). If you use this tool, please [cite us](LINK_TO_PAGE_WITH_CITATION_FORMATS).

### Read the Docs!
Navigate to [individual pages below](#documentation) or [start at the overview](bio/theory.md).

### General Overview
#### The `evSeq` workflow
![Workflow](assets/figure2.png)
**A)** Laboratory procedure. **B)** Computational procedure.

#### Construct and visualize sequence-function pairs
![SeqFunc](assets/figure3.png)
Sequencing eight site-saturation libraries (768 wells) in a single `evSeq` run and combining this with activity data to create low-cost sequence-function data. **A)** Enzyme and active-site structure highlighting mutated residues. **B)** Heatmap of the number of identified variants/mutations ("counts") for each position mutated ("library") from processed `evSeq` data. **C)** Heatmap of the average activity ("normalized rate") for each variant/mutation in each library. **D)** Counts for a single library, also showing the number of unidentified wells. **E)** Activity for a single library, showing biological replicates. (Inset displays the mutated residue in this library.)

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
- [Post Installation](comp/usage.md#post-installation)
- [Using `evSeq` from the command line or GUI](comp/usage.md#using-evseq-from-the-command-line-or-gui)
- [Required Arguments](comp/usage.md#required-arguments)
  - [The `refseq` file](comp/usage.md#the-refseq-file)
  - [`folder`](comp/usage.md#folder)
- [Optional Arguments](comp/usage.md#optional-arguments)
#### [Understanding the Outputs](comp/outputs.html)
- [`Qualities`](comp/outputs.html#qualities)
- [`OutputCounts`](comp/outputs.html#outputcounts)
- [`Platemaps`](comp/outputs.html#platemaps)
- [`evSeqLog`](comp/outputs.html#evSeqLog)
- [`ParsedFilteredFastqs`](comp/outputs.html#parsedfilteredfastqs)
- [`Alignments`](comp/outputs.html#alignments)

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
