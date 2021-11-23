# evSeq
### *No sequence-function pair left behind.*

Every Variant Sequencing (`evSeq`) is a library preparation and analysis protocol that **slots neatly into existing workflows to enable extremely low-cost, massively parallel sequencing of protein variants**. Designed for heterologously expressed protein variants arrayed in 96-well plates (or similar), this workflow enables **sequencing all variants** from targetted mutagenesis libraries produced during a protein engineering or biochemical mutagenesis experiment at a cost of **cents per variant**, even for labs that do not have expertise in or access to next-generation sequencing (NGS) technology.

### Read the Paper!
This repository accompanies the work ["evSeq: Cost-Effective Amplicon Sequencing of Every Variant in Protein Mutant Libraries"](https://doi.org/10.1101/2021.11.18.469179).
<details>
<summary>If you use this tool, please cite us:</summary>

```
TY  - JOUR  
T1  - evSeq: Cost-Effective Amplicon Sequencing of Every Variant in a Protein Library  
JF  - bioRxiv  
DO  - 10.1101/2021.11.18.469179  
SP  - 2021.11.18.469179  
AU  - Wittmann, Bruce J.  
AU  - Johnston, Kadina E.  
AU  - Almhjell, Patrick J.  
AU  - Arnold, Frances H.  
Y1  - 2021/01/01  
UR  - http://biorxiv.org/content/early/2021/11/19/2021.11.18.469179.abstract  
N2  - Widespread availability of protein sequence-fitness data would revolutionize both our biochemical understanding of proteins and our ability to engineer them. Unfortunately, even though thousands of protein variants are generated and evaluated for fitness during a typical protein engineering campaign, most are never sequenced, leaving a wealth of potential sequence-fitness information untapped. This largely stems from the fact that sequencing is unnecessary for many protein engineering strategies; the added cost and effort of sequencing is thus unjustified. Here, we present every variant sequencing (evSeq), an efficient protocol for sequencing a variable region within every variant gene produced during a protein engineering campaign at a cost of cents per variant. Execution of evSeq is simple, requires no sequencing experience to perform, relies only on resources and services typically available to biology labs, and slots neatly into existing protein engineering workflows. Analysis of evSeq data is likewise made simple by its accompanying software (found at github.com/fhalab/evSeq, documentation at fhalab.github.io/evSeq), which can be run on a personal laptop and was designed to be accessible to users with no computational experience. Low-cost and easy to use, evSeq makes collection of extensive protein variant sequence-fitness data practical.
```

</details>

### Read the Docs!
Find detailed documentation at the [individual pages linked below](#documentation) or [start at the overview](0-theory.md).

#### Quick links to common resources:

##### Biology
- [Inner Primer Design](1-lib_prep.md#inner-primer-design)
- [PCR Protocol](1-lib_prep.md#pcr-protocol)
- [PCR Product Purification](1-lib_prep.md#pcr-product-purification)

##### Computation
- [Installation](3-installation.md)
- [The `refseq` file](4-usage.md#the-refseq-file)
- [`OutputCounts`](5-outputs.html#OutputCounts)
- [Running `evSeq` in a Jupyter Notebook](8-full_demo.html)

##### [Troubleshooting](9-troubleshooting.md)

### General Overview
#### The `evSeq` workflow
![Workflow](assets/figure1.png)
**A)** `evSeq` amplifies out a region of interest that contains variability, attaches well-specific barcodes and adapters, and is ready for NGS.

**B)** All that's required to perform the `evSeq` laboratory procedure is:

- a 96-well thermalcycler
- standard PCR reagents and materials
- access to an NGS provider
- two 96-well plates of `evSeq` barcoding ("outer") primers
- a pair of region-specific `evSeq`-compatible ("inner") primers
- 96-well plate(s) of cultures containing DNA encoding protein variants
- a 12-channel 10 µL pipette is also helpful

**That's it.**

Due to the two-primer, culture-based PCR methodology employed by `evSeq`, only a new pair of inner primers needs to be ordered when targeting new regions/sequences and no DNA isolation needs to be performed.

**C)** Once the sequences are returned by the NGS provider, the computational workup can be performed on a standard laptop by users with little-to-no computational experience.

The amplicons prepared with `evSeq` can yield nearly 1000 high-quality protein variant sequences for the just cost of the multiplexed NGS run (typically ~$100 from commercial sequencing providers, likely lower for in-house providers).

#### Construct and visualize sequence-function pairs
![SeqFunc](assets/figure2.png)
Sequencing eight site-saturation libraries (768 wells) in a single `evSeq` run and combining this with activity data to create low-cost sequence-function data. **A)** Enzyme and active-site structure highlighting mutated residues. **B)** Heatmap of the number of identified variants/mutations ("counts") for each position mutated ("library") from processed `evSeq` data. **C)** Heatmap of the average activity ("normalized rate") for each variant/mutation in each library. **D)** Counts for a single library, also showing the number of unidentified wells. **E)** Activity for a single library, showing biological replicates. (Inset displays the mutated residue in this library.)

## Documentation
### Biology
#### [Theoretical overview](0-theory.md)

#### [Library preparation](1-lib_prep.md)
- [Dual-Index Barcode Plates](1-lib_prep.md#dual-index-barcode-plates)
- [Inner Primer Design](1-lib_prep.md#inner-primer-design)
- [Inner Primer Test PCR](1-lib_prep.md#inner-primer-test-pcr)
- [PCR Protocol](1-lib_prep.md#pcr-protocol)
- [PCR Product Purification](1-lib_prep.md#pcr-product-purification)

### Computation
#### [Computational basics](2-basics.md)
#### [Installation](3-installation.md)
- [Installing from GitHub with the `conda` environment](3-installation.md#installing-from-github-with-the-conda-environment)
- [Standard `pip` Install and Dependencies](3-installation.md#standard-pip-install-and-dependencies)

#### [Running `evSeq`](4-usage.md)
- [Post Installation](4-usage.md#post-installation)
- [Using `evSeq` from the command line or GUI](4-usage.md#using-evseq-from-the-command-line-or-gui)
- [Required Arguments](4-usage.md#required-arguments)
  - [The `refseq` file](4-usage.md#the-refseq-file)
  - [`folder`](4-usage.md#folder)
- [Optional Arguments](4-usage.md#optional-arguments)

#### [Understanding the Outputs](5-outputs.html)
- [`Qualities`](5-outputs.html#Qualities)
- [`OutputCounts`](5-outputs.html#OutputCounts)
- [`Platemaps`](5-outputs.html#Platemaps)
- [`evSeqLog` files](5-outputs.html#evSeqLog-files)
- [`ParsedFilteredFastqs`](5-outputs.html#ParsedFilteredFastqs)
- [`Alignments`](5-outputs.html#Alignments)

### Additional Examples
Below are a collection of Jupyter Notebooks (rendered as documents) with examples on how to get the most out of `evSeq`. If you want to run them on your own, they can be found in the [examples](https://github.com/fhalab/evSeq/blob/master/examples/) directory of the `evSeq` repository.
#### [Using `evSeq` data](6-using_evseq_data.html)
- [Importing and viewing `evSeq` data](6-using_evseq_data.html#Importing-and-viewing-evSeq-data)
- [Pairing sequence to function](6-using_evseq_data.html#Pairing-sequence-to-function)
- [Analyzing single-site-saturation libraries](6-using_evseq_data.html#Analyzing-single-site-saturation-libraries)
- [`evSeq` for multisite libraries](6-using_evseq_data.html#evSeq-for-multisite-libraries)
- [Generating our visualizations](6-using_evseq_data.html#Generating-our-visualizations)
- [Submitting to Protαβank and other databases](6-using_evseq_data.html#Submitting-to-Protαβank-and-other-databases)
- [Miscellaneous visualizations](6-using_evseq_data.html#Miscellaneous-visualizations)

#### [Creating barcode/index pairs](7-index_mapping.html)
- [Using new barcode primers](7-index_mapping.html#Using-new-barcode-primers)
- [Creating new index pair mappings](7-index_mapping.html#Creating-new-index-pair-mappings)

#### [Running `evSeq` in a Jupyter Notebook](8-full_demo.html)
- To run this notebook on your own, open it in the `evSeq` repository and run it from its current location (found as [evSeq/examples/8-full_demo.ipynb](https://github.com/fhalab/evSeq/blob/master/examples/8-full_demo.ipynb)).

### Troubleshooting
- [Poor reverse read quality](9-troubleshooting.md#poor-reverse-read-quality)
- [Poor results but good quality sequencing](9-troubleshooting.md#poor-results-but-good-quality-sequencing)
- [Windows: `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'`](9-troubleshooting.md#windows-commandnotfounderror-your-shell-has-not-been-properly-configured-to-use-conda-activate)
- [macOS: `PermissionError: [Errno 1] Operation not permitted`](9-troubleshooting.md#macos-permissionerror-errno-1-operation-not-permitted)
- [Linux: `CondaEnvException: Pip failed`](9-troubleshooting.md#linux-condaenvexception-pip-failed)
