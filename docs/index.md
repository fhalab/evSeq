## Every Variant Sequencing: No sequence-function pair left behind.

Every Variant Sequencing (`evSeq`) provides resources for extremely low cost massively parallel sequencing of protein/enzyme variants arrayed in 96-well plates. This library preparation technique and data processing tool enables sequencing all variants produced during a protein engineering or biochemical mutagenesis experiment at a cost of cents per variant, even for labs that do not have expertise in or access to next-generation sequencing (NGS) technology.

### [Read the Paper!](one_day...)

### Read the Docs!
Navigate to [individual pages below](#documentation) or [start at the overview](bio/theory.md).

### General Overview
#### Workflow
![Workflow](assets/overview.png)
#### Example Output (1 of 8 plates)
![Output](assets/platemap.png)

#### Sequence-Function Data
## fig 3

## Documentation
### Biology
#### [Theoretical overview](bio/theory.md)

#### [Library preparation](bio/lib_prep.md)
- [Inner Primer Design](bio/lib_prep.md#inner-primer-design)
- [PCR Protocol](bio/lib_prep.md#pcr-protocol)
- [PCR Product Purification](bio/lib_prep.md#pcr-product-purification)

### Computation
#### [Computational basics](comp/basics.md)
#### [Installation](comp/installation.md)
- [Installing from GitHub](comp/installation.md#installing-from-github)
- [Installing with PyPI (pip)](comp/installation.md#installing-from-pypi)
- [Using the `evSeq` environment](comp/installation.md#using-the-evseq-environment)
#### [Running `evSeq`](comp/usage.md)
- [Post Installation](comp/usage.md#post-installation)
- [The `refseq` file](comp/usage.md#the-refseq-file)
- [Using `evSeq` from the command line or GUI](comp/usage.md#using-evseq-from-the-command-line-or-gui)
- [Optional Arguments](comp/usage.md#optional-arguments)
#### [Understanding the outputs](comp/outputs.md)
- [Qualities](comp/outputs.md#qualities)
- [OutputCounts](comp/outputs.md#outputcounts)
- [Platemaps](comp/outputs.md#platemaps)


#### [Using `evSeq` tools in a Python environment](comp/additional.html)
- [Running `evSeq` functions]()
- [Mapping sequence to function]()

### Troubleshooting
- [Poor reverse read quality](troubleshooting.md#poor-reverse-read-quality)
- [Poor alignments but good quality](troubleshooting.md#poor-alignments-but-good-quality)
