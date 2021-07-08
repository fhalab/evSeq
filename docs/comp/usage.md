# Running `evSeq`

## Post Installation
### Confirming your installation
[shell script that runs lots of things]
### Running example files
[refer to GUI and command line or in a python environment sections, and use [path_to_file] to test.]

## The `refseq` file
The primary user inputs that are required are contained in the `refseq` file, which contains information that allows the `evSeq` software to know how to process each well.
### Default `refseq`
This form of the file assumes the same reference sequence in each well of the analyzed plates and requires eight columns: `PlateName`, `IndexPlate`, and `FPrimer`, `RPrimer`, `VariableRegion`, `FrameDistance`, `BpIndStart`, and `AaIndStart`. These columns are detailed below:

| Column | Description |
|:-------|-------------|
| `PlateName` | This is a nickname given to the plate. For instance, if you performed `evSeq` on a plate that you called "Plate1", you would put "Plate1" in this column. |
| `IndexPlate` | This is the `evSeq` index plate used for library preparation corresponding to the plate in `PlateName`. For instance, if I prepared "Plate1" using index plate 2, I would write "DI02" in the `IndexPlate` column. Allowed barcode names are `DI01` through `DI08`. |
| `FPrimer` | This is the sequence-specific forward primer you used to create the amplicon for attaching `evSeq` barcodes, including the `evSeq` adapter regions. It should be input exactly as ordered from your oligo supplier, 5' - 3'. |
| `RPrimer` | This is the sequence-specific reverse primer you used to create the amplicon for attaching `evSeq` barcodes, including the `evSeq` adapter regions. It should be input exactly as ordered from your oligo supplier, 5' - 3'. Do *not* use the reverse complement. |
| `VariableRegion` | This is the entire region between the 3' ends of each primer used to generate the `evSeq` amplicon, 5' - 3'. This is called the "variable region" as it is the sequenced region that can vary meaningfully, since the primers should be invariant. If a specific codon has been specifically mutagenized, this codon may be replaced by "NNN". See below for more details. |
| `FrameDistance` | For a given `evSeq` run, 2 in 3 times the sequence used will be out of reading frame with the full amplicon. **To allow for accurate translation of sequences, the distance to the first base in the `VariableRegion` that is in the reading frame of the full gene must be provided.**  For example, if the 3' end of your `FPrimer` ends on the last base of a codon, your `VariableRegion` is in-frame and this is argument should be `0`. If `FPrimer` ends on the second base of a codon (e.g., is shifted back 1 bp), then your first in-frame base is 1 base away, so this argument should be `1`. Note that if "NNN" is used in the `VariableRegion`, `evSeq` will double check that you correctly defined this argument — with no "NNN" it will assume you are correct. |
| `BpIndStart` | This argument tells the program what index the first base in the variable region belongs to. This is useful for formatting the outputs, as any variation identified in `evSeq` can be output at the index corresponding to the full gene, rather than the amplicon. |
| `AaIndStart` | This argument tells the program what index the first *in-frame* amino acid in the variable region belongs to. This means that if your `FrameDistance` argument is not `0`, `AaIndStart` should not be the position of the codon your variable region starts in but rather the next one, since the first codon is not in frame. |

As currently deployed, up to 8 plates (`DI01`–`DI08`) can be input in a single `evSeq` run. No more than 8 rows should thus ever be filled in this form of `refseq` file. An example Default `refseq` format is given in the `evSeq` GitHub repository [here](../../examples/refseqs/DefaultRefSeqs.csv).

#### `VariableRegion` containing "NNN"
(OPTIONAL) You may replace the bases at the known mutagenized positions with "NNN" as the codon. Doing so forces `evSeq` to return the sequence identified at these positions (e.g., from a site-saturation mutagenesis library), whether or not it matches the parent. If you know where your mutations will occur, this is the recommended way to use `evSeq`; any off-target mutations not given by "NNN" will still be identified and reported.

### Detailed `refseq`
This form of the file allows for a different reference sequence in each well of the analyzed plates. In addition to the column headers given in [Default `refseq`](#default-refseq), this form of the file has a required `Well` column, enabling specification of a different `FPrimer`, `RPrimer`, and `VariableRegion` for each well in the input plates. As currently deployed, up to 8 plates (`DI01`–`DI08`) can be input in a single `evSeq` run, so no more than 768 rows should ever be filled in this form of `refseq` file. An example Detailed `refseq` format is given in the `evSeq` GitHub repository [here](../../examples/refseqs/DetailedRefSeqs.csv).

**When using this form of `refseq`, the `detailed_refseq` flag must be set** (see next sections for details).
## Using `evSeq` from the command line or GUI

## Using `evSeq` in a Python environment


---

*Back to the [main page](../index.md).*
