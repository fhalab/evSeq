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

| Column | Data Type | Description |
|:-------|:----------|-------------|
| `PlateName` | `str` | This is a nickname given to the plate. For instance, if you performed `evSeq` on a plate that you called "Plate1", you would put "Plate1" in this column. |
| `IndexPlate` | `DI0X`, `X=[1,8]` | This is the `evSeq` index plate used for library preparation corresponding to the plate in `PlateName`. For instance, if "Plate1" weere prepared using index plate 2, `IndexPlate` would be `DI02`. Allowed barcode names are `DI01` through `DI08`. |
| `FPrimer` | `str`, DNA Sequence | This is the sequence-specific forward primer you used to create the amplicon for attaching `evSeq` barcodes, including the `evSeq` adapter regions. It should be input exactly as ordered from your oligo supplier, 5' - 3'. |
| `RPrimer` | `str`, DNA Sequence | This is the sequence-specific reverse primer you used to create the amplicon for attaching `evSeq` barcodes, including the `evSeq` adapter regions. It should be input exactly as ordered from your oligo supplier, 5' - 3'. Do *not* use the reverse complement. |
| `VariableRegion` | `str`, DNA Sequence | This is the entire region between the 3' ends of each primer used to generate the `evSeq` amplicon, 5' - 3'. This is called the "variable region" as it is the sequenced region that can vary meaningfully, since the primers should be invariant. If a specific codon has been specifically mutagenized, this codon may be replaced by "NNN". See below for more details. |
| `FrameDistance` | `int`, `[0,1,2]` | Distance (in base pairs) from the 3' end of `FPrimer` to the first in-frame codon in your `VariableRegion`. **This is required for accurate translation of sequences.** For a given `evSeq` run, 2 in 3 times the sequence used will be out of reading frame with the full amplicon, so this is important to check. As an example, if the 3' end of your `FPrimer` ends on the last base of a codon, your `VariableRegion` is in-frame and this is argument should be `0`. If `FPrimer` ends on the second base of a codon (e.g., is shifted back 1 bp), then your first in-frame base is 1 base away and this argument should be `1`. Note that if "NNN" is used in the `VariableRegion`, `evSeq` will double check that you correctly defined this argument — with no "NNN" it will assume you are correct. |
| `BpIndStart` | `int`, `[0,inf]` | This argument tells the program what index the first base in the variable region belongs to. This is useful for formatting the outputs, as any variation identified in `evSeq` can be output at the index corresponding to the full gene, rather than the amplicon. |
| `AaIndStart` | `int`, `[0,inf]` | This argument tells the program what index the first *in-frame* amino acid in the variable region belongs to. This means that if your `FrameDistance` argument is not `0`, `AaIndStart` should not be the position of the codon your variable region starts in but rather the next one, since the first codon is not in frame. |

#### Example
```
   -------FPrimer------>
5'-AAAAAAAAAAGGGGGGGGGGG-3'
             |||||||||||--------VariableRegion------->
5'-CCCCCCCCCCGGGGGGGGGGGTNNNTTTTTTTT...TTTTTTTTTTTTTTTGGGGGGGGGGGGCCCCCCCCCC-3'
   FrameDistance = 1 -> x                             ||||||||||||
                                                   3'-CCCCCCCCCCCCAAAAAAAAAA-5'
                                                      <-------RPrimer-------
```
In this simple example, the `FPrimer` sequence is `AAAAAAAAAAGGGGGGGGGGG` and the `RPrimer` sequence is `AAAAAAAAAACCCCCCCCCCCC`, with the `A` regions corresponding to the `evSeq` adapter regions and the `G/C` regions corresponding to the directional sequence-specific regions. The `VariableRegion` is the area between them, `TNNNTTTTTTTT...TTTTTTTTTTTTTTT`. The `NNN` sequence specifies the correct reading frame, and the 3' base of `FPrimer` is 1 base away from this in-frame codon (indicated with the `x`), therefore `FrameDistance` is `1`. If the `VariableRegion` started as base 35 in the gene, `BpIndStart` would be 35 and `AaIndStart` would correspond to the position of the `NNN` codon, which would be position 12 (the first base of this codon is base 36, which is the 12th codon overall).

(Note that if any of your sequences look anything like this example, your library preparation will not work. This is to be interpreted as a simple example only.)

#### `VariableRegion` containing "NNN"
(OPTIONAL) You may replace the bases at the known mutagenized positions with "NNN" as the codon. Doing so forces `evSeq` to return the sequence identified at these positions (e.g., from a site-saturation mutagenesis library), whether or not it matches the parent. If you know where your mutations will occur, this is the recommended way to use `evSeq`; any off-target mutations not given by "NNN" will still be identified and reported.

#### `evSeq` Index Plates
As currently deployed, up to 8 plates (`DI01`–`DI08`) can be input in a single `evSeq` run. No more than 8 rows should thus ever be filled in this form of `refseq` file. An example Default `refseq` format is given in the `evSeq` GitHub repository [here](../../examples/refseqs/DefaultRefSeqs.csv).

### Detailed `refseq`
This form of the file allows for a different reference sequence in each well of the analyzed plates. In addition to the column headers given in [Default `refseq`](#default-refseq), this form of the file has a required `Well` column, enabling specification of a different `FPrimer`, `RPrimer`, and `VariableRegion` for each well in the input plates. As currently deployed, up to 8 plates (`DI01`–`DI08`) can be input in a single `evSeq` run, so no more than 768 rows should ever be filled in this form of `refseq` file. An example Detailed `refseq` format is given in the `evSeq` GitHub repository [here](../../examples/refseqs/DetailedRefSeqs.csv).

**When using this form of `refseq`, the `detailed_refseq` flag must be set** (see next sections for details).
## Using `evSeq` from the command line or GUI
### Command line
Thanks to `setuptools` `entry_points`, the main functionality of `evSeq` can be accessed from the command line as if it were added to `PATH` by running:
```
evSeq refseq data_dir --OPTIONAL_ARGS [ARG_VALUE] --FLAGS
```
where `evSeq` has been set up in PATH, `refseq` is the .csv file containing information about the experiment described above, and `data_dir` is the directory that contains the raw .fastq files for the experiment.

For information on optional arguments and flags, run
```
evSeq -h
```
or see [below](#optional-arguments).

### GUI
# Describe how to launch; fix with `py2app` and `py2exe`?
[need to fix it for Mac still?]

Once opened, you will see two required arguments — the `refseq` and `folder` args — at the top of the GUI. Details on the `refseq` argument is described [above](#the-refseq-file), and the GUI should provide a description of what the `folder` contains. For more advanced use, other arguments can be accessed by scrolling down. (These additional arguments are detailed in [Optional Arguments](#optional-arguments)). You will typically not need these arguments, however, and the standard `evSeq` run can be started by clicking `Start` once `refseq` and `folder` are populated. Once started, the progress of the program will be printed to the GUI along with any encountered warnings and errors.

## Optional Arguments
There are a number of flags and optional arguments that can be passed for `evSeq`, all detailed in the table below:

| Argument | Type | Description |
|:---------|------|-------------|
| `fastq_r` | Argument | This argument is only available for command line use. If a case arises where, for whatever reason, `evSeq` cannot auto-identify the forward and reverse read files, this option acts as a failsafe. Instead of passing the folder containing the forward and reverse files in to the `folder` required argument, pass in the forward read file as the `folder` argument and the reverse read file as this optional argument. |
| `output` | Argument | By default, `deSeq` will save to the current working directory (command line) or the `evSeq` Git repository folder (GUI). The default save location can be overwritten with this argument. |
| `detailed_refseq` | Flag | Set this flag (check the box in the GUI) when passing in a detailed reference sequence file. See [Detailed refseq](#detailed-refseq) for more information. |
| `analysis_only` | Flag | Set this flag (check the box in the GUI) to only perform Q-score analysis on the input fastq files. The only output in this case will be the [quality score histograms](outputs.md#qualities).|
| `only_parse_fastqs` | Flag | Set this flag to stop `evSeq` after generation of parsed, well-filtered fastq files. Counts, platemaps, and alignments will not be returned in this case. Used in case the well-specific fastq sequences are desired but not the entire `evSeq` analysis. |
| `keep_parsed_fastqs` | Flag | Set this flag to save parsed, well-filtered fastq files as in `only_parse_fastqs`  but to also finish the regular `evSeq` run. |
| `return_alignments` | Flag | Set this flag to return alignments along with the `evSeq` output. Note that this flag is ignored if either `analysis_only` or `stop_after_fastq` are used. |
| `average_q_cutoff` | Argument | During initial sequencing QC, `evSeq` will discard any sequence with an average quality score below this value. The default value is 25. |
| `bp_q_cutoff` | Argument | Bases with a q-score below this value are ignored when counting the number of sequences aligned at each position. For the coupled outputs (see below), counts are only returned if all bases in the combination pass. The default value is 30. |
| `length_cutoff` | Argument | During initial sequencing QC, `evSeq` will discard any sequence with an read with total length below `length_cutoff * read_length`. The default value is `0.9`. |
| `variable_thresh` | Argument | This argument sets the threshold that determines whether or not a position is variable. In other words, if a position contains a non-reference sequence sequence at a given position at a fraction greater than `variable_thresh`, then it is a variable position. The default is `0.2`. Setting this value lower makes `evSeq` more sensitive to variation, while setting it higher makes it less sensitive. A value of 1, for instance, would find no variable positions. |
| `variable_count` | Argument | This sets the count threshold for identifying "dead" wells. If a well has fewer sequences that pass QC than this value, then it is considered "dead". The default value is 10 (meaning only wells with fewer than 10 sequences are dead). |
| `jobs` | Argument | This is the number of processors used by deSeq for data processing. By default, `evSeq` uses 1 less processor than are available on your computer. As with all multiprocessing programs, it is typically not recommended to use all available processors unless you are okay devoting all computer resources to the task (e.g. you don't want to be concurrently checking email, playing music, running another program, etc.). The number of jobs can be lowered to reduce the memory demands of `evSeq`. |
| `read_length` | Argument | By default, `evSeq` will attempt to determine the read length from the fastq files. If this process is failing (e.g., due to heavy primer-dimer contamination), the read length can be manually set using this argument. |

---
*Next page: [Understanding the outputs](outputs.md).*

*Back to the [main page](../index.md).*
