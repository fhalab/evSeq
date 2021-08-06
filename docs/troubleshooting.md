# Troubleshooting
## Poor reverse read quality
If your forward read quality is excellent but your reverse read quality is very poor, this can ruin your `Coupled` output files. While you can still get useful information from the forward reads in the `Decoupled` output files, this issue usually comes from library preparation and should be fixed as soon as possible.

This can occur due to certain primer designs. While the `evSeq` adapters have been optimized (both computationally and through trial-and-error) to work with most inner primer sequences, we can't expect this to always be perfect. If possible, re-design your inner primers in slightly new locations (re-designing both and testing all combinations is usually a safe bet).

This also occurs due to too much primer-dimer in your submitted sample. This is a subset of the above problem (better primer design means less primer-dimer) but can be resolved by running a longer agarose gel to separate the primer-dimer (~75-150 bp) from the actual barcoded `evSeq` amplicon. This is especially important when the amplicon is small.
## Poor alignments/platemap but good quality
If your plate seems to contain good quality reads but mostly dead wells, and your average read length in the fastqs is not unreasonably low, it is likely that your `refseq` file information is not correct for the given sequence. Often this comes from having the wrong sequence information (`FPrimer/RPrimer`s and/or `VariableRegion`, setting the alignment process off.

#### Tips to solve
- Ensure that your reference sequence information accurately reflects the sequences you sent your NGS provider. Small mistakes in the reference frame can ruin the analysis.
- Run `evSeq` with the `--keep_parsed_fastqs` or `--only_parse_fastqs` flag to confirm that they match your reference sequence information.
- Run `evSeq` with the `--return_alignments` flag to see where the alignment of each read to your reference sequence is failing.

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

## macOS: `PermissionError: [Errno 1] Operation not permitted`
Any time after upgrading to a newer macOS, you might randomly find that `evSeq` has stopped working and gives `PermissionError: [Errno 1] Operation not permitted`.

This is a security issue for software being run on macOS, and can be solved by following [these directions](https://stackoverflow.com/questions/58479686/permissionerror-errno-1-operation-not-permitted-after-macos-catalina-update) in your System Preferences/Security and Privacy settings.

---

*Back to the [main page](../index.md).*