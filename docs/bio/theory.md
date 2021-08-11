# Theoretical Overview
![Overview](../assets/figure2.png)

`evSeq` is a sub-barcoding procedure designed to make the data in a routine multiplexed NGS experiment go further. This is useful for situations that only call for hundreds of reads to confidently describe protein variant populations in a single well, where the hundreds of thousands of reads from the typical experiment are unnecessary.

In a typical commericial multiplexed NGS experiment, amplicons from each researcher are submitted with specific adapters that allow the attachment of barcodes. Different barcodes are attached to different researchers' samples, which then allows the company running the NGS experiment to pool the samples, sequence them together, and then de-multiplex the sequencing data to send each researcher their (and only their) data. For a typical MiSeq NGS run, this spreads the typical ~10 million reads into 96 groups of ~100,000 reads, as the full ~10 million reads is not as frequently necessary.

For protein engineering, protein variants that are arrayed into 96-well plates are typically (and ideally) monoclonal, meaning that minimal sequencing information is required to determine the identity of the variant in each well. In a world with perfect sequencing and perfectly monoclonal variants, only a single read would be required. These steps are imperfect though, so more than one read is necessary. However, ~100,000 reads is still very much overkill. Therefore, *`evSeq` was created to further spread these reads over ~100–1000 protein variants*, returning enough reads per well for the confident assessment of even polyclonal cultures.

`evSeq` provides hundreds of protein sequences for the same price as a handful of Sanger sequencing runs and can fully identify each variant in polyclonal cultures. (This is not always possible with Sanger sequencing, which just gives bulk information on the population and not individual DNA molecules.) See the page on [library preparation](lib_prep.md) for details on how to use `evSeq` in the laboratory, and the [Computation](../index.md#computation) documentation on how to use the data.

---

*Next page: [Library Preparation](lib_prep.md).*

*Back to the [main page](../index.md).*