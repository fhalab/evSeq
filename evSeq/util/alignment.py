from Bio import pairwise2
import warnings

# Redefine the BioPython aligment function so that we can quicky change
# parameters later
def evseq_align(reference, query, 
                match = None, mismatch = None,
                open_penalty = None, extend = None):
    
    # Redefine biopython aligment function (this is just for code neatness)
    return pairwise2.align.globalms(reference, query, 
                                    match, mismatch,
                                    open_penalty, extend,
                                    one_alignment_only=True, 
                                    penalize_end_gaps = False)[0]