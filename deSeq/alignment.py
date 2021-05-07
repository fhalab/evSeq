from Bio import pairwise2

# Redefine the BioPython aligment function so that we can quicky change
# parameters later
def deseq_align(reference, query):
    
    # Redefine biopython aligment function (this is just for code neatness)
    return pairwise2.align.globalxs(reference, query, open = -2, extend = -1,
                                    one_alignment_only=True, penalize_end_gaps = False)[0]