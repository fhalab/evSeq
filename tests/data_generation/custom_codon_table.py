# Import 3rd party modules
from Bio.Data import CodonTable

class CustomCodonTable():
    """
    Because the BioPython implementation of a codon table is all over the place,
    fixing all the flaws in this class. Not inheriting because I just don't trust it.
    """
    def __init__(self):
        
        # Get the standard codon table from BioPython. Their docs do not say what this is,
        # but it looks right from a spot check. As it is, we wrote our own for evSeq, so it
        # is good to have a separate codon table used for validation.
        bio_codon_table = CodonTable.standard_dna_table
        
        # Build the forward dictionary. Add the stop codon translations (why
        # aren't these included in the first place???).
        self.codon_to_aa = bio_codon_table.forward_table.copy()
        for stop_codon in bio_codon_table.stop_codons:
            self.codon_to_aa[stop_codon] = "*"
            
        # Build the reverse dictionary. BioPython reverse dictionary includes just one 
        # codon per aa. Fix this.
        self.aa_to_codon = {}
        for codon, aa in self.codon_to_aa.items():
            if aa in self.aa_to_codon:
                self.aa_to_codon[aa].append(codon)
            else:
                self.aa_to_codon[aa] = [codon]
                
# Instantiate a default codon table and allowed amino acid characters
CODON_TABLE = CustomCodonTable()
ALLOWED_AAS = tuple(sorted(list(CODON_TABLE.aa_to_codon.keys())))
INT_TO_AA = dict(enumerate(ALLOWED_AAS))
N_AAS = len(ALLOWED_AAS)