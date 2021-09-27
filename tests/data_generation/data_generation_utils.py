# Import 3rd party
from Bio.Seq import Seq

# Import globals module
import tests.data_generation.globals as test_glob

# Import specific variables needed from globals that WILL NEVER BE UPDATED
from tests.data_generation.globals import MAX_QUAL_ALLOWED, Q_SCORE_BASE

class QualityGenerator():
    """
    Utility for generating quality score arrays
    """
    def __init__(self, min_q_allowed):
        self.min_q_allowed = min_q_allowed
        
    def generate_qualities(self, size):
        return test_glob.NP_RNG.integers(self.min_q_allowed, MAX_QUAL_ALLOWED,
                                         size = size)
        
# Use qualities as an ordinal encoding
def ord_to_chr(qualities):
    rebased_qs = qualities + Q_SCORE_BASE
    return "".join(chr(qual) for qual in rebased_qs)

                
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())