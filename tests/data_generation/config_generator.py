# Import globals module
import tests.data_generation.globals as test_glob

# Import testing globals that never change
from tests.data_generation.globals import (
    N_INDICES, N_PLATES, MAX_READLENGTH, MIN_READLENGTH, MIN_GLOBAL_QUAL_CUTOFF,
    MAX_GLOBAL_QUAL_CUTOFF, MIN_BP_QUAL_CUTOFF, MAX_BP_QUAL_CUTOFF, 
    MIN_SEQLEN_CUTOFF, MAX_SEQLEN_CUTOFF, MAX_VARIABLE_THRESH, 
    MIN_VARIABLE_THRESH, MIN_VARIABLE_COUNT, MAX_VARIABLE_COUNT
)

# Import refseq generator
from .refseq_generator import FakeRefseq

# Class that holds parameters that will be needed by evSeq
class Config():
    def __init__(self, detailed = True):
        """
        Builds a set of conditions that might be passed into an evSeq run.
        """
        # Record whether or not this is a detailed run
        self.detailed = detailed
        
        # Build as many reference sequences as needed
        n_refs_needed = N_INDICES if detailed else N_PLATES
        self.refseqs = [FakeRefseq() for _ in range(n_refs_needed)] 
        
        # Decide on the readlength that we will be using. The maximum allowed
        # readlength cannot be longer than the full sequencable length of
        # the shortest a reference sequence
        min_refseq_len = min(refseq.readable_window_len for refseq in self.refseqs)
        max_readlength = min(min_refseq_len, MAX_READLENGTH)
        self.readlength = test_glob.NP_RNG.integers(MIN_READLENGTH, max_readlength)
                
        # Decide on the input variables that will define the run, including 
        # average_q_cutoff, bp_q_cutoff, length_cutoff, variable_thresh, and
        # variable_count
        self.bp_q_cutoff = test_glob.NP_RNG.integers(MIN_BP_QUAL_CUTOFF, MAX_BP_QUAL_CUTOFF)
        max_global_qual_cutoff = max(self.bp_q_cutoff, MAX_GLOBAL_QUAL_CUTOFF)
        self.average_q_cutoff = test_glob.NP_RNG.integers(MIN_GLOBAL_QUAL_CUTOFF, max_global_qual_cutoff)
        self.length_cutoff = test_glob.NP_RNG.uniform(MIN_SEQLEN_CUTOFF, MAX_SEQLEN_CUTOFF)
        self.variable_thresh = test_glob.NP_RNG.uniform(MIN_VARIABLE_THRESH, MAX_VARIABLE_THRESH)
        self.variable_count = test_glob.NP_RNG.integers(MIN_VARIABLE_COUNT, MAX_VARIABLE_COUNT)
        
        # Build base quality scores for the reference sequences. Also define the allowed
        # mutagenesis windows for the amino acids, which encompasses the last amino acids
        # captured in full by the readlength. 
        for refseq in self.refseqs:
            refseq.assign_qualities(max(self.average_q_cutoff, self.bp_q_cutoff) + 1)
            refseq.define_refseq_windows(self.readlength)