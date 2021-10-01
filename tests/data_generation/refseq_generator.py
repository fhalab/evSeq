# Import globals module
import tests.data_generation.globals as test_glob

# Import testing globals that never change
from tests.data_generation.globals import (
    MIN_REFSEQ_LEN, MAX_REFSEQ_LEN, FRAMESHIFT_MIN, FRAMESHIFT_MAX,
    ALLOWED_NUCLEOTIDES, PRIMER_MAX_LEN, PRIMER_MIN_LEN, AA_IND_START_MIN,
    AA_IND_START_MAX, BP_IND_START_MAX, BP_IND_START_MIN
)
from tests.data_generation.custom_codon_table import (
    N_AAS, INT_TO_AA, CODON_TABLE
)

# Import evSeq globals needed for testing
from evSeq.util.globals import (
    ADAPTER_F, ADAPTER_LENGTH_F, ADAPTER_LENGTH_R, ADAPTER_R, BARCODE_LENGTH
)

# Import utils
from .data_generation_utils import QualityGenerator, reverse_complement

# Import 3rd party modules
import numpy as np

class FakeRefseq():
    def __init__(self):
        
        # Randomly create an amino acid reference sequence. Record the length.
        self.refseq_len = test_glob.NP_RNG.integers(MIN_REFSEQ_LEN, MAX_REFSEQ_LEN)
        aa_ints = test_glob.NP_RNG.choice(N_AAS, size = self.refseq_len)
        self.aa_refseq = [INT_TO_AA[aa_int] for aa_int in aa_ints]

        # Assign a codon to each amino acid.
        self.codon_refseq = [test_glob.RANDOM_RNG.choice(CODON_TABLE.aa_to_codon[aa]) 
                             for i, aa in enumerate(self.aa_refseq)] 
        self.codon_refseq_len = self.refseq_len * 3
        
        # Decide on the frameshift of the reference sequence
        self.frameshift_front = test_glob.NP_RNG.integers(FRAMESHIFT_MIN, FRAMESHIFT_MAX)
        self.frameshift_back = test_glob.NP_RNG.integers(FRAMESHIFT_MIN, FRAMESHIFT_MAX)        
        
        # Determine the bases for the frameshift
        self.frameshift_bp_front = test_glob.RANDOM_RNG.choices(ALLOWED_NUCLEOTIDES, 
                                                      k = self.frameshift_front)
        self.frameshift_bp_back = test_glob.RANDOM_RNG.choices(ALLOWED_NUCLEOTIDES, 
                                                     k = self.frameshift_back)
                
        # Create primer seeds
        self.primer_seed_len_f = test_glob.NP_RNG.integers(PRIMER_MIN_LEN, PRIMER_MAX_LEN)
        self.primer_seed_len_r = test_glob.NP_RNG.integers(PRIMER_MIN_LEN, PRIMER_MAX_LEN)
        self.primer_seed_f = test_glob.RANDOM_RNG.choices(ALLOWED_NUCLEOTIDES, k = self.primer_seed_len_f)
        self.primer_seed_r = test_glob.RANDOM_RNG.choices(ALLOWED_NUCLEOTIDES, k = self.primer_seed_len_r)
        
        # Assign the total readable window
        self.readable_window_len = (self.codon_refseq_len + 
                                    self.frameshift_front +
                                    self.frameshift_back + 
                                    self.primer_seed_len_f + 
                                    self.primer_seed_len_r +
                                    ADAPTER_LENGTH_F +
                                    ADAPTER_LENGTH_R +
                                    2 * BARCODE_LENGTH)
        
        # Set the bp and aa ind starts
        self.bp_ind_start = test_glob.NP_RNG.integers(BP_IND_START_MIN, BP_IND_START_MAX)
        self.aa_ind_start = test_glob.NP_RNG.integers(AA_IND_START_MIN, AA_IND_START_MAX)
    
    def define_refseq_windows(self, readlength):
        """
        Defines the region within the mutable refseq region that can be modified.
        """
        # Choose legal positions for making mutations. For forward, the first readable comes after
        # the barcode, adapter, and seed region; this is the last readable for reverse. For forward,
        # the last readable is the last amino acid  captured in full by the readlength, considering
        # the primer binding region, the adapter machinery, and the barcode.
        effective_readlength = readlength - BARCODE_LENGTH
        readable_f_window = (effective_readlength - self.frameshift_front -
                             ADAPTER_LENGTH_F - self.primer_seed_len_f) // 3
        readable_r_window = ((effective_readlength - self.frameshift_back - 
                              ADAPTER_LENGTH_R - self.primer_seed_len_r) // 3)
        
        # Make sure we don't break the bounds of the readable region. We don't
        # allow mutations at the first and last positions during testing (this is
        # because insertions end up randomly placed as the aligned prefs to add
        # a base or two to get better alignment)
        self.max_readable_aa_ind_f = min(readable_f_window, self.refseq_len)
        self.min_readable_aa_ind_r = max(self.refseq_len - readable_r_window, 0)
        
        # Define the readable positions in each direction
        self.forward_readable_aas = list(range(0, self.max_readable_aa_ind_f))
        self.reverse_readable_aas = list(range(self.min_readable_aa_ind_r, self.refseq_len))
                
        # Get the set of indices that is captured by both the forward and reverse reads
        self.double_count_inds = set(self.forward_readable_aas) & set(self.reverse_readable_aas)
        
        # Get the og mutable inds
        self.og_mutable = np.array(
            list(set(self.forward_readable_aas + self.reverse_readable_aas)),
            dtype = int
        )
                
        # Remove the last element from the forward and first from the reverse.
        last_f_el = self.forward_readable_aas.pop()
        first_r_el = self.reverse_readable_aas.pop(0)
        
        # Remove the last and next-to-last elements of f from r. This makes sure
        # that neither forward nor reverse captures the edge cases. We do this
        # to not worry about writing separate code for counting bases.
        if last_f_el in self.reverse_readable_aas:
            self.reverse_readable_aas.remove(last_f_el)
        if first_r_el in self.forward_readable_aas:
            self.forward_readable_aas.remove(first_r_el)
        
        # Get all mutable positions
        mutable_aa_inds = self.forward_readable_aas + self.reverse_readable_aas
                
        # Get only unique mutable inds
        mutable_aa_set = set(mutable_aa_inds)
        self.mutable_aa_inds = np.array(list(mutable_aa_set), dtype = int)
        self.mutable_aa_inds.sort()
        
        # Checks on the mutable inds
        all_possible_inds = set(range(self.refseq_len))
        assert mutable_aa_set.issubset(all_possible_inds)
        
    def assign_qualities(self, min_q_allowed):
        """
        Assigns base quality scores to all sequence components. This includes:
        1. Primer seeds
        2. Frameshifts
        3. The codon refseq
        """
        # Create a quality generator
        q_generator = QualityGenerator(min_q_allowed)
        
        # Primer seed qualities
        self.primer_qualities_f = q_generator.generate_qualities(self.primer_seed_len_f)
        self.primer_qualities_r = q_generator.generate_qualities(self.primer_seed_len_r)
        
        # Barcode qualities
        self.fbc_qualities = q_generator.generate_qualities(BARCODE_LENGTH)
        self.rbc_qualities = q_generator.generate_qualities(BARCODE_LENGTH)
        
        # Adapter qualities
        self.adapter_qualities_f = q_generator.generate_qualities(ADAPTER_LENGTH_F)
        self.adapter_qualities_r = q_generator.generate_qualities(ADAPTER_LENGTH_R)
        
        # Frameshift qualities
        self.frameshift_front_qualities = q_generator.generate_qualities(self.frameshift_front)
        self.frameshift_back_qualities = q_generator.generate_qualities(self.frameshift_back)
        
        # Variable region base qualities
        self.base_variable_qualities = q_generator.generate_qualities(self.codon_refseq_len)
        
    @staticmethod
    def aa_seq_to_codon_seq(aa_seq):
        return [test_glob.RANDOM_RNG.choice(CODON_TABLE.aa_to_codon[aa]) 
                for i, aa in enumerate(aa_seq)]         
        
    @property
    def refseq_bp_seq(self):
        """
        Returns the sequence of the refseq that will be fed into evSeq
        """
        return "".join(
            self.primer_seed_f +
            self.frameshift_bp_front +
            self.codon_refseq +
            self.frameshift_bp_back +
            self.primer_seed_r
        )
        
    @property
    def f_primer(self):
        """
        Returns the forward primer for the reference sequence.
        """
        return ADAPTER_F + "".join(self.primer_seed_f)
    
    @property
    def r_primer(self):
        """
        Returns the reverse primer for the reference sequence.
        """
        return ADAPTER_R + reverse_complement("".join(self.primer_seed_r))
    
    @property
    def frame_distance(self):
        """
        Where we expect the first in-frame base to be. 
        """
        return self.frameshift_front
    
    @property
    def variable_region(self):
        return "".join(
            self.frameshift_bp_front +
            self.codon_refseq +
            self.frameshift_bp_back
        )