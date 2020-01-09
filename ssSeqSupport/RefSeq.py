# Import third party modules
import numpy as np

# Import ssSeqSupport objects
from . import Translate, FindNNN
from . import LogError
from . import AdapterLengthF, AdapterLengthR

# Define a class that represents a reference sequence
class RefSeq():

    # Define the initialization. This takes the reference sequence, identifies
    # variable sites, and translates where possible.
    def __init__(self, reference_sequence, variable_positions, n_var_sites):

        # Assign the reference sequence as a class attribute as well as its length
        self._seq = reference_sequence
        self._seq_length = len(reference_sequence)

        # Identify the translation start site
        self._trans_start = variable_positions[0] % 3

        # Translate the reference sequence
        translated_seq = Translate(reference_sequence, self.trans_start)

        # Record the length of the translation
        self._trans_length = len(translated_seq)

        # Record where we can find the amino acids of interest
        self._var_aa_sites = [int(np.floor(site / 3)) for site in variable_positions]
        
        # Write the bp and aa sequences as lists
        self._aas_as_list = [char for char in translated_seq]
        self._bps_as_list = [char for char in self._seq]
        
    # Set properties
    @property
    def seq(self):
        return self._seq
    
    @property
    def seq_length(self):
        return self._seq_length
    
    @property
    def trans_start(self):
        return self._trans_start
    
    @property
    def trans_length(self):
        return self._trans_length
    
    @property
    def var_aa_sites(self):
        return self._var_aa_sites
    
    @property
    def aas_as_list(self):
        return self._aas_as_list
    
    @property
    def bps_as_list(self):
        return self._bps_as_list
    
# Write a function that builds individual reference sequences from the larger
# reference sequence passed in
def BuildRefSeqs(reference_sequence, read_length):
    
    # Make sure the reference sequence is entirely uppercase
    reference_sequence = reference_sequence.upper()
    
    # Identify the variable positions, the total number of variable positions,
    # and whether or not we have found our "N" values in codon format for the 
    # input reference sequence
    _, full_n_sites, full_codon_check = FindNNN(reference_sequence)
    
    # If we did not find this in codon format, throw an error
    if not full_codon_check:
        LogError("""The variable sites must be specified as blocks of 'NNN'.
                 The current reference sequence has 'N's not in blocks of 3.""")
    
    # Construct the raw forward and reverse reference sequences. 
    readable_ref_length_f = read_length - AdapterLengthF
    readable_ref_length_r = read_length - AdapterLengthR
    raw_f = reference_sequence[:readable_ref_length_f]
    raw_r = reference_sequence[-readable_ref_length_r:]
    
    # Make sure that our forward and reverse reference sequences obey the codon
    # rules set. Also calcuate the number of variable positions in them.
    f_var_sites, f_n_sites, _ = FindNNN(raw_f)
    r_var_sites, r_n_sites, _ = FindNNN(raw_r)
    
    # Make sure that we are capturing at least as many variable sites between
    # the f and r raw sequences as we do in the input reference sequence
    if f_n_sites + r_n_sites < full_n_sites:
        LogError("""Your read length does not capture all variable positions in the reference sequence provided.
                 Your forward read captures {} variable positions and the reverse {} variable positions, while
                 your input sequence has {} variable positions.""".format(f_n_sites, r_n_sites, full_n_sites))
        
    # Build a reference sequence object for both the forward and reverse sequences
    f_ref_seq = RefSeq(raw_f, f_var_sites, f_n_sites)
    r_ref_seq = RefSeq(raw_r, r_var_sites, r_n_sites)
    
    # Return both the forward and reverse reference sequence as well as the total
    # number of variable sites
    return f_ref_seq, r_ref_seq, full_n_sites
    