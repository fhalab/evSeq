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
    # ref_ind_start used with ssSeq_tile 
    def __init__(self, raw_ref_seq, variable_positions, n_var_sites, 
                 args, full_ref_seq, ref_ind_start = 0, forward = True):
        """
        

        Parameters
        ----------
        raw_ref_seq : string
            raw_ref_seq is the forward or reverse substring of the full 
            reference sequence, of length read_length - adapter. For forward
            reference sequence it is that length from start of full_ref_seq and 
            for reverse it is that length from end of full_ref_seq
        variable_positions : TYPE
            DESCRIPTION.
        n_var_sites : TYPE
            DESCRIPTION.
        args : TYPE
            DESCRIPTION.
        full_ref_seq : TYPE
            The full reference sequence given. Used here to determine the 
            reverse ref_seq reading frame.
        ref_ind_start : TYPE, optional
            DESCRIPTION. The default is 0.
        forward : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """

        # Assign the reference sequence as a class attribute as well as its length
        self._seq = raw_ref_seq
        self._seq_length = len(raw_ref_seq)

        # Identify the translation start site
        # For tiles, based on input reference index start
        if args["ssSeq_tile"]:
            if forward:
                # set reading frame based on ref_ind_start frame
                self._trans_start = ref_ind_start % 3
                
                # record the bp index on which reference index starts
                # for use in determining mutation positions.
                self._ref_ind_start = ref_ind_start
                
            if not forward:
                # set reading frame based on ref_ind_start frame and seq len
                # second modulo to get correct reading frame from difference
                self._trans_start = (3 - 
                    (ref_ind_start + len(full_ref_seq) - self._seq_length) % 3) % 3
                
                # record the bp index on which reference index starts
                # for use in determining mutation positions.
                self._ref_ind_start = (ref_ind_start + len(full_ref_seq) 
                                       - self._seq_length)
        else:
            self._trans_start = variable_positions[0] % 3

        # Translate the reference sequence
        translated_seq = Translate(raw_ref_seq, self.trans_start)

        # Record the length of the translation
        self._trans_length = len(translated_seq)

        if args["ssSeq_tile"]:
            self._var_aa_sites = [0,]
        else: # Record where we can find the amino acids of interest
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
    
    @property
    def ref_ind_start(self):
        return self._ref_ind_start
    
# Write a function that builds individual reference sequences from the larger
# reference sequence passed in
def BuildRefSeqs(reference_sequence, read_length, ref_ind_start, args):
    """
    In RefSeq.py
    
    Builds individual reference sequences from larger reference sequence 
    passed in.
    
    AMK updated for tiles

    Parameters
    ----------
    reference_sequence : TYPE
        DESCRIPTION.
    read_length : TYPE
        DESCRIPTION.
    ref_ind_start: where the reference index starts from reference sequence
        csv; used to determine translation start point.

    Returns
    -------
    f_ref_seq : TYPE
        DESCRIPTION.
    r_ref_seq : TYPE
        DESCRIPTION.
    full_n_sites : TYPE
        DESCRIPTION.

    """

    
    # Make sure the reference sequence is entirely uppercase
    reference_sequence = reference_sequence.upper()
    
    if not args["ssSeq_tile"]:
        # Identify the variable positions, the total number of variable positions,
        # and whether or not we have found our "N" values in codon format for the 
        # input reference sequence
        _, full_n_sites, full_codon_check = FindNNN(reference_sequence)
        
        # If we did not find this in codon format, throw an error
        if not full_codon_check:
            LogError("""The variable sites must be specified as blocks of 'NNN'.
                     The current reference sequence has 'N's not in blocks of 3.""")
    else:
        full_n_sites = 0
    
    # Construct the raw forward and reverse reference sequences. 
    readable_ref_length_f = read_length - AdapterLengthF
    readable_ref_length_r = read_length - AdapterLengthR
    raw_f = reference_sequence[:readable_ref_length_f]
    raw_r = reference_sequence[-readable_ref_length_r:]
    
    # No NNNs are given for tile sequencing; skip checking for NNN
    # Build RefSeq with 0 variable sites, feed in reference index start
    # to determine reading frame
    if args["ssSeq_tile"]:
        f_ref_seq = RefSeq(raw_f, 0, 0, args, 
                           reference_sequence, ref_ind_start)
        
        r_ref_seq = RefSeq(raw_r, 0, 0, args, 
                           reference_sequence, ref_ind_start, forward=False)
     
        
    # Make sure that our forward and reverse reference sequences obey the codon
    # rules set. Also calcuate the number of variable positions in them.  
    else:
        f_var_sites, f_n_sites, _ = FindNNN(raw_f)
        r_var_sites, r_n_sites, _ = FindNNN(raw_r)
    
        # Make sure that we are capturing at least as many variable sites between
        # the f and r raw sequences as we do in the input reference sequence
        if f_n_sites + r_n_sites < full_n_sites:
            LogError("""Your read length does not capture all variable positions in the reference sequence provided.
                     Your forward read captures {} variable positions and the reverse {} variable positions, while
                     your input sequence has {} variable positions.""".format(f_n_sites, r_n_sites, full_n_sites))
        
        # Build a reference sequence object for both the forward and reverse sequences
        f_ref_seq = RefSeq(raw_f, f_var_sites, f_n_sites, 
                           args, reference_sequence)
        
        r_ref_seq = RefSeq(raw_r, r_var_sites, r_n_sites, 
                           args, reference_sequence, forward=False)
    
    # Return both the forward and reverse reference sequence as well as the total
    # number of variable sites
    return f_ref_seq, r_ref_seq, full_n_sites
    