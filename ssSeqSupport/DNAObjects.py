# Import third party packages
import numpy as np

# Import local packages
from . import CodonTable, BarcodeLength
from . import GetBlockInfo, ReverseComplement, Translate

class Translation():

    # Initialize
    def __init__(self, seq, start_ind, low_quality_chars = []):

        # Get the number of codons in the sequence
        n_codons = np.floor((len(seq) - start_ind)/3)

        # Identify all codons
        codons = [seq[int(start_ind + i*3): int(start_ind + (i+1)*3)] for i in range(int(n_codons))]

        # Link bp to codon
        bp_ranges = [list(range(int(start_ind + i*3), int(start_ind + (i+1)*3)))
                     for i in range(int(n_codons))]
        bp_to_codon = {bp: codon for codon, bp_list in enumerate(bp_ranges)
                       for bp in bp_list}

        # Identify low quality codons
        self.low_quality_codons = {bp_to_codon[bp] for bp in low_quality_chars
                                   if bp in bp_to_codon.keys()}

        # Translate all codons
        translation = []
        for codon in codons:

            # If this is "NNN" we give a question mark
            if "N" in codon:
                translation.append("?")
            else:
                translation.append(CodonTable[codon])


        # Store the translation
        self.translation  = "".join(translation)

class SeqPair():

    # Define the initialization. This takes the sequence information and
    # organizes it
    def __init__(self, info_block, start_f = True):

        # Parse the first-seen block to build the object. Set the forward
        # or reverse sequence as appropriate
        if start_f:
            raw_id, self.f_seq, _, f_ASCII_qual_scores = info_block
            self.f_qual_scores = [ord(char)-33 for char in f_ASCII_qual_scores]

        else:
            raw_id, self.r_seq, _, r_ASCII_qual_scores = info_block
            self.r_qual_scores = [ord(char)-33 for char in r_ASCII_qual_scores]

        # Get the id information for the block
        (instrument_name, _, _, lane, tile, x_coord, y_coord,
         _, _, _, sample_number) = GetBlockInfo(raw_id)

        # Create a unique id for the pair that can pair ends
        self.id = (instrument_name, lane, tile, x_coord,
                   y_coord, sample_number)

    # Write a function to attach the partner information
    def attach_partner(self, info_block, r_partner = True):

        # Depending on the partner that's being appended, attach different values
        if r_partner:
            _, self.r_seq, _, r_ASCII_qual_scores = info_block
            self.r_qual_scores = [ord(char)-33 for char in r_ASCII_qual_scores]
        else:
            _, self.f_seq, _, f_ASCII_qual_scores = info_block
            self.f_qual_scores = [ord(char)-33 for char in f_ASCII_qual_scores]

        # Get the forward and reverse lengths
        self.f_len = len(self.f_seq)
        self.r_len = len(self.r_seq)

        # Force everything to be uppercase
        self.f_seq = self.f_seq.upper()
        self.r_seq = self.r_seq.upper()

        # Identify the forward barcode and its prob scores
        self.f_barcode = self.f_seq[:BarcodeLength]

        # Identify the reverse barcode and its prob scores
        self.r_barcode = self.r_seq[:BarcodeLength]

        # Define sequences without the barcodes
        self.barcodeless_f = self.f_seq[BarcodeLength:]
        self.barcodeless_f_q = self.f_qual_scores[BarcodeLength:]
        self.reversed_barcodeless_r = ReverseComplement(self.r_seq[BarcodeLength:])
        self.reversed_barcodeless_r_q = list(reversed(self.r_qual_scores[BarcodeLength:]))

    # Write a function to determine if this is a sequence without a partner
    def is_orphan(self):

        # Test to see if we have forward and reverse sequences
        if hasattr(self, "f_seq") and hasattr(self, "r_seq"):
            return False
        else:
            return True