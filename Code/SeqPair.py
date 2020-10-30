# Import deSeq dependencies
from .Globals import (BARCODE_LENGTH, ADAPTER_LENGTH_F, ADAPTER_LENGTH_R,
                      BP_TO_IND, AA_TO_IND, CODON_TABLE)
from .Alignment import deseq_align

# Import other required modules
import numpy as np

# Define an object that holds BioPython SeqRecords
class SeqPair():
    
    # Record that we don't have forward or reverse information yet. Record that
    # we don't have alignment information yet
    def __init__(self):
        
        # Did sequence pass QC?
        self._use_f = False
        self._use_r = False
        
        # Did alignments pass QC?
        self._use_f_alignment = False
        self._use_r_alignment = False
            
    # Assign forward reads
    def assign_f(self, f_record):
        
        # Build summary stats
        self._f_barcode, self._f_len, self._f_average_q = self.calculate_read_stats(f_record)
                
        # Assign the forward barcode and the adapterless sequence
        self._f_adapterless = f_record[(BARCODE_LENGTH + ADAPTER_LENGTH_F):]
        
        # Note that we have forward information
        self._use_f = True
    
    # Assign a paired reverse read
    def assign_r(self, r_record):
        
        # Build summary stats
        self._r_barcode, self._r_len, self._r_average_q = self.calculate_read_stats(r_record)
        
        # Assign the reverse barcode and adapterless sequence. 
        # We want to have the reverse complement of this sequence to match
        # the reference sequence
        self._sliced_r = r_record[(BARCODE_LENGTH + ADAPTER_LENGTH_R):]
        self._r_adapterless = self.sliced_r.reverse_complement(id = True, description = True)
        
        # Note that we have reverse information
        self._use_r = True
        
    # Calculate summary stats for a read
    def calculate_read_stats(self, record):
        
        # Calculate relevant summary stats needed by both forward and reverse methods
        barcode = str(record.seq)[:BARCODE_LENGTH]
        length = len(record)
        average_quality = np.mean(record.letter_annotations["phred_quality"])
        
        return barcode, length, average_quality
    
    # Write a function that performs QC on reads
    def qc_reads(self, length_filter, average_q_cutoff):
        
        # Make sure forward reads pass the length and quality cutoff. Only run QC
        # if we actually recorded a forward read.
        if self.use_f:
            if self.f_len < length_filter or self.f_average_q < average_q_cutoff:
                self._use_f = False
                
        # Make sure reverse reads pass the length and quality cutoff. Only run QC
        # if we actually recorded a forward read.
        if self.use_r:
            if self.r_len < length_filter or self._r_average_q < average_q_cutoff:
                self._use_r = False
    
    # Align forward and reverse reads to a reference sequence
    def align(self, reference):
        
        # Make a pairwise alignment. Only align the reads we are using.
        if self.is_paired():
            self._f_alignment = deseq_align(reference, self.f_adapterless.seq)
            self._r_alignment = deseq_align(reference, self.r_adapterless.seq)
            
        # If we are only using forward read, handle this here
        elif self.use_f:
            self._f_alignment = deseq_align(reference, self.f_adapterless.seq)
            self._r_alignment = None
            
        # If we are only using reverse read, handle this here
        elif self.use_r:
            self._f_alignment = None
            self._r_alignment = deseq_align(reference, self.r_adapterless.seq)
            
        else:
            raise AssertionError("No reads to align in reference.")
        
    # Write a function that runs QC on an alignment. We automatically discard an alignment
    # with an insertion or deletion. 
    def qc_alignment(self, forward_check):
        
        # By default, this is a good alignment. If it fails the qc tests, then it
        # will become a bad alignment
        good_alignment = True
        
        # Pull the appropriate alignment 
        test_alignment = self.f_alignment if forward_check else self.r_alignment
        
        # If the alignment is None, this is a bad alignment (the original sequence
        # failed qc) and we cannot conitnue.
        if test_alignment is None:
            return False, -1
        
        # Pull the reference sequence and alignmed sequence
        refseq = test_alignment.seqA
        aligned_seq = test_alignment.seqB

        # If there are any dashes in the reference sequence, we have a bad alignment
        if "-" in refseq:
            good_alignment = False

        # Get the stripped down aligned sequences
        lstripped = aligned_seq.lstrip("-")
        rstripped = aligned_seq.rstrip("-")
        
        # If this is a forward check, dashes in the middle or to the left of the aligned
        # sequence indicate a deletion or insertion
        if forward_check:

            # Check to see if we have an insertion or deletion
            if len(lstripped) < len(aligned_seq) or "-" in rstripped:
                good_alignment = False

            # Get the first instance of a dash in the full sequence. This indicates the
            # first character after the alignment ends
            first_dash = lstripped.find("-")

            return good_alignment, first_dash

        # If this is a reverse check, dashes in the middle or to the right of the aligned
        # sequence indicate a deletion or insertion
        else:

            # Check to see if we have an insertion or deletion
            if len(rstripped) < len(aligned_seq) or "-" in lstripped:
                good_alignment = False

            # Get the last instance of a dash in the full sequence. This indicates the last
            # character index before the aligned sequence begins.
            last_dash = rstripped.rfind("-")

            return good_alignment, last_dash

    # Write a function that runs QC on a pair of alignments. This will set flags for whether
    # or not an alignment is usable
    def qc_alignments(self):
        
        # Run QC on the forward and reverse alignments
        self._use_f_alignment, self._first_dash = self.qc_alignment(True)
        self._use_r_alignment, self._last_dash = self.qc_alignment(False)
        
    # Build a composite alignment for paired ends
    def build_paired_composite_alignment(self):
        
        # Both forward and reverse reads must pass aligment qc to enable this 
        assert self.is_paired_post_alignment_qc(), "Cannot build composite from 1 read."
        
        # Grab the reference sequence, the aligned sequences, 
        # and the quality scores
        refseq = self.f_alignment.seqA
        reflength = len(refseq)
        forward_seq = self.f_alignment.seqB
        reverse_seq = self.r_alignment.seqB
        forward_qual = np.array(self.f_adapterless.letter_annotations["phred_quality"])
        reverse_qual = np.array(self.r_adapterless.letter_annotations["phred_quality"])

        # Get the end of the f read. If it goes all the way to the end of the reference
        # sequence, then the first non-f character is the length of the sequence
        post_forward_dash_ind = reflength if self.first_dash == -1 else self.first_dash

        # Get the beginning of the r read. If it starts from the beginning of the ference
        # sequence, then the first non-r character is -1, so we don't actually need to
        # make any adjustments
        pre_reverse_dash_ind = self.last_dash
        first_r_char_ind = pre_reverse_dash_ind + 1

        # See if the forward and reverse overlap. If they don't overlap. Then the composite
        # is just the forward DNA + dashes + reverse DNA
        if post_forward_dash_ind <= pre_reverse_dash_ind:

            # Calculate the number of dashes needed
            n_dashes = pre_reverse_dash_ind - post_forward_dash_ind + 1

            # Build the composite sequence between the two
            composite_seq = "".join((forward_seq[:post_forward_dash_ind], 
                                   "-" * n_dashes,
                                   reverse_seq[first_r_char_ind:]))

            # Build the composite quality. The quality scores are not extended for the
            # alignment, and so map directly to the pulled sequences.
            composite_qual = np.concatenate((forward_qual,
                                             np.full(n_dashes, np.inf),
                                             reverse_qual))

        # Otherwise, take the sequence with the highest quality in the overlapping region
        else:

            # Pull the forward up to the start of the reverse sequence
            only_f_seq = forward_seq[:first_r_char_ind]
            only_f_qual = forward_qual[:first_r_char_ind]

            # Pull the reverse after the end point of forward. Quality scores only cover 
            # sequence (not alignment gaps), so we need to calculate where the qualities
            # end for the reverse sequence.
            only_r_seq = reverse_seq[post_forward_dash_ind:]
            reverse_qual_break = len(reverse_qual) - len(only_r_seq)
            only_r_qual = reverse_qual[reverse_qual_break:]

            # Now compare the middle parts. Take the one with the higher sequence quality. 
            # The middle characters all fall 
            middle_f_seq = forward_seq[first_r_char_ind:post_forward_dash_ind]
            middle_f_qual = forward_qual[first_r_char_ind:]
            middle_r_seq = reverse_seq[first_r_char_ind:post_forward_dash_ind]
            middle_r_qual = reverse_qual[:reverse_qual_break]

            # The middle sequences should be equal in length (they might differ in sequence
            # due to sequencing errors.  The quality scores should have the same length as well
            middle_size = len(middle_f_seq)
            assert middle_size == len(middle_r_seq)
            assert middle_size == len(middle_f_qual)
            assert middle_size == len(middle_r_qual)

            # Build the composite middle sequence and quality
            middle_seq = [None] * middle_size
            middle_qual = np.zeros(middle_size, dtype = int)
            quality_comparison = np.greater(middle_f_qual, middle_r_qual).astype(int)
            for i in range(middle_size):

                # If the reverse read has better quality, use that
                if quality_comparison[i]:
                    middle_seq[i] = middle_r_seq[i]
                    middle_qual[i] = middle_r_qual[i]

                # If the forward read has better quality, use that
                else:
                    middle_seq[i] = middle_f_seq[i]
                    middle_qual[i] = middle_f_qual[i]

            # Build the overall composite sequence and qualities. 
            composite_seq = "".join((only_f_seq, "".join(middle_seq), only_r_seq))
            composite_qual = np.concatenate((only_f_qual, middle_qual, only_r_qual))
            
        # Check to be sure lengths are correct
        assert reflength == len(composite_seq)
        assert reflength == len(composite_qual)
            
        return composite_seq, composite_qual
    
    # Build a pairwise composite alignment for non-paired ends
    def build_unpaired_composite_alignment(self):
        
        # First make sure that we are calling this function appropriately
        assert not self.is_paired_post_alignment_qc(), "This function only works for unpaired reads"
        
        # Determine if it is forward or reverse reads
        if self.use_f_alignment:
            
            # Get the length of the reference sequence
            refseq = self.f_alignment.seqA
            composite_length = len(refseq)
            
            # The composite sequence is just the aligned sequence
            composite_seq = self.f_alignment.seqB
            
            # The qualities continue after the alignment. Add as many zeros as 
            # there are differences between existing qualities and the end of
            # the sequence
            forward_qual = self.f_adapterless.letter_annotations["phred_quality"]
            composite_qual = np.concatenate((forward_qual, 
                                             np.full(composite_length - len(forward_qual), np.inf)))
            
        else:
            
            # Get the length of the reference sequence
            refseq = self.r_alignment.seqA
            composite_length = len(refseq)
            
            # The composite sequence is just the aligned sequence
            composite_seq = self.r_alignment.seqB
            
            # The qualities must be before the alignment. Prepend as many zeros
            # as there are differences between existing qualities and the end of 
            # the sequence
            reverse_qual = self.r_adapterless.letter_annotations["phred_quality"]
            composite_qual = np.concatenate((np.full(composite_length - len(reverse_qual), np.inf),
                                             reverse_qual))
            
        # Assert that everything is the expected length
        assert composite_length == len(composite_seq)
        assert composite_length == len(composite_qual)
            
        return composite_seq, composite_qual
    
    # Write a function that builds a composite sequence regardless of alignment type
    def build_composite_alignment(self):
        
        # Complicated composite if this is paired end
        if self.is_paired_post_alignment_qc():
            return self.build_paired_composite_alignment()
        
        # Simple composite if this is not paired end
        else:
            return self.build_unpaired_composite_alignment()
        
    # Write a function for extracting information from the alignment
    def analyze_alignment(self, inframe_ind, ref_len, n_aas, qual_thresh):
        
        # Pull the composite alignment for the sequence
        composite_sequence, composite_qual = self.build_composite_alignment()

        # Create matrices in which to store counts
        bp_counts = np.zeros([6, ref_len], dtype = int)
        aa_counts = np.zeros([23, n_aas], dtype = int)

        # Loop over the composite sequence up to the in-frame part
        base_ind = -1 # Initilaize for the case where inframe_ind is 0
        for base_ind, (bp, qual) in enumerate(zip(composite_sequence[:inframe_ind],
                                                  composite_qual[:inframe_ind])):

            # Only record counts if we meet a quality threshold
            if qual >= qual_thresh:
                bp_counts[BP_TO_IND[bp], base_ind] += 1

        # Initialize variables for holding codon information
        aa_counter = 0
        record_aa = True
        codon = [None] * 3
        codon_counter = 0
        
        # Loop over the remaining sequence that is in frame
        for inframe_counter, (bp, qual) in enumerate(zip(composite_sequence[inframe_ind:],
                                                         composite_qual[inframe_ind:])):

            # Update the base ind (this continues from our previous loop)
            base_ind += 1

            # Only record counts if we meet a quality threshold
            if qual >= qual_thresh:
                bp_counts[BP_TO_IND[bp], base_ind] += 1
                codon[codon_counter] = bp

            # If we don't meet a quality threshold, then throw a flag to
            # not record the aa in this codon
            else:
                record_aa = False

            # Increment the codon counter
            codon_counter += 1
            
            # If this is the third character in a codon reset the codon counter
            # and other codon-related variables
            if (inframe_counter + 1) %3 == 0:

                # If all members of the codon passed quality control record
                if record_aa:
                    
                    # Join the characters
                    joined_codon = "".join(codon)
                    
                    # If this is in a gap, record gap
                    if "-" in joined_codon:
                        aa = "-"
                    
                    # If it isn't in the codon table, record question mark
                    elif joined_codon not in CODON_TABLE:
                        aa = "?"
                    
                    else:
                        aa = CODON_TABLE[joined_codon]
                    
                    # Add to counts
                    aa_counts[AA_TO_IND[aa], aa_counter] += 1

                # Reset all codon related variables and increment the aa counter
                aa_counter += 1
                record_aa = True
                codon = [None] * 3
                codon_counter = 0
            
        # Run a check on the count. A sum across the 0th axis should
        # return all ones and zeros, as we should never count two bases or two
        # amino acids in one position
        bp_test = np.sum(bp_counts, axis = 0)
        aa_test = np.sum(aa_counts, axis = 0)
        assert np.all(np.logical_or(bp_test == 1, bp_test == 0)), "Double counting bases"
        assert np.all(np.logical_or(aa_test == 1, aa_test == 0)), "Double counting amino acids"
            
        # Return the filled out count matrices
        return bp_counts, aa_counts
    
    # Write a function that returns read lengths
    def read_lengths(self):
        if self.is_paired():
            return [self.f_len, self.r_len]
        elif self.use_f:
            return [self.f_len, np.nan]
        elif self.use_r:
            return [np.nan, self.r_len]
        else:
            raise AssertionError("No reads for which to return lengths.")
    
    # Check to see if we are using both sequences
    def is_paired(self):
        if self.use_r and self.use_f:
            return True
        else:
            return False
        
    # Check to see if we have no sequences aligned
    def is_dud(self):
        if not (self.use_r or self.use_f):
            return True
        else:
            return False
        
    # Check to see if both alignments pass QC
    def is_paired_post_alignment_qc(self):
        if self.use_f_alignment and self.use_r_alignment:
            return True
        else:
            return False
        
    # Check to see if we have no alignments that pass
    def is_dud_post_alignment_qc(self):
        if not(self.use_f_alignment or self.use_r_alignment):
            return True
        else:
            return False
        
    # Make all the properties
    @property
    def use_f(self):
        return self._use_f
    
    @property
    def use_r(self):
        return self._use_r
    
    @property
    def use_f_alignment(self):
        return self._use_f_alignment
    
    @property
    def use_r_alignment(self):
        return self._use_r_alignment
        
    @property
    def f_barcode(self):
        return self._f_barcode
    
    @property
    def f_len(self):
        return self._f_len
    
    @property
    def f_average_q(self):
        return self._f_average_q
    
    @property
    def f_adapterless(self):
        return self._f_adapterless
    
    @property
    def r_barcode(self):
        return self._r_barcode
    
    @property
    def r_len(self):
        return self._r_len
    
    @property
    def r_average_q(self):
        return self._r_average_q
    
    @property
    def sliced_r(self):
        return self._sliced_r
    
    @property
    def r_adapterless(self):
        return self._r_adapterless
    
    @property
    def f_alignment(self):
        return self._f_alignment
    
    @property
    def r_alignment(self):
        return self._r_alignment
    
    @property
    def first_dash(self):
        return self._first_dash
    
    @property
    def last_dash(self):
        return self._last_dash