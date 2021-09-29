# Import globals module
import tests.data_generation.globals as test_glob

# Import testing globals that never change
from tests.data_generation.globals import (
    ALLOWED_NUCLEOTIDES, MAX_N_READS, MIN_N_READS, MAX_N_VARIANTS, 
    MIN_N_VARIANTS, MIN_DUD_READS, MAX_DUD_READS, MIN_INDELS_ADDED,
    MAX_INDELS_ADDED, MAX_QUAL_ALLOWED, 
)

# Import evSeq globals needed for testing
from evSeq.util.globals import (
    ADAPTER_LENGTH_F, ADAPTER_LENGTH_R, BARCODE_LENGTH, ADAPTER_F, ADAPTER_R
)

# Import custom objects
from .data_generation_utils import reverse_complement, ord_to_chr
from .variant_generator import FakeVariant

# Import 3rd party modules
import numpy as np
import warnings

# Class that holds information for a test well
class FakeWell():
    def __init__(self, config, reference_sequence):
    
        # Assign the config and reference sequence objects as instance variables
        self.config = config
        self.refseq = reference_sequence
        
        # Create plate, well, and barcode variables as a placeholder.
        # This will be filled when the well is passed with others to
        # a FakeRun instance.
        self.platename = None
        self.wellname = None
        self.f_barcode = None
        self.r_barcode = None
                
        # Get the total number of reads in the well. 
        self.total_reads = test_glob.NP_RNG.integers(MIN_N_READS, MAX_N_READS)
        
        # Decide on how many variants we want in the well. 
        self.assign_n_variants()
        
        # Decide the relative abundance of each variant, then build variants
        abundances = self.calculate_variant_abundances()
        
        # Only continue if there are any variants
        if abundances is None:
            self.dud_well = True
        else:
            self.dud_well = False
            variant_abundances, minimum_reads_per_variant = abundances
            self.variants = [FakeVariant(self, abundance, minimum_reads_per_variant) for
                             abundance in variant_abundances]
         
    def assign_n_variants(self):
        """
        Determines how many variants are in the well.
        """
        # Decide on how many variants we want in the well. The maximum allowed is
        # the minimum of (1) how many variants we can spread reads over to get above
        # the `variable_count` threshold, (2) the most allowed with the given 
        # `variable_thresh`, and (3) the `MAX_N_VARIANTS` global.
        count_divisor = 1 if self.config.variable_count == 0 else self.config.variable_count
        max_n_variants = min(
            self.total_reads // count_divisor,
            int(1 // self.config.variable_thresh) - 1,
            MAX_N_VARIANTS
        )
        
        # Make it so that we always have at least one variant
        if max_n_variants == 0:
            max_n_variants = 1
        self.n_variants = test_glob.NP_RNG.integers(MIN_N_VARIANTS, max_n_variants + 1)
    
    def calculate_variant_abundances(self):
        """
        Decide how many counts go to each variant.
        """
        # Each variant must be more abundant than both the variable threshold
        # and the variable counts
        minimum_reads_per_variant = max(
            int(np.ceil(self.config.variable_thresh * self.total_reads)),
            self.config.variable_count
        )
        
        # If we cannot get reads for all variants, this is a dud well
        if minimum_reads_per_variant * self.n_variants > self.total_reads:
            return None
        
        # Now assign read counts to each variant. If there is only 1 variant, then it gets
        # all reads
        if self.n_variants == 1:
            variant_counts = [self.total_reads]

        # If there are more than 1 variants, for each variant, we sample
        # from the range of minimum reads per variant to total reads 
        # remaining, considering that some reads have already been
        # assigned to variants
        else:
            total_reads_available = self.total_reads
            variants_remaining = self.n_variants
            variant_counts = [None] * self.n_variants
            for i in range(self.n_variants):

                # Get the maximum number of reads that we can sample
                max_reads_available_ind = FakeWell.calculate_maximum_reads_ind(total_reads_available, 
                                                                               variants_remaining,
                                                                               minimum_reads_per_variant)

                # Sample to get the number of variants
                if minimum_reads_per_variant == max_reads_available_ind:
                    sampled_variants = minimum_reads_per_variant
                else:
                    sampled_variants = test_glob.NP_RNG.integers(minimum_reads_per_variant, max_reads_available_ind)
                variant_counts[i] = sampled_variants

                # Update the number of reads available for sampling and the number of variants
                # still needing samples
                total_reads_available -= sampled_variants
                assert total_reads_available > 0
                variants_remaining -= 1

                # If this is the breakpoint, assign the remaining reads to the remaining variant
                if variants_remaining == 1:
                    variant_counts[-1] = total_reads_available
                    break

        # Checks to make sure the counts were assigned appropriately
        assert not any(count is None for count in variant_counts)
        assert sum(variant_counts) == self.total_reads   
        
        return variant_counts, minimum_reads_per_variant
    
    def build_corrupted_reads(self):
        """
        Assigns `corrupted_reads` and `corrupted_bases` variables to the instance.
        """
        # Choose how many sequences to add of each flavor.
        corruption_levels = test_glob.NP_RNG.integers(MIN_DUD_READS,
                                            MAX_DUD_READS,
                                            size = 3)
        
        # Add reads with indels. 
        indel_reads = self.build_indel_reads(corruption_levels[0])
        
        # Add sequences filtered out by the average_q_cutoff. 
        low_q_reads = self.build_low_q_reads(corruption_levels[1])
        
        # Add sequences filtered out by the length_cutoff. This is added to `corrupted_reads`.
        short_reads = self.build_short_reads(corruption_levels[2])
        
        return indel_reads, low_q_reads, short_reads
        
    def build_indel_reads(self, n_indels):
        
        # Get the base sequences and qualities
        f_indel_reads = [list(self.base_refseq) for _ in range(n_indels)]
        f_indel_qs = [self.refseq.base_variable_qualities.copy() for _ in range(n_indels)]
        r_indel_reads = [list(self.base_refseq) for _ in range(n_indels)]
        r_indel_qs = [self.refseq.base_variable_qualities.copy() for _ in range(n_indels)]

        # Create lists for storing fastq entries
        r1s = [None] * n_indels
        r2s = [None] * n_indels
        target_f_read = [True, False]

        # Add indels. 
        for read_target in target_f_read:

            # Set targets
            if read_target:
                mutable_positions = self.refseq.forward_readable_aas
                indel_reads = f_indel_reads
                indel_qs = f_indel_qs
            else:
                mutable_positions = self.refseq.reverse_readable_aas
                indel_reads = r_indel_reads
                indel_qs = r_indel_qs

            # Add indels
            for i in range(n_indels):

                # Decide if adding an insertions or deletions and the number
                # to add
                insertion = True if 0.5 < test_glob.NP_RNG.uniform() else False
                n_muts = test_glob.NP_RNG.integers(MIN_INDELS_ADDED, MAX_INDELS_ADDED)

                # Decide on locations of indels
                indel_aa_locs = test_glob.NP_RNG.choice(mutable_positions,
                                              size = n_muts,
                                              replace = False)
                indel_codon_locs = test_glob.NP_RNG.integers(0, 3, size = n_muts)
                indel_locs = indel_aa_locs * 3 + indel_codon_locs
                indel_locs = np.sort(indel_locs)[::-1] # So that largest indices are removed first
                
                # Add indels to sequences
                for indel_loc in indel_locs:

                    # If insertion, add a random character
                    if insertion:
                        indel_reads[i].insert(indel_loc, test_glob.RANDOM_RNG.choice(ALLOWED_NUCLEOTIDES))
                    else:
                        del(indel_reads[i][indel_loc])

                # Add indels to qualities
                if insertion:
                    new_qs = test_glob.NP_RNG.integers(self.config.bp_q_cutoff, MAX_QUAL_ALLOWED, size = n_muts)
                    indel_qs[i] = np.concatenate((indel_qs[i], new_qs))
                else:
                    indel_qs[i] = indel_qs[i][:-n_muts]

                # Convert reads to strings
                indel_reads[i] = "".join(indel_reads[i])

        # Store results
        for i in range(n_indels):
            r1s[i], r2s[i] = self.build_fastq_entry(f_indel_reads[i],
                                                    r_indel_reads[i],
                                                    f_indel_qs[i],
                                                    r_indel_qs[i],
                                                    f"{self.platename}_{self.wellname}_indel_{i}")
        
        return r1s, r2s        
        
    def build_low_q_reads(self, n_low_q):
        
        # Create as many copies of the reference sequence as we want bad sequences
        bad_seq_copies = [self.base_refseq] * n_low_q
        bad_seq_quals = test_glob.NP_RNG.integers(0, self.config.average_q_cutoff,
                                        size = (n_low_q, len(self.base_refseq)))

        # Build the fastq entries
        forward_bad_q = [None] * n_low_q
        reverse_bad_q = [None] * n_low_q
        for i, (bad_seq_copy, bad_seq_qual) in enumerate(zip(bad_seq_copies, bad_seq_quals)):
            forward_bad_q[i], reverse_bad_q[i] = self.build_fastq_entry(bad_seq_copy,
                                                                        bad_seq_copy,
                                                                        bad_seq_qual,
                                                                        bad_seq_qual,
                                                                        f"{self.platename}_{self.wellname}_lowq_{i}")
            
        return forward_bad_q, reverse_bad_q
    
    def build_short_reads(self, n_short):

        # Identify the longest possible short read
        max_allowed_length = int(self.config.length_cutoff * self.config.readlength)
        
        # Determine the maximum allowed readable window that gives the longest
        # possible short read
        forward_readable_max = (max_allowed_length - 
                                self.refseq.frameshift_front - 
                                self.refseq.primer_seed_len_f -
                                ADAPTER_LENGTH_F - BARCODE_LENGTH)
        reverse_readable_max = (max_allowed_length -
                                self.refseq.frameshift_back -
                                self.refseq.primer_seed_len_r -
                                ADAPTER_LENGTH_R - BARCODE_LENGTH)
        
        # If we cannot add short reads, just return empty lists
        if forward_readable_max <= 0 or reverse_readable_max <= 0:
            return ([], [])
        
        # Sample the lengths of the bad fragments
        bad_fragment_lengths_f = test_glob.NP_RNG.integers(0, forward_readable_max, 
                                                 size = n_short)
        bad_fragment_lengths_r = test_glob.NP_RNG.integers(0, reverse_readable_max, 
                                                 size = n_short)

        # Build the fragments
        r1 = [None] * n_short
        r2 = [None] * n_short
        for i, (frag_len_f, frag_len_r) in enumerate(zip(bad_fragment_lengths_f, bad_fragment_lengths_r)):

            # Slice the base reference sequence and the base qualities
            min_r = self.refseq.codon_refseq_len - frag_len_r
            newseq_f = self.base_refseq[:frag_len_f]
            newseq_r = self.base_refseq[min_r:]
        
            newq_f = self.refseq.base_variable_qualities[:frag_len_f]
            newq_r = self.refseq.base_variable_qualities[min_r:]

            # Create entries
            r1[i], r2[i] = self.build_fastq_entry(newseq_f, newseq_r,
                                                  newq_f, newq_r, 
                                                  f"{self.platename}_{self.wellname}_short_{i}")
            
        return r1, r2  
            
    # Build a fastq entry
    def build_fastq_entry(self,
                          full_forward_bp,
                          full_rev_bp,
                          full_forward_q, 
                          full_rev_q, 
                          read_id):
        """
        Sequences and qs should be in the order we expect to see them in the fastq (e.g,
        the reverse complement should be taken of the forward seq before going into
        this function)
        """
        # Get the start sequences and qualities
        start_f_seq, start_f_q = self.build_start_seq_f()
        start_r_seq, start_r_q = self.build_start_seq_r()
        
        # Build the start sequences for the forward and reverse reads.
        # We need the reverse complement of the reverse primer
        available_read_f = self.config.readlength - len(start_f_seq)
        available_read_r = self.config.readlength - len(start_r_seq)
        
        # Get the forward and reverse readable sequences
        forward_readable = full_forward_bp[:available_read_f]
        reverse_readable = reverse_complement(full_rev_bp)[:available_read_r]

        # Get the forward and reverse readable qualities
        forward_qs = ord_to_chr(full_forward_q[:available_read_f])
        reverse_qs = ord_to_chr(np.flip(full_rev_q)[:available_read_r])
                
        # Complete the forward and reverse sequences
        complete_f_seq = start_f_seq + forward_readable
        complete_r_seq = start_r_seq + reverse_readable
        complete_f_qual = start_f_q + forward_qs
        complete_r_qual = start_r_q + reverse_qs

        # Record fastq entries
        r1 = f"@{read_id}\n{complete_f_seq}\n+\n{complete_f_qual}"
        r2 = f"@{read_id}\n{complete_r_seq}\n+\n{complete_r_qual}"
        
        return r1, r2
    
    def build_start_seq_f(self):
        """
        Returns the sequence and qualities for the forward primer, adapter, barcode,
        and frameshift as we would see them in a fastq file
        """
        # Build the sequence
        forward_seed_and_shift = "".join(self.refseq.primer_seed_f + 
                                         self.refseq.frameshift_bp_front)
        forward_read_start = (self.f_barcode + ADAPTER_F + forward_seed_and_shift)

        # Build qualities
        quals = np.concatenate((
            self.refseq.fbc_qualities,
            self.refseq.adapter_qualities_f,
            self.refseq.primer_qualities_f,
            self.refseq.frameshift_front_qualities
        ))
        assert len(quals) == len(forward_read_start)
        
        return forward_read_start, ord_to_chr(quals)
        
    def build_start_seq_r(self):
        """
        Returns the sequence and qualities for the reverse primer, adapter, barcode,
        and frameshift as we would see them in a fastq file
        """
        # Build the sequence
        reverse_seed_and_shift = reverse_complement("".join(self.refseq.frameshift_bp_back +
                                                            self.refseq.primer_seed_r))
        reverse_read_start = (self.r_barcode + ADAPTER_R + reverse_seed_and_shift)
        
        # Build qualities
        quals = np.concatenate((
            self.refseq.rbc_qualities,
            self.refseq.adapter_qualities_r,
            self.refseq.frameshift_back_qualities,
            self.refseq.primer_qualities_r
        ))
        
        assert len(quals) == len(reverse_read_start)
        
        return reverse_read_start, ord_to_chr(quals)
    
    def build_all_reads(self):
        """
        Generates reads for a well.
        """
        # Build lists for holding outputs
        forward_reads = []
        reverse_reads = []

        # Build perfect reads
        for i, variant in enumerate(self.variants):
            f_perfect, r_perfect = variant.build_perfect_reads(i)
            forward_reads.extend(f_perfect)
            reverse_reads.extend(r_perfect)

        # Now augment with bad reads
        for f_bad, r_bad in self.build_corrupted_reads():
            forward_reads.extend(f_bad)
            reverse_reads.extend(r_bad)
        
        return forward_reads, reverse_reads
    
    def build_refseq_entry(self, include_nnn):
        """
        Builds a row in the refseq.csv file for the well. Includes "NNN" for the
        first variant if requested.
        """
        # Get the variable region
        if self.dud_well:
            variable_region = self.refseq.variable_region
        else:
            variable_region = self.variants[0].build_variable_region(include_nnn)
        
        return (
            self.platenickname,
            self.platename,
            self.wellname,
            self.refseq.f_primer,
            self.refseq.r_primer,
            variable_region,
            self.refseq.frame_distance,
            self.refseq.bp_ind_start,
            self.refseq.aa_ind_start
        )
        
    
    def build_output_counts(self):
        """
        Builds output files for the different `OutputCounts`
        """
        pass
        
    
    def murder_well(self):
        """
        With some frequency, screws up the variants in a well to force a DEAD variant.
        """
        pass

    
    @staticmethod
    def calculate_maximum_reads_ind(total_reads_available, total_variants, minimum_reads_per_variant):
        
        # Get the upper bound of sampling maximum reads
        return total_reads_available - ((total_variants - 1) * minimum_reads_per_variant) + 1
    
    @property
    def base_refseq(self):
        return "".join(self.refseq.codon_refseq)
    
    @property
    def platenickname(self):
        return f"TestPlate{self.platename[2:]}"
    