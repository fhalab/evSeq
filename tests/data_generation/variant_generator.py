# Import globals module
import tests.data_generation.globals as test_glob

# Import testing globals that never change
from tests.data_generation.globals import (
    MIN_PERC_MUTATED, MAX_PERC_MUTATED, ALLOWED_AAS, CODON_TABLE, 
    MIN_NOISE_PERC, MAX_NOISE_PERC, RESCUE_FREQ, ALLOWED_NUCLEOTIDES
)

# Import 3rd party modules
import warnings
import numpy as np
from copy import deepcopy

# Create a class that holds all relevant information for a variant
class FakeVariant():
    def __init__(self, well, counts, minimum_reads_allowed):
        
        # The well is an instance variables
        self.well = well
        
        # Set the total number of counts for the variant and the
        # minimum reads that are allowed for it
        self.total_counts = counts
        self.minimum_reads_allowed = minimum_reads_allowed
        
        # Don't move forward if there are no counts
        self.no_counts = (counts == 0)
        
        if not self.no_counts:
            
            # Choose the variable positions and build variants
            self.choose_variable_positions()
            self.build_variants()

            # Assign Q scores. Get expected counts.
            self.incorporate_noisy_positions()
        
    def choose_variable_positions(self):
        """
        Chooses the variable positions for the variant and assigns mutant
        amino acids/codons.
        """
        warnings.warn("You need to make sure that the amino acid returned is different"
                      " from the existing one. Same goes for the nucleotides.")
        # Get the number of mutations to make
        n_mutable_positions = len(self.well.refseq.mutable_aa_inds)
        min_n_muts = int(MIN_PERC_MUTATED * n_mutable_positions)
        max_n_muts = int(MAX_PERC_MUTATED * n_mutable_positions) + 1
        self.n_variable_positions = test_glob.NP_RNG.integers(min_n_muts, max_n_muts)
        
        # Decide on the positions that will be varied within the variant
        # This must be in the readable region. 
        self.mutated_positions = test_glob.NP_RNG.choice(self.well.refseq.mutable_aa_inds, 
                                               size = self.n_variable_positions,
                                               replace = False)
        self.mutated_positions.sort()
        
        # Choose variable amino acids
        self.variable_aas = test_glob.RANDOM_RNG.choices(ALLOWED_AAS, k = self.n_variable_positions)
        
        # Choose a random codon for encoding said amino acids. It cannot be the
        # same as the existing codon.
        self.variable_codons = [None] * self.n_variable_positions
        for i, aa in enumerate(self.variable_aas):
            
            # Get the options for the alternate codons. Again, we cannot
            # reuse the codon that is already present in the refseq
            existing_codon = self.well.refseq.codon_refseq[self.mutated_positions[i]]
            codon_opts = [codon for codon in CODON_TABLE.aa_to_codon[aa]
                          if codon != existing_codon]
            
            # If there are no alternate codons (e.g., if we selected methionine as
            # our mutant aa and refseq was already methionine), then just sample from
            # leucine
            if len(codon_opts) == 0:
                self.variable_aas[i] = "L"
                codon_opts = CODON_TABLE.aa_to_codon["L"]
                
            # Select the new codon
            self.variable_codons[i] = test_glob.RANDOM_RNG.choice(codon_opts)
            
    def build_variants(self):
        """
        Builds sequence variants based on the refseq sequence and chosen replacement
        amino acids/positions
        """
        # Create base mutant sequences based on the refseq
        self.base_mut_aa_seq = self.well.refseq.aa_refseq.copy()
        base_mut_codon_seq = self.well.refseq.codon_refseq.copy()

        # Loop over the mutated positions and add to the base sequence
        assert len(self.mutated_positions) == len(self.variable_aas)
        assert len(self.mutated_positions) == len(self.variable_codons)
        for i, mutated_pos in enumerate(self.mutated_positions):
            self.base_mut_aa_seq[mutated_pos] = self.variable_aas[i]
            base_mut_codon_seq[mutated_pos] = self.variable_codons[i]

        # Create forward and reverse copies of the codon sequences. Make as many copies
        # as there are variants for the refseq
        self.mut_codon_seqs_f = [base_mut_codon_seq.copy() for _ in range(self.total_counts)]
        self.mut_codon_seqs_r = deepcopy(self.mut_codon_seqs_f)
        
    def id_noisy_positions(self):
        """
        Adds noise to a number of mutated positions.
        """
        # Decide how many reads for combos of amino acids will have noise
        # added to them. This number must fall above the minimum counts required.
        buffer_region = self.total_counts - self.minimum_reads_allowed
        n_combos_destroyed = test_glob.NP_RNG.integers(buffer_region) if buffer_region > 0 else 0

        # Get the expected number of both bp and aa combinations
        self.expected_combo_counts = self.total_counts - n_combos_destroyed 

        # Decide which reads will have noise added.
        read_inds = np.arange(self.total_counts)
        noisy_reads = test_glob.NP_RNG.choice(read_inds, 
                                    size = n_combos_destroyed,
                                    replace = False)
        noisy_reads.sort()

        # Determine options for noisy positions.
        if len(self.well.refseq.double_count_inds) != 0:
            forbidden_noisy_aas = {min(self.well.refseq.double_count_inds), 
                                   max(self.well.refseq.double_count_inds)}
        else:
            forbidden_noisy_aas = set()
        noisened_position_options = np.array([pos for pos in self.mutated_positions
                                              if pos not in forbidden_noisy_aas])

        # Determine how many noisy positions per read. 
        n_noisey_per_read = int(len(noisened_position_options) * 
                                test_glob.NP_RNG.uniform(MIN_NOISE_PERC, MAX_NOISE_PERC))

        # Decide which amino acids within each read will have noise added.
        # Do not add noise to codons within 1 of the double overlap region.
        noisy_positions = np.array([test_glob.NP_RNG.choice(noisened_position_options, 
                                                  size = n_noisey_per_read,
                                                  replace = False)
                                    for _ in range(n_combos_destroyed)])

        # Determine which bases are noisy for each noisy amino acid. 
        n_noisy_bases_by_noisy_pos = test_glob.NP_RNG.integers(1, 4, size = (n_combos_destroyed, n_noisey_per_read))
        noisy_bases_by_noisy_pos = [[test_glob.NP_RNG.choice(3,
                                                  size = n_noisy_bases,
                                                  replace = False) 
                                     for n_noisy_bases in noisy_base_array]
                                    for noisy_base_array in n_noisy_bases_by_noisy_pos]
        
        # Checks
        assert len(noisy_reads) == len(noisy_positions)
        
        return noisy_reads, noisy_positions, noisy_bases_by_noisy_pos
                            
    def build_expected_count_arrays(self):

        # Create arrays of counts.
        expected_aa_counts = np.full(self.well.refseq.refseq_len, 
                                     self.total_counts)
        expected_bp_counts = np.full(self.well.refseq.codon_refseq_len, 
                                     self.total_counts)

        # Double positions in the counts where we have overlap 
        for i, mutant_pos in enumerate(self.mutated_positions):
            if mutant_pos in self.well.refseq.double_count_inds:

                # Double aa counts
                expected_aa_counts[i] *= 2

                # Double bp counts
                bp_start_ind = i * 3
                for bp_ind in range(bp_start_ind, bp_start_ind + 3):
                    expected_bp_counts[bp_ind] *= 2

        return expected_aa_counts, expected_bp_counts
        
    def incorporate_noisy_positions(self):
        """
        Assigns counts to different bases and amino acids, then assigns quality scores
        to get them to these counts
        """
        # Identify noisy reads, positions, and nucleotides
        noisy_reads, noisy_positions, noisy_bases_by_noisy_pos = self.id_noisy_positions()

        # Build expected output counts for amino acids and bases
        self.expected_aa_counts, self.expected_bp_counts = self.build_expected_count_arrays()
        
        # Create two quality score arrays. One is for the forward read and the other
        # is for the reverse reads
        self.f_quals = np.tile(self.well.refseq.base_variable_qualities.copy(),
                          (self.total_counts, 1))
        self.r_quals = self.f_quals.copy()
        min_existing_q = self.f_quals.min()
        bad_qual_q = min_existing_q - 2

        # Add noise to positions. Adjust counts and qualities accordingly.
        for noisy_read, noisy_position_array, noisy_bp_array in \
            zip(noisy_reads, noisy_positions, noisy_bases_by_noisy_pos):

            # Loop over all positions and adjust basepair quality as appropriate
            for noisy_pos, noisy_base_set in zip(noisy_position_array, noisy_bp_array):

                # Determine count adjustment for the position
                double_count_pos = (noisy_pos in self.well.refseq.double_count_inds)
                count_adj = 2 if double_count_pos else 1

                # If a position is in the double-read region, with some probability
                # one read will rescue the other. If we rescue, then we set only 1 
                # codon/amino acid as low quality. The other one is kept fine. We also
                # mutate the low-quality codon again (this should never be counted, providing
                # a test to make sure that we are appropriately ignoring codons)
                rescue_check = (double_count_pos and (test_glob.NP_RNG.uniform() < RESCUE_FREQ))
                rescue = True if rescue_check else False

                # Get the base index for the noisy position
                noisy_base_index_zero = noisy_pos * 3

                # Loop over the base set and add noise
                for noisy_base in noisy_base_set:

                    # Calculate the actual index
                    actual_base_ind = noisy_base_index_zero + noisy_base

                    # If rescuing, choose one of the quality score arrays to update. Also add 
                    # an error to the low-quality sequence. This should be ignored the evSeq software.
                    if rescue:

                        # Define options for qualities and sequences
                        qual_opts = [self.f_quals, self.r_quals]
                        codon_opts = [self.mut_codon_seqs_f, self.mut_codon_seqs_r]

                        # Choose which quality and sequence array will be update
                        target_ind = 0 if test_glob.NP_RNG.uniform() < 0.5 else 1
                        target_qual_array = qual_opts[target_ind]
                        target_codon_list = codon_opts[target_ind]

                        # Update the qualities in the chosen array
                        target_qual_array[noisy_read, actual_base_ind] = bad_qual_q

                        # Choose a replacement base
                        replacement_codon = list(target_codon_list[noisy_read][noisy_pos])
                        replacement_codon[noisy_base] = test_glob.RANDOM_RNG.choice(ALLOWED_NUCLEOTIDES)
                        target_codon_list[noisy_read][noisy_pos] = "".join(replacement_codon)

                        # Set the count adjuster. It is only 1 here, because we lost a base.
                        count_adj = 1

                    # Otherwise, we adjust both quality arrays and make no changes to sequence
                    else:

                        # Set quality to be 2 less than the minimum existing quality in
                        # the array of q-scores.
                        self.f_quals[noisy_read, actual_base_ind] = bad_qual_q
                        self.r_quals[noisy_read, actual_base_ind] = bad_qual_q

                    # Adjust the counts. 
                    self.expected_aa_counts[noisy_pos] -= count_adj
                    self.expected_bp_counts[actual_base_ind] -= count_adj
        
    def build_perfect_reads(self):
        """
        Assigns a `perfect_reads` variable to the instance.
        """
        # Confirm that we have as many qualities as mutant sequences
        assert len(self.f_quals) == self.total_counts
        assert len(self.f_quals) == len(self.r_quals)
        assert len(self.f_quals) == len(self.mut_codon_seqs_f)
        assert len(self.f_quals) == len(self.mut_codon_seqs_r)

        # Now loop over the remaining sequence variants and build entries.
        fastq_r1 = [None] * self.total_counts
        fastq_r2 = fastq_r1.copy()
        for i in range(self.total_counts):

            # Confirm that mutant sequence length matches with quality score length
            full_forward_bp = "".join(self.mut_codon_seqs_f[i])
            full_forward_q = self.f_quals[i]
            assert len(full_forward_bp) == len(full_forward_q)

            full_rev_bp = "".join(self.mut_codon_seqs_r[i])
            full_rev_q = self.r_quals[i]
            assert len(full_rev_bp) == len(full_rev_q)
            
            assert len(full_forward_bp) == len(full_rev_bp)    
            
            # Record fastq entries
            fastq_r1[i], fastq_r2[i] = self.well.build_fastq_entry(full_forward_bp,
                                                                   full_rev_bp,
                                                                   full_forward_q,
                                                                   full_rev_q,
                                                                   i)

        return fastq_r1, fastq_r2        
    
    def build_output_counts(self):
        """
        Builds output files for the different `OutputCounts`
        """
        pass  
