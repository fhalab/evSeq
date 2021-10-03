"""
Contains the functions needed for running and evaluating stress tests.
"""
# Import 3rd party code
import itertools
import numpy as np

# Calculates expected counts for a parent well
def calculate_parent_counts(well):
    
    # Get the complete expected counts across the readable region
    total_counts = np.zeros(well.refseq.refseq_len)
    for variant in well.variants:
        total_counts += variant.expected_aa_counts    
    
    # Get the number of codons added due to the length of the forward and
    # reverse primers/frameshifts
    added_len_f = well.refseq.primer_seed_len_f + well.refseq.frameshift_front
    added_len_r = well.refseq.primer_seed_len_r + well.refseq.frameshift_back
    n_added_codons = added_len_f // 3 + added_len_r // 3    
        
    # Get total counts
    usual_counts = np.sum(total_counts[well.refseq.og_mutable])
    added_codon_counts = n_added_codons * sum(variant.total_counts for variant in well.variants)
    total_counts = usual_counts + added_codon_counts
    
    # Get the mean number of counts
    return int(total_counts / (n_added_codons + len(well.refseq.og_mutable)))

def check_well_is_parent(well, all_mutated_positions, nnn_positions):
    
    # If there are nnn positions, this cannot be a parent seq
    if len(nnn_positions) > 0:
        return False
    
    # If there are no mutated positions, this is a parent well
    if len(all_mutated_positions) == 0:
        return True
        
    # Get the reference sequence for the variants
    parent_seq = well.refseq.aa_refseq

    # Loop over all variants and check to see if all mutated
    # positions match the parent
    parent_checks = [all(variant.base_mut_aa_seq[pos] == parent_seq[pos]
                         for pos in all_mutated_positions)
                    for variant in well.variants]
    
    # If all variants are parent, this is a parent well
    if all(parent_checks):
        return True
    
    # Not a parent if nothing else was triggered
    return False    

def test_decoupled_aa(expected_out, true_out):
        
    # Go row by row and make sure the two dataframe align
    for expected_row, true_row in itertools.zip_longest(expected_out.itertuples(index = False),
                                                        true_out.itertuples(index = False)):

        # Convert to dicts
        expected_row = expected_row._asdict()
        true_row = true_row._asdict()

        # Check if the rows are equal.
        row_passes = True
        for key, true_val in true_row.items():

            # Get the expected value
            expected_val = expected_row[key]

            # Special cases 
            if key == "AlignmentFrequency":
                is_equiv = (np.isclose(true_val, expected_val))
            elif key == "Flags" and (true_val is np.nan):
                is_equiv = (true_val is expected_val)
            elif key == "AaPosition":
                is_equiv = (str(true_val) == str(expected_val))
            else:
                is_equiv = (true_val == expected_val)

            # Break the loop if the true and expected values are not equivalent
            if not is_equiv:
                row_passes = False
                break

        # Break the loop if the rows are not equivalent
        if not row_passes:
            print("FAILURE FOUND")
            return False, expected_out, true_out, expected_row, true_row

    # If we get to this point and row passes, we are successful
    if row_passes:
        print("SUCCESS!!!!")
        return True, expected_out, true_out, expected_row, true_row