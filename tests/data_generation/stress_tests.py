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

def run_aa_stress_test(expected_out, true_out):
    
    # Get the unique plates and wells between the two. 
    expected_platewell = {tuple(platewell) for platewell in 
                          expected_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    true_platewell = {tuple(platewell) for platewell in 
                      true_out.loc[:, ["IndexPlate", "Well"]].values.tolist()}
    
    # Find the differences between the true and expected. These are our first
    # bad plate-well combos.
    bad_platewells = list(expected_platewell ^ true_platewell)
    good_platewells = list(expected_platewell & true_platewell)
    
    # Create a report for the bad platewells
    error_reports = ["Missing well" for _ in range(len(bad_platewells))]
    
    # For all potentially good platewells, test to be sure that we see the same
    # rows coming out of each
    for plate, well in good_platewells:
        
        # Get the limited dataframe
        limited_expected_out = expected_out.loc[(expected_out.IndexPlate == plate)&
                                                (expected_out.Well == well)]
        limited_true_out = true_out.loc[(true_out.IndexPlate == plate)&
                                        (true_out.Well == well)]
        
        # Assert the dataframes are the same size
        len_expected_out = len(limited_expected_out)
        len_true_out = len(limited_true_out)
        assert (len_expected_out > 0) and (len_true_out > 0)
        if len_expected_out != len_true_out:
            error_reports.append(f"Mismatched number of entries for {plate}-{well}")
            bad_platewells.append((plate, well))
            continue
            
        # Go row by row and make sure the two dataframes align
        for expected_row, true_row in itertools.zip_longest(limited_expected_out.itertuples(index = False),
                                                            limited_true_out.itertuples(index = False)):

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

            # Report if the rows are not equivalent
            if not row_passes:
                bad_platewells.append((plate, well))
                error_reports.append([limited_expected_out.copy(),
                                      limited_true_out.copy()])
                break
                
    # Evaluate how many bad platewells were found.
    successful_test = (len(bad_platewells) == 0)
    return successful_test, bad_platewells, error_reports