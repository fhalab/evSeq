"""
Tests for refseq file parsing and sequence construction.
"""
import pandas as pd

from evSeq.util.globals import ADAPTER_F, ADAPTER_R
from evSeq.util.input_processing import construct_ref_seq

def test_default_construct_ref_seq():

    f_primer_seed = 'AAA'
    r_primer_seed = 'TTT'
    variable_region = 'GGGGGGGGGGGGG'
    f_dist = 0
    bp_idx = 0
    aa_idx = 0

    default_refseq_df = pd.DataFrame({
        'PlateName': 'TEST',
        'IndexPlate': 'DI01',
        'FPrimer': ADAPTER_F + f_primer_seed,
        'RPrimer': ADAPTER_R + r_primer_seed,
        'VariableRegion': variable_region,
        'FrameDistance': f_dist,
        'BpIndStart': bp_idx,
        'AaIndStart': aa_idx,
    }, index=[0])

    refseqs, new_frame, new_bp, new_aa = construct_ref_seq(default_refseq_df)

    def rev_comp(seq):
        comp_dict = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
        }
        return ''.join([comp_dict[base] for base in seq[::-1]])
    
    assert refseqs == [f_primer_seed+variable_region+rev_comp(r_primer_seed)]
    # Seed region = 3 bp, so f_dist does not change
    assert new_frame == [f_dist]
    assert new_bp == [idx - len(f_primer_seed) for idx in [bp_idx]]
    # f_dist is zero, so no adjustment needed here
    assert new_aa == [idx - len(f_primer_seed) // 3 for idx in [aa_idx]]


def test_default_construct_ref_seq_with_offsets():

    # If you change these, you will need to change some things below.
    # These are the simplest combos that test all possible change types,
    # So they should not need to be changed.
    offset_f_primer_seeds = ['AAA', 'AAAA', 'AAAAA']
    offset_f_dists = [0, 1, 2]
    
    # Pair together in known way, then unpack
    combos = [(seeds, dists)
              for seeds in offset_f_primer_seeds
              for dists in offset_f_dists]
    f_primer_seed = [combo[0] for combo in combos]
    f_dist = [combo[1] for combo in combos]

    # Adjust length of extra components
    r_primer_seed = ['TTT']*len(f_dist)
    variable_region = ['GGGGGGGGGGGGG']*len(f_dist)

    # These are appropriate bp_idx and aa_idx values for the given
    # combos of f_seeds and frame_distances.
    # each line is starting (bp, aa), # -> expected (bp, aa)
    bp_aa_pairs = [
        # f_seed = 'AAA', f_dists = [0, 1, 2]
        (3, 1), # -> (0, 0)
        (2, 1), # -> (-1, 0)
        (1, 1), # -> (-2, 0)
        # f_seed = 'AAAA', f_dists = [0, 1, 2]
        (3, 1), # -> (-1, 0)
        (2, 1), # -> (-2, 0)
        (1, 1), # -> (-3, -1)
        # f_seed = 'AAAAA', f_dists = [0, 1, 2]
        (3, 1), # -> (-2, 0)
        (3, 1), # -> (-3, -1)
        (3, 1), # -> (-4, -1)
    ]

    # Unpack
    bp_idx = [idx[0] for idx in bp_aa_pairs]
    aa_idx = [idx[1] for idx in bp_aa_pairs]

    # Create refseq_df
    default_refseq_df = pd.DataFrame({
        'PlateName': 'TEST',
        'IndexPlate': 'DI01',
        'FPrimer': [ADAPTER_F + seed for seed in f_primer_seed],
        'RPrimer': [ADAPTER_R + seed for seed in r_primer_seed],
        'VariableRegion': variable_region,
        'FrameDistance': f_dist,
        'BpIndStart': bp_idx,
        'AaIndStart': aa_idx,
    })

    # Construct new lists
    refseqs, new_frame, new_bp, new_aa = construct_ref_seq(default_refseq_df)

    def rev_comp(seq):
        comp_dict = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
        }
        return ''.join([comp_dict[base] for base in seq[::-1]])

    assert refseqs == [f_seed+var_reg+rev_comp(r_seed)
                       for f_seed, var_reg, r_seed
                       in zip(f_primer_seed, variable_region, r_primer_seed)]
    assert new_frame == [(dist + len(seed)) % 3
                         for dist, seed in zip(f_dist, f_primer_seed)]
    assert new_bp == [idx - len(seed)
                      for idx, seed in zip(bp_idx, f_primer_seed)]

    # Test from comments in bp_aa_pairs above; expectation = pair[1] values
    assert new_aa == [0, 0, 0, 0, 0, -1, 0, -1, -1]
    # Test from function
    assert new_aa == [idx - ((len(seed)+dist) // 3)
                      for idx, dist, seed
                      in zip(aa_idx, f_dist, f_primer_seed)]

def test_detailed_construct_ref_seq():
    pass
