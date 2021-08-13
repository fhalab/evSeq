"""
Test to confirm that evSeq/util/index_map.csv maps correctly to the
sequences listed in lib_prep_tools/evSeq_barcode_primer_seqs.csv
"""
from pathlib import Path
import os
import pandas as pd

import evSeq
from evSeq.util.globals import ADAPTER_F, ADAPTER_R
from evSeq.util.index_plate_mapper import index_plate_maker

# Define Nextera adapter seqs
NXT_I5_F = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
NXT_I7_R = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'

# Get the path to evSeq repo as root
root = str(Path(os.path.abspath(evSeq.__file__)).parents[1])
primer_seqs_path = os.path.join(root,
                            'lib_prep_tools',
                            'evSeq_barcode_primer_seqs.csv')
index_path = os.path.join(root, 'evSeq', 'util', 'index_map.csv')
    

def test_barcode_primer_file():
    """Confirms that sequences in barcode primer file are suitable
    for Nextera NGS prep and for evSeq analysis.
    """
    # Import seqs and names
    primer_df = pd.read_csv(primer_seqs_path)

    # Create empty list to check
    bad_seqs = []
    for seq, name in zip(primer_df['Sequence'].values,
                         primer_df['Name'].values):
        # Check forward seqs
        if name[-1] == 'f':
            if seq[:len(NXT_I5_F)] != NXT_I5_F:
                bad_seqs.append(f'{name}, Nextera_i5')
            if seq[-len(ADAPTER_F):] != ADAPTER_F:
                bad_seqs.append(f'{name}, evSeq_f')

        # Check reverse seqs
        if name[-1] == 'r':
            if seq[:len(NXT_I7_R)] != NXT_I7_R:
                bad_seqs.append(f'{name}, Nextera_i7')
            if seq[-len(ADAPTER_R):] != ADAPTER_R:
                bad_seqs.append(f'{name}, evSeq_r')
    
    assert not bad_seqs, f'Found bad adapters: {bad_seqs}'


def test_pairings():
    """Confirms no duplicated barcode pairs in index_map.csv"""
    pairings_df = pd.read_csv(index_path)

    pairings_df['Combo'] = pairings_df['FBC'] + pairings_df['RBC']
    pairings_df['Name'] = pairings_df['IndexPlate'] + '-' + pairings_df['Well']
    counts = pairings_df['Combo'].value_counts()
    bad_combos = counts[counts > 1].index.to_list()
    if bad_combos:
        bad_wells = pairings_df[pairings_df['Combo'].isin(bad_combos)]['Name']
        assert False, \
        f'Degenerate well pairings in index_map.csv! Check {list(bad_wells)}'


def test_mappings():
    """Confirms that primer combinations are mapped correctly to barcodes."""
    # Import data
    primer_df = pd.read_csv(primer_seqs_path)
    index_df = pd.read_csv(index_path)

    # If test_barcode_file() passed, create new column with barcodes
    try:
        test_barcode_primer_file()
    except AssertionError:
        raise AssertionError('Failed due to test_barcode_file() failure.')

    def trim_seqs_to_bcs(row):
        """Trims a barcode primer sequence to just its barcode"""
        # Trim forward seqs
        if row.Name[-1] == 'f':
            bc = row.Sequence[len(NXT_I5_F):-len(ADAPTER_F)]

        # Trim reverse seqs
        if row.Name[-1] == 'r':
            bc = row.Sequence[len(NXT_I7_R):-len(ADAPTER_R)]

        return bc

    primer_df['BC'] = primer_df.apply(trim_seqs_to_bcs, axis=1)

    # Set Name as index
    primer_df = primer_df.set_index('Name')

    # Assume wells are cycled the default way...
    cycle_df = index_plate_maker()
    cycle_df = cycle_df.set_index(['IndexPlate', 'Well'])

    def check_barcodes(index_row):
        """Apply with axis=1 to index_df"""
        # Get info
        plate = index_row.IndexPlate
        well = index_row.Well
        index_FBC = index_row.FBC
        index_RBC = index_row.RBC

        # Map Plate+Well to cycle_df
        f_well = cycle_df.loc[(plate, well), 'FBCSource']
        r_well = cycle_df.loc[(plate, well), 'RBCSource']

        # Get the primer barcodes (index = well)
        primer_FBC = primer_df.loc[f'{f_well}_f', 'BC']
        primer_RBC = primer_df.loc[f'{r_well}_r', 'BC']

        if (index_FBC == primer_FBC) and (index_RBC == primer_RBC):
            match = True
        else:
            match = False
        
        return match

    index_df['Match'] = index_df.apply(check_barcodes, axis=1)
    index_df['Name'] = index_df['IndexPlate'] + index_df['Well']
    
    bad_wells = index_df[~index_df['Match']]['Name'].to_list()
    if bad_wells:
        if len(bad_wells) == len(index_df):
            assert False, \
                f"""Barcode primers and index map do not match for ANY wells!
    This assumes index map has been generated in the default (row-wise shift)
    way to yield 8 dual-index plates. If your method was different, you should
    ignore this failure and check on your own.
    """
        elif len(bad_wells) > 10:
            n = len(bad_wells)
            assert False, \
                f"""Barcode primers and index map do not match for {n} wells!
    This assumes index map has been generated in the default (row-wise shift)
    way to yield 8 dual-index plates. If your method was different, you should
    ignore this failure and check on your own.
    """
        else:
            assert False, \
                f"""Barcode primers and index map do not match for wells:
    {bad_wells}
    """
