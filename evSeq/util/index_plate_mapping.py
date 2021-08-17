"""
Code for mapping barcode plate sequences to index plates.
"""
from pathlib import Path

import pandas as pd
from evSeq.util.globals import ADAPTER_F, ADAPTER_R

# Nextera adapters
NXT_I5_F = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
NXT_I7_R = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'


def _validate_index_inputs(
    plate_prefix,
    plate_start,
    axis,
    axis_list,
    alt_axis_list,
    FBC_map,
    RBC_map,
    axis_start,
    axis_end,
    n_digits,
):
    """Validates inputs to index_mapping()"""
    if not isinstance(plate_prefix, str):
        raise ValueError('plate_prefix argument must be a string.')

    if not isinstance(plate_start, int):
        raise ValueError('plate_start argument must be int')

    if axis not in ('row', 'column', 0, 1, None):
        raise ValueError(
            'axis argument must be either "row" (or 0), "column" (or 1), or '
            'None (if RBC = stamp).'
        )

    if axis is not None and axis_list is not None:
        if not hasattr(axis_list, '__iter__'):
            raise ValueError(
                'axis_list argument is not iterable.'
            )

    if axis is None:
        axis_list = None

    if alt_axis_list is not None:
        if not hasattr(alt_axis_list, '__iter__'):
            raise ValueError(
                'alt_axis_list argument is not iterable.'
            )

        # Confirm string
        alt_axis_list = [str(i) for i in alt_axis_list]

    if alt_axis_list is None:
        if (axis_list is None) and (axis in ('row', 0)):
            # Default to columns of 96-well plate
            alt_axis_list = [str(i) for i in range(1, 12+1)]

    if axis_list is None:
        if axis in ('row', 0):
            axis_list = 'ABCDEFGH'
        if axis in ('column', 1):
            axis_list = [str(i) for i in range(1, 12+1)]

    axis_len = len(axis_list) if axis_list is not None else 0

    if not isinstance(axis_start, int):
        raise ValueError('start argument must be int')

    if not isinstance(n_digits, int):
        raise ValueError('n_digits argument must be int')

    if axis_end is None:
        axis_end = axis_start + axis_len
    else:
        axis_list = axis_list[axis_start:axis_end]
        axis_len = len(axis_list)

    if not isinstance(axis_end, int):
        raise ValueError('end argument must be int')

    plate_end = len(str(plate_start + axis_len))
    if plate_end > n_digits:
        n_digits = plate_end

    if FBC_map != 'stamp':
        if len(FBC_map.split('->')) != 2:
            raise ValueError(
                'FBC_map argument incorrectly specified. See the docstring.'
            )
        f_src, f_dst = FBC_map.split('->')

        if (f_src not in axis_list) and (f_src != 'N'):
            raise ValueError(
                'Cannot find FBC_map source in axis_list.'
            )

        if f_src == f_dst:
            f_src, f_dst = ([axis_list.index(f_src)]*axis_len,)*2
        elif f_src != 'N':
            f_src = [axis_list.index(f_src)]*axis_len
            f_dst = [axis_start+i for i, _
                     in enumerate(_cycle(axis_list)[axis_start:axis_end])]
        # Swap things if the opposite is given
        elif f_src == 'N':
            f_dst = [axis_list.index(f_dst)]*axis_len
            f_src = [axis_start+i for i, _
                     in enumerate(_cycle(axis_list)[axis_start:axis_end])]
        else:
            raise NotImplementedError('Mapping format not supported.')

    else:
        f_src, f_dst = [0]*axis_len, [0]*axis_len

    if RBC_map == 'stamp':
        if (FBC_map == 'stamp') and (axis is not None):
            raise ValueError(
                'Only one index plate can be made with these inputs. axis must'
                ' be None or RBC_map arg must be changed.'
            )

    if RBC_map != 'stamp':
        if len(RBC_map.split('->')) != 2:
            raise ValueError(
                'RBC_map argument incorrectly specified. See the docstring.'
            )
        r_src, r_dst = RBC_map.split('->')

        if (r_src not in axis_list) and (r_src != 'N'):
            raise ValueError(
                'Cannot find RBC_map source in axis_list.'
            )

        if r_src == r_dst:
            r_src, r_dst = ([axis_list.index(r_src)]*axis_len,)*2
        elif r_src != 'N':
            r_src = [axis_list.index(r_src)]*axis_len
            r_dst = [axis_start+i for i, _
                     in enumerate(_cycle(axis_list)[axis_start:axis_end])]
        # Swap things if the opposite is given
        elif r_src == 'N':
            r_dst = [axis_list.index(r_dst)]*axis_len
            r_src = [axis_start+i for i, _
                     in enumerate(_cycle(axis_list)[axis_start:axis_end])]
        else:
            raise NotImplementedError('Mapping format not supported.')

    else:
        r_src, r_dst = [0]*axis_len, [0]*axis_len

    returned_args = (
        axis_list,
        alt_axis_list,
        f_src,
        f_dst,
        r_src,
        r_dst,
        n_digits,
    )

    return returned_args


def _cycle(cycle):
    # Interesting note: time for cycle*1e3 ~= cycle*1e4 (~300 ns),
    # while cycle*1e5 is ~10x slower than either (~3 µs)
    return cycle*1000


def index_plate_maker(
    plate_prefix='DI',
    plate_start=1,
    axis='row',
    axis_list=None,
    alt_axis_list=None,
    FBC_map='stamp',
    RBC_map=None,
    axis_start=0,
    axis_end=None,
    n_digits=2,
    hide_wells=False,
):
    """Code to generate a general index plate map, providing directions
    for how to arrange forward and reverse primer stock wells for a
    number of dual-index plates. The most useful ways to use this are:
        1. Run with only the axis='column' argument, to generate 12
        index plates.
        2. Change this axis_list/alt_axis_list arguments, for using non-
        96-well barcode plates.

    Inputs:
    -------
    plate_prefix: str, default 'DI' for 'dual-index'
        Prefix for the plate name.
    plate_start: int, default 1
        Number appended to `plate_prefix`.
    axis: 'row', 'column', or 0, 1, or None, default 'row'
        The axis along the plate in which the wells are cycled.
    axis_list: iterable, default None
        Either row-like ('ABC...') or column-like ([1, 2, 3, ...]). If
        None, defaults to axis argument for a 96-well plate.
    alt_axis_list: iterable, default None
        A list of values for the alternate axis, e.g., if `axis='row'`,
        this should be a list of the columns: [1, 2, 3, etc.]
    FBC_map: str, default 'stamp'
        Must be 'stamp' or of form '{source index}->{destination index}',
        where one of the index values is 'N' and the other is in the
        `axis_list`, e.g., 'A->N' or 'N->1'.
    RBC_map: str, default 'A->N'
        Must be 'stamp' or of form '{source index}->{destination index}',
        where one of the index values is 'N' and the other is in the
        `axis_list`, e.g., 'A->N' or 'N->1'.
    axis_start: int, default 0
        Where in the axis_list to start (alternative use a different
        axis argument).
    axis_end: int, default None
        Similar to axis_start.
    n_digits: int, default 2
        For padding plate names; if 2, plates with numbers in the single
        digits (e.g., plate 1) will be written with zero to pad the name
        as in `DI01`, etc.
    hide_wells: bool, default False
        Whether or not to hide the 'Well' column.

    Returns:
    --------
    df: pd.DataFrame
        Contains the columns ['IndexPlate', 'FBCSource', 'RBCSource']
        as well as ['Well'] if `alt_axis_list` is given/assumed. The
        'Well' column is needed for downstream generate_index_map()
        function, but without it (`hide_wells=True`) it is a good 
        overview/guide.
    """
    # Set some defaults
    if axis in ('row', 0):
        if (axis_list is None) and (RBC_map is None):
            RBC_map = 'A->N'
            alt_axis_list = [str(i) for i in range(1, 13)]

    if axis in ('column', 1):
        if (axis_list is None) and (RBC_map is None):
            RBC_map = '1->N'
            alt_axis_list = 'ABCDEFGH'

    # Validate inputs
    (axis_list, alt_axis_list, f_src, f_dst, r_src, r_dst, n_digits) \
        = _validate_index_inputs(
            plate_prefix,
            plate_start,
            axis,
            axis_list,
            alt_axis_list,
            FBC_map,
            RBC_map,
            axis_start,
            axis_end,
            n_digits
    )

    # Create plate names
    n_plates = len(f_src)
    plate_end = plate_start + n_plates
    plate_strings = [str(i) for i in range(plate_start, plate_end)]
    plate_strings = [f'0{i}' if len(
        i) < n_digits else i for i in plate_strings]
    plate_names = [f'{plate_prefix}{string}' for string in plate_strings]

    # Create a list to store the DataFrames
    dfs = []

    # For each plate
    for i, plate in enumerate(plate_names):

        # The source/destination lists are the *offsets* for a given plate
        f_src_offset = f_src[i]
        f_dst_offset = f_dst[i]
        r_src_offset = r_src[i]
        r_dst_offset = r_dst[i]

        # For each value along the given axis, grab the appropriate values given the offset
        # E.g., for each row, adjust the source/destination to the appropriate row for that plate
        working_f_src = [_cycle(axis_list)[i+f_src_offset]
                         for i, _ in enumerate(axis_list)]
        working_f_dst = [_cycle(axis_list)[i+f_dst_offset]
                         for i, _ in enumerate(axis_list)]
        working_r_src = [_cycle(axis_list)[i+r_src_offset]
                         for i, _ in enumerate(axis_list)]
        working_r_dst = [_cycle(axis_list)[i+r_dst_offset]
                         for i, _ in enumerate(axis_list)]

        # Create a source-destination DataFrame for the Forward Barcodes
        f_df = pd.DataFrame({
            'IndexPlate': plate,
            'Destination': working_f_dst,
            'FBCSource': working_f_src,
        })

        # Do the same for the Reverse Barcodes
        r_df = pd.DataFrame({
            'IndexPlate': plate,
            'Destination': working_r_dst,
            'RBCSource': working_r_src,
        })

        # Merge these (along Destination) to create a FBC+RBC source pair for a given destination
        working_df = f_df.merge(r_df)

        # Add to the list
        dfs.append(working_df)

    # Create final full DataFrame
    df = pd.concat(dfs).reset_index(drop=True)

    # If alt_axis_list is given, use this to fill in all wells
    if (alt_axis_list is not None) and (not hide_wells):
        if axis in ('row', 0):
            axis = 'row'
            rows = axis_list
            columns = alt_axis_list
        else:
            axis = 'column'
            rows = alt_axis_list
            columns = axis_list

        wells = [f'{row}{column}' for row in rows for column in columns]
        well_df = pd.DataFrame({
            'row': [well[0] for plate in plate_names for well in wells],
            'column': [well[1:] for plate in plate_names for well in wells]
        }).rename(columns={axis: 'Destination'})

        merged = df.merge(well_df, on='Destination').drop_duplicates()

        def pad(x): return str(x) if len(str(x)) == 2 else f'0{x}'

        if axis == 'row':
            merged['Well'] = merged['Destination'] + \
                merged['column'].apply(pad)
            merged['FBCSource'] = merged['FBCSource'] + \
                merged['column'].apply(pad)
            merged['RBCSource'] = merged['RBCSource'] + \
                merged['column'].apply(pad)
        else:
            merged['Well'] = merged['row'] + merged['Destination'].apply(pad)
            merged['FBCSource'] = merged['row'] + \
                merged['FBCSource'].apply(pad)
            merged['RBCSource'] = merged['row'] + \
                merged['RBCSource'].apply(pad)

        # Only take relevant columns
        df = merged[
            ['IndexPlate', 'Well', 'FBCSource', 'RBCSource']
        ].sort_values(['IndexPlate', 'Well']).reset_index(drop=True)

    return df


def generate_index_map(
    mapping,
    barcode_plate_seqs,
    NGS_adapter_f=NXT_I5_F,
    NGS_adapter_r=NXT_I7_R,
):
    """Generates a new index map in the style of evSeq/util/index_map.csv
    in case new barcode layouts are needed. The simplest and safest way
    to use this is to generate column-wise dual-index plates, e.g.,
    making DI01 – DI12 where the standard evSeq barcode primers (found
    at lib_prep_tools/evSeq_barcode_primer_seqs.csv) are cycled by
    columns as opposed to rows. This will reduce sequencing depth, so it
    is not recommended unless you are sure your sequencing depth will
    still be sufficient with 12 plates as opposed to 8.

    Inputs:
    -------
    mapping: pd.DataFrame or path to csv file
        Path to a mapping file, as generated by cycle_plate_axis().
    barcode_plate_seqs: path to csv file
        Path to a barcode primer sequences, with sequences entered the
        same way as in 'lib_prep_tools/evSeq_barcode_primer_seqs.csv'.
    NGS_adapter_f: str DNA sequence, default Nextera i5
        Which forward adapter is used for NGS chemistry. Should not need
        to change this.
    NGS_adapter_r: str DNA sequence, default Nextera i7
        Which reverse adapter is used for NGS chemistry. Should not need
        to change this.

    Returns:
    --------
    index_map: pd.DataFrame
        New index_map data
    """

    if isinstance(mapping, str):
        mapping = pd.read_csv(mapping)

    if isinstance(barcode_plate_seqs, str):
        primer_df = pd.read_csv(barcode_plate_seqs)
    else:
        primer_df = barcode_plate_seqs.copy()

    def trim_seqs_to_bcs(row):
        """Trims a barcode primer sequence to just its barcode"""
        # Trim forward seqs
        if row.Name[-1] == 'f':
            bc = row.Sequence[len(NGS_adapter_f):-len(ADAPTER_F)]

        # Trim reverse seqs
        if row.Name[-1] == 'r':
            bc = row.Sequence[len(NGS_adapter_r):-len(ADAPTER_R)]

        return bc

    primer_df['BC'] = primer_df.apply(trim_seqs_to_bcs, axis=1)

    # Set Name as index
    primer_df = primer_df.set_index('Name')

    def pull_barcodes(row):
        """Apply with axis=1 to mapping"""
        # Get info
        plate = row.IndexPlate
        well = row.Well
        mapping_FBC = row.FBCSource
        mapping_RBC = row.RBCSource

        # Get the primer barcodes (index = well)
        primer_FBC = primer_df.loc[f'{mapping_FBC}_f', 'BC']
        primer_RBC = primer_df.loc[f'{mapping_RBC}_r', 'BC']

        return primer_FBC, primer_RBC

    mapping['FBC'], mapping['RBC'] = zip(*mapping.apply(pull_barcodes, axis=1))

    index_df = mapping.reset_index()[['IndexPlate', 'Well', 'FBC', 'RBC']]

    return index_df


def save_csv(df, filename):
    """Saves a csv of the input DataFrame with the given filename."""

    df.to_csv(filename, index=False)
    path = str(Path(filename).resolve())
    print(f"DataFrame saved at '{path}'.")


def check_barcode_pairings(
    index_map,
    FBC_col=None,
    RBC_col=None,
):
    """Confirms no duplicated barcode pairs in index_map.csv"""
    if isinstance(index_map, str):
        index_df = pd.read_csv(index_map)
    else:
        index_df = index_map.copy()

    # Guess columns
    if FBC_col is None:
        if 'FBC' in index_df.columns:
            FBC_col = 'FBC'
        elif 'FBCSource' in index_df.columns:
            FBC_col = 'FBCSource'
    if FBC_col not in index_df.columns:
        raise ValueError(
            'Could not find column corresponding to Forward Barcode (FBC_col)'
        )
    if RBC_col is None:
        if 'RBC' in index_df.columns:
            RBC_col = 'RBC'
        elif 'RBCSource' in index_df.columns:
            RBC_col = 'RBCSource'
    if RBC_col not in index_df.columns:
        raise ValueError(
            'Could not find column corresponding to Reverse Barcode (RBC_col)'
        )

    index_df['Combo'] = index_df[FBC_col] + index_df[RBC_col]
    index_df['Name'] = index_df['IndexPlate'] + '-' + index_df['Well']
    counts = index_df['Combo'].value_counts()
    bad_combos = counts[counts > 1].index.to_list()
    if bad_combos:
        bad_wells = index_df[index_df['Combo'].isin(bad_combos)]['Name']
        assert False, \
            f'Degenerate well pairings in index_map.csv! Check {list(bad_wells)}'
