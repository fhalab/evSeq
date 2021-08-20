"""
Data visualization functions for evSeq output data.
"""

# Import required modules
import numpy as np
import pandas as pd
import os

import holoviews as hv
import colorcet as cc
import bokeh.io
from bokeh.layouts import row
from bokeh.models import HoverTool
import ninetysix as ns

from .util.logging import log_warning

# because the ninetysix package sets bokeh as the backend, we don't set it here
# otherwise we would need:
# hv.extension('bokeh')
hv.renderer('bokeh')

#### Heatmap ####
def generate_platemaps(max_combo_data, cmap=None):
    """Saves a plate heatmap html generated from from evSeq data.
    
    Input:
    ------
    max_combo_data: path (str) or DartaFrame
        Path to 'Combos_Coupled_Max.csv' from an evSeq experiment or
        a pandas DataFrame of that file.
    cmap: list-like or str, default None
        The colormap to use for the well outline indicating alignment
        frequency. If None, defaults to a Plasma-like (colorcet.bmy)
        colormap. If 'stoplight', uses a green-yellow-red colormap (not 
        the most colorblind friendly, but highly intuitive). Otherwise
        you may pass any list -like object containing four colors (e.g.,
        ['#337D1F', '#94CD35', '#FFC300', '#C62C20'] for 'stoplight').
    
    Returns:
    --------
    hm_holomap: an interactive Platemap
    """

    # Convert to dataframe if necessary
    if isinstance(max_combo_data, str):
        max_combo_df = pd.read_csv(max_combo_data)
    else:
        max_combo_df = max_combo_data.copy()
    
    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()
    
    # dictionary for storing plots
    hm_dict = {}

    # Generate plots for each plate
    for plate in unique_plates:
        
        # Split to just the information of interest
        df = max_combo_df.loc[max_combo_df.Plate == plate].copy()
    
        # generate a holoviews plot
        hm_dict[plate] = _make_platemap(df, title=plate, cmap=cmap)

    # make logseqdepth column
    max_combo_df['logseqdepth'] = np.log(
        max_combo_df['WellSeqDepth'], 
        out=np.zeros_like(
            max_combo_df['WellSeqDepth'], 
            dtype=float
        ),
        where=max_combo_df['WellSeqDepth'] != 0
    )

    # Set the center
    center = np.log(10)

    # Adjust if it is greater than max of data (avoids ValueError)
    if max_combo_df['logseqdepth'].max() <= center:

        # Log a warning
        log_warning(f"All wells associated with {plate} have a read depth <=10. "
                    "Be careful comparing heatmaps between this plate and others. "
                    "Be careful using this data; sequencing was not good.")

        # Adjust the center
        center = max_combo_df['logseqdepth'].median()

    # set color levels
    color_levels = ns.viz._center_colormap(df['logseqdepth'], center)

    # Uniform color levels
    for _hm in hm_dict.values():
        _hm.opts({'HeatMap': {'color_levels': color_levels}})
    
    # plot from the dictionary
    hm_holomap = hv.HoloMap(
        hm_dict, 
        kdims=['Plate']
    )

    return hm_holomap

def save_platemap_to_file(heatmaps, outputdir):
    
    file_path = os.path.join(outputdir, "Platemaps", "Platemaps")
    hv.renderer('bokeh').save(heatmaps, file_path)

def _make_platemap(df, title, cmap=None):
    """Generates a plate heatmap from evSeq data using Holoviews with
    bokeh backend.

    Called via `generate_platemaps`; see docs there.
    """

    # Convert SeqDepth to log for easier visualization.
    df['logseqdepth'] = np.log(
        df['WellSeqDepth'],
        out=np.zeros_like(df['WellSeqDepth'], dtype=float),
        where=(df['WellSeqDepth'] != 0)
    )

    # Create necessary Row and Column values and sort
    df['Row'] = df.apply(lambda row: row['Well'][0], axis=1)
    df['Column'] = df.apply(lambda row: int(row['Well'][1:]), axis=1)
    df = df.sort_values(['Column', 'Row'])
    df['Column'] = df['Column'].astype('str')

    # Set some base opts
    opts = dict(invert_yaxis=True, title=title, show_legend=True)

    # logseqdepth heatmap
    seq_depth_cmap = list(reversed(cc.CET_D9))

    # Set the center
    center = np.log(10)

    # Adjust if it is greater than max of data (avoids ValueError)
    if df['logseqdepth'].max() <= center:

        # Log a warning
        log_warning(f"All wells associated with {title} have a read depth <=10. "
                    "Be careful comparing heatmaps between this plate and others. "
                    "Be careful using this data; sequencing was not good.")

        # Adjust the center
        center = df['logseqdepth'].median()

    # center colormap
    color_levels = ns.viz._center_colormap(df['logseqdepth'], center)

    # Get heights
    n_rows = len(df['Row'].unique())
    n_cols = len(df['Column'].unique())
    height = int(50* n_rows)
    width = height * n_cols // n_rows

    # add tooltips
    tooltips = [
        ('Well', '@Well'),
        ('Mutations', '@VariantCombo'),
        ('Seq. depth', '@WellSeqDepth'),
        ('Align. freq.', '@AlignmentFrequency')
    ]

    hover = HoverTool(tooltips=tooltips)

    # generate the heatmap
    hm = hv.HeatMap(
        df,
        kdims=['Column', 'Row'],
        vdims=[
            'logseqdepth',
            'VariantCombo',
            'AlignmentFrequency',
            'WellSeqDepth',
            'Well'
        ],
    ).opts(
        **opts,
        colorbar=True,
        cmap=seq_depth_cmap,
        height=height,
        width=width,
        line_width=4,
        clipping_colors={'NaN': '#DCDCDC'},
        color_levels=color_levels,
        tools=[hover],
        colorbar_opts=dict(
            title='LogSeqDepth',
            background_fill_alpha=0
        )
    )

    # function to bin the alignment frequencies into more relevant groupings
    def bin_align_freq(value):
        if value > 0.95:
            bin_vals = '0.95+'
        if value <= 0.95 and value > 0.9:
            bin_vals = '0.90-0.95'
        if value <= 0.9 and value > 0.8:
            bin_vals = '0.80-0.90'
        # anything below 0.8 should really be discarded
        if value <= 0.8:
            bin_vals = '<0.80'

        return bin_vals
    
    # Bin alignment frequencies for easier viz
    bins = ['0.95+', '0.90-0.95', '0.80-0.90','<0.80']
    if cmap is None:
        cmap = [cc.bmy[int((1.1-i)*len(cc.bmy))]
                for i in [0.95, 0.9, 0.8, 0.4]]
    if 'stoplight' in cmap:
        cmap = ['#337D1F', '#94CD35', '#FFC300', '#C62C20']
    else:
        # Validate colormap
        if not isinstance(cmap, (list, tuple)):
            raise ValueError('cmap argument must be a list or tuple')
        if len(cmap) > 4:
            raise ValueError(
                'cmap argument has too many entries; only 4 should be passed'
            )
    cmap = {bin: color for bin, color in zip(bins, cmap)}

    # apply binning function to the AlignmentFrequency
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(
        bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows*1.2

    # alignment frequency heatmap for edges around wells
    boxes = hv.Points(
        df.sort_values(['AlignmentFrequency'], ascending=False),
        ['Column', 'Row'],
        'AlignmentFrequencyBinned'
    ).opts(
        **opts,
        marker='square',
        line_color='AlignmentFrequencyBinned',
        line_join='miter',
        cmap=cmap,
        line_width=6,
        fill_alpha=0,
        line_alpha=1,
        legend_position='right',
        size=box_size,
    )

    # USe in apply statement for residue labels
    def split_variant_labels(mutation_string):
        
        num_mutations = len(mutation_string.split('_'))

        if  num_mutations > 4:
            return str(num_mutations)+' muts'

        mutation_string = mutation_string.replace('?','')
        new_line_mutations = mutation_string.replace('_','\n')
        
        return new_line_mutations
    
    _df = df.copy()
    _df['Labels'] = _df['VariantCombo'].apply(split_variant_labels)

    # Set the font size based on if #PARENT# is in a well and num of mutations
    max_num_mutations = _df['Labels'].apply(lambda x: len(x.split('\n'))).max()
    has_parent = ('#PARENT#' in _df['Labels'])
    
    if max_num_mutations > 3 or has_parent:
        label_fontsize = '8pt'
    else:
        label_fontsize = '10pt'

    labels = hv.Labels(
        _df,
        ['Column', 'Row'],
        'Labels',
    ).opts(
        text_font_size=label_fontsize,
        **opts
    )

    # return formatted final plot
    return (hm*boxes*labels).opts(frame_height=550,
                                  frame_width=550 * 3 // 2,
                                  border=50,
                                  show_legend=True)

#### Read quality chart ####
def generate_qualplot(seqpairs, output_dir):
    """Makes histograms of read qualities and saves to designated
    location.
    """
    
    # Generate forward and reverse qualities
    all_qualities = [seqpair.read_quals() for seqpair in seqpairs]
    f_qual_counts = np.unique(
        [int(qual[0]) for qual in all_qualities if not np.isnan(qual[0])],
        return_counts=True
    )
    r_qual_counts = np.unique(
        [int(qual[1]) for qual in all_qualities if not np.isnan(qual[1])],
        return_counts=True
    )
        
    # Plot counts
    p_f = plot_read_qual(f_qual_counts).opts(title='Forward Read Quality')
    p_r = plot_read_qual(r_qual_counts).opts(title='Reverse Read Quality')
    
    # Render to bokeh and combine into single chart
    p = row(
        hv.render(p_f),
        hv.render(p_r)
    )

    # Output as html
    filename = os.path.join(output_dir, "Qualities", "QualityPlot.html")
    bokeh.io.output_file(filename)
    bokeh.io.save(p)

def plot_read_qual(counts):
    """Given an array, creates a holoviews histogram."""
    p = hv.Histogram(counts).opts(
        xlabel='Mean quality score of sequence',
        ylabel='Counts',
        xlim=(0, 40),
        height=300,
        width=400,
        yformatter='%f',
    )

    return p

### SEQUENCE-FUNCTION PLOTTING ###
def combine_seq_func_data(
    data_df,
    evseq_output_dir
):
    """Function to combine evSeq data with function data to generate
    sequence-function pairs.

    Inputs:
    -------
    data_df: DataFrame
        DataFrame containing function/activity data.
    evseq_output_dir: path (str)
        Path to the location of the evSeqOutput directory corresponding
        to the variants used to obtain function data.

    Returns:
    --------
    merged_data: DataFrame
        Cleaned, merged sequence-function data containing all columns
        from data_df and the evSeq output file.
    """

    # If the well column of the data_df is not zero-padded, pad it
    data_df['Well'] = data_df['Well'].apply(ns.parsers.pad)

    # Read in the results of evSeq
    seq_df = pd.read_csv(evseq_output_dir+'AminoAcids_Decoupled_Max.csv')
    
    # Find wells that have multiple mutations since this only works for
    # single-mutant libraries
    count_df = pd.DataFrame(seq_df.groupby(
        ['Plate', 'Well'])['AaPosition'].count())
    mut_count_dict = count_df.to_dict()['AaPosition']

    # Get the number of mutations in a sequence for later filtering
    seq_df['MutCount'] = seq_df.apply(
        lambda x: mut_count_dict[x['Plate'], x['Well']], 
        axis=1
    )

    # Rename columns for downstream use
    seq_df = seq_df.rename(
        columns={'AaPosition': 'Position', 'Aa': 'AA'}
    )
    
    # Remove dead sequencing wells
    seq_df = seq_df.loc[(seq_df['Flag'] != '#DEAD#')].copy()

    # Remove sequencing wells with a flag
    # allow wells w/ 'Unexpected Variation' if MutCount is 1
    inds = (seq_df['Flag'].isna()) | ((seq_df['Flag'] == 'Unexpected Variation') & (seq_df['MutCount'] < 2))

    seq_df = seq_df.loc[inds].copy()

    # Now that the #DEAD# wells are removed, set to int
    seq_df = seq_df.astype({'Position': 'int64'})

    # Check for unexpected shared columns between seq_df and data_df
    shared_cols = set(seq_df.columns).intersection(set(data_df.columns))
    expected_cols = {'Position', 'Well', 'Plate'}

    if shared_cols != expected_cols:
        if len(shared_cols) > len(expected_cols):
            raise AssertionError(
                "data_df has too many shared columns. "
                "You should have 'Position', 'Well', and 'Plate'. "
                f"Your shared columns are: {shared_cols}"
            )

        elif len(shared_cols) < len(expected_cols):
            raise AssertionError(
                "data_df has too few shared columns. "
                "You should have 'Position', 'Well', and 'Plate'. "
                f"Your shared columns are: {shared_cols}"
            )

        else:
            raise AssertionError(
                "data_df and seq_df may share the wrong columns. "
                "You should have 'Position', 'Well', and 'Plate'. "
                f"Your shared columns are: {shared_cols}"
            )
    
    # Merge data_df and seq_df on shared columns and fill missing w/ NaN
    merged_data = data_df.merge(seq_df, how='outer')
    return merged_data

def _make_SSM_activity_plot(
    df, 
    value,
    missing,
    activity_range,
    title,
    sort_counter,
    center_cmap,
    known,
    standard,
    jitter,
    unknown_jitter,
):
    """Plots activities of each mutation assuming the data corresponds
    to a site-saturation mutagenesis (SSM) library.

    Called via `plot_SSM_activities`; see docs there.
    """

    # list of AAs
    AAs = list('ACDEFGHIKLMNPQRSTVWY')

    # Copy and sort the sequence/function data
    temp_df = df.copy()    
    temp_df = temp_df.sort_values('Sort')
    
    # Opts for plotting
    opts = dict(
        width=600,
        height=500,
        xrotation=45
    )

    # Set up jitter
    if jitter is None:
        jitter = 0
    if unknown_jitter is None:
        unknown_jitter = jitter

    # Create tooltips for plotting
    tooltips = [
        ('Well', '@Well'),
        ('Amino Acid', '@AA'),
        ('Seq. depth', '@WellSeqDepth'),
        ('Align. freq.', '@AlignmentFrequency'),
        ('Mutations', '@MutCount'),
        # Can't do f-string here, weird syntax for tooltips
        (value, '@{'+value+'}')
    ]
    hover = HoverTool(tooltips=tooltips)

    # Plot the data for which the AA is known
    _df = temp_df[temp_df['Residue'] != 'Unknown'].copy()
    
    p = ns.viz.plot_bar(
        _df,
        'Residue',
        value_name=value,
        color=value,
        cmap='coolwarm',
        jitter=jitter,
    ).opts(
        {'Scatter': dict(tools=[hover])}
    )

    # Plot the data w/ an unknown AA
    _df = temp_df[temp_df['Residue'] == 'Unknown'].copy()
    
    # Create an empty DataFrame if there are no unknowns
    if len(_df) != 0:
        unknowns = True
    else:
        unknowns = False
        _df = pd.DataFrame({
            'Residue': ['Unknown'],
            value: [np.nan],
            'Sort': sort_counter
        })
    
    # Plot unknown activities if they exist and overlay with knowns
    if unknowns:
        p_unknown = ns.viz.plot_bar(
            _df,
            'Residue',
            value_name=value,
            color=value,
            cmap='coolwarm',
            jitter=unknown_jitter,
        ).opts({
            'Bars': dict(line_alpha=0, fill_alpha=0),
            'Scatter': dict(tools=[hover])
        })

        p = hv.Overlay([p, p_unknown])

    # Place 'n.d.' for missing amino acids
    # Get the height for the text
    if activity_range is None:
        max_value = temp_df[value].max()
    else:
        max_value = activity_range[1]
    height = max_value / 20

    # Plot the 'n.d.' in the empty sections
    for AA in AAs:
        if AA in missing:
            p = p*hv.Text(AA, height, 'n.d.', rotation=90)
        else:
            p = p*hv.Text(AA, height, ' ', rotation=90)
    
    # Plot 'Unknown' if there are none!
    if not unknowns:
        p = p*hv.Text('Unknown', height*1.25, 'none', rotation=90)

    # Add opts for the size of the plot
    p = p.opts(**opts)

    # Relabel from AA
    if title is not None:
        if not isinstance(title, str):
            raise ValueError(
                f'Title must be of type str, you passed {type(title)}.')
        p.opts(title=title)
    if activity_range is not None:
        p.opts(ylim=activity_range)

    if center_cmap:
        center = temp_df[temp_df[known] == standard][value].mean()
        color_levels = ns.viz._center_colormap(
            temp_df[value].dropna(), center)
        p.opts(
            {'Bars': {'color_levels':color_levels}}
        )

    return p

def _make_count_plot(
    temp_df,
    missing,
    sort,
    counts_range
):
    """Creates a bar plot of the number of occurences of each mutation
    in a site-saturation mutagenesis (SSM) library.

    Called via `plot_SSM_activities`; see docs/code there.
    """
    
    # Get observed AA counts
    counts = temp_df['Residue'].value_counts()

    # Create DataFrame of counts, accounting for missing values
    count_dict = {
        'Residue': list(counts.index),
        'Counts': [value if AA not in missing else 0
                    for value, AA in zip(counts.values, counts.index)],
    }

    for AA in missing:
        count_dict['Residue'].append(AA)
        count_dict['Counts'].append(0)

    df_counts = pd.DataFrame(count_dict)

    # Sort meaningfully, again...
    df_counts['Sort'] = df_counts['Residue'].replace(sort)
    # return df_counts
    df_counts = df_counts.sort_values('Sort')

    # Make chart
    p_counts = hv.Bars(
        df_counts,
        'Residue',
        'Counts',
    ).opts(
        width=600,
        height=200,
        xrotation=45,
        color='#DCDCDC'
    )

    if counts_range is not None:
        p_counts = p_counts.opts(ylim=counts_range)
    
    return p_counts


def plot_SSM_activities(
    seq_func_data,
    value,
    min_align_freq=0.8,
    min_seq_depth=10,
    max_muts=1,
    counts=True,
    title=None,
    center_cmap=True,
    activity_range=None,
    counts_range=None,
    known=None,
    variant=None,
    standard=None,
    jitter=None,
    unknown_jitter=None,
):
    """Takes a DataFrame of single mutant sequence/fitness data (see
    `combine_seq_func_data`) and plots the activities of all sequences
    as well as a frequency plot of each mutant.

    Inputs:
    -------
    seq_func_data: DataFrame
        Output of (or of the same form as) `combine_seq_func_data`.
        Effectively a list of sequence-function pairs + extra info.
    value: str
        The name of the column that corresponds to the function data.
    min_align_freq: float, default 0.8
        The minimum alignment frequency of the variants to plot. You
        will likely want to increase this, not lower it. If it is
        lowered, you will have fewer in the 'Unknown' column at the
        expense of accuracy, as mixed populations will be more likely to
        count for the plotted activity values. You can hover and
        discover points with bade alignment frequency, however.
    min_seq_depth: int, default 10
        The minimum sequencing depth (non-log scale) of the variants to
        plot. You will liekly want to increase this, not lower it.
    max_muts: int, default 1
        The maximum allowable mutations found in the variant. For single
        SSM libraries you should only have one. Can be increased if your
        mutations are synonymous. This does not apply to known controls,
        since the code assumed these are not variable. If you have known
        bad controls, filter these beforehand.
    counts: bool, default True
        Whether or not to create a bar plot for the observations of each
        amino acid (the variant counts) below the main activity plot.
    title: str, default None
        Title of the plot.
    center_cmap: bool, default True
        Whether or not to center the colormap around `standard`. This
        will switch to False if a `standard` is not passed.
    activity_range: tuple, default None
        The range of the activity plot. Defaults to (0, max(data)*1.2).
        Useful to check when data has some values below 0, but bars will
        not look as nice. Keep 0 as min unless your negative values are
        significant/important.
    counts_range: tuple, default None
        The range of the counts bar plot. Same as `activity_range` but
        negative values are not possible.
    known: str, default None
        The name of the column that contains control information. E.g.,
        if a column 'Controls' specifies Parent v. Negative v. Variant,
        known='Controls'. Requires that the `variant` argument also be
        filled out.
    variant: str, default None
        If `known` is not None, then this value specifies the name of
        the variant entries in the `known` column. E.g., from the 
        example in `known` above, variant='Variant'.
    standard: str, default None
        Similar to `variant`, but not necessary (although very useful).
        If you have a positive control/parent protein standard for
        comparison, this string specifies the rows containing these. 
        E.g., from the example in `known` above, standard='Parent'.
    jitter: float [0,1], default None
        How much jitter to add to the points. Defaults to zero.
    unknown_jitter: float [0,1], default None
        How much jitter to add to the points in the Unknown column,
        since there may be more there than any other column. If `jitter`
        is not None and this is None, defaults to the value of `jitter`.

    Returns:
    --------
    activity plot: Interactive Holoviews HoloMap
    """
    
    df = seq_func_data.copy()

    # Remove rows w/o a known activity
    df = df.loc[df[value].notna()]

    # Add a 'Residue' column for tracking AAs and controls together
    df['Residue'] = df['AA'].copy()

    # list of AAs
    AAs = list('ACDEFGHIKLMNPQRSTVWY*')

    # Set unknowns vs controls
    def set_knowns(df):
        
        # Filter based on alignment frequency and sequencing depth
        if df['AlignmentFrequency'] < min_align_freq:
            df['Residue'] = 'Unknown'
        elif df['WellSeqDepth'] < min_seq_depth:
            df['Residue'] = 'Unknown'

        # If the AA is NaN then set to unknown
        elif pd.isna(df['AA']):
            df['Residue'] = 'Unknown'

        # If a sequence has unacceptable number of mutations set to unknown
        elif df['MutCount'] > max_muts:
            df['Residue'] = 'Unknown'

        # Set controls if a type column passed
        if known is not None:
            if df[known] == variant:
                pass
            elif df[known] == standard:
                df['Residue'] = standard_string
            else:
                df['Residue'] = df[known]
        return df
    
    # Set up dictionary for sorting meaningfully within the plotting loop
    sort = {}

    # If user passes a standard, use it as key for value 0 so that it
    # displays first on the plot
    if standard is not None:
        standard_string = standard
        sort[standard_string] = 0
        
    # Default behavior given a standard is to center the cmap
    if standard is None:
        center_cmap = False

    # Set the sort counter so that controls follow the variant library
    sort_counter = 22

    # If a column name for data type labels is passed, add the
    # other controls/values in it to the sort dict in order of appearance.
    # The order will be the same for all subplots of the HoloMap
    if known is not None:
        if variant is None:
            raise ValueError(
                'Unspecified `variant` argument. \n\n'
                'The `variant` argument should not be None if `known` '\
                'is not None. Of the following values in your `known` '\
                f'column "{known}", which one corresponds to the wells'\
                f' with protein variants? \n\n'\
                f'Options: {df[known].unique()}'
            )
        for _type in df[known].unique():
            if _type not in [standard, variant]:

                sort[_type] = sort_counter
                sort_counter += 1

    # Add an 'Unknown' sort as the last item in the sort
    sort['Unknown'] = sort_counter

    # Set known variants/AAs in sort: 1-21
    sort.update({AA: i for i, AA in enumerate(AAs, 1)})
    
    # Initialize dictionaries for storing activity and AA count plots
    activity_plot_dict = {}
    count_plot_dict = {}

    # Loop through the plates slice out one plate at a time
    for plate in df['Plate'].unique():
        temp_df = df.loc[df['Plate'] == plate].copy()

        # Loop through the positions that appear on this plate and slice out
        for position in temp_df['Position'].unique():
            temp_df = temp_df.loc[temp_df['Position'] == position].copy()

            # Requires standard_string and sets controls, variants, and unknowns
            temp_df = temp_df.apply(set_knowns, axis=1)

            # Find which AAs are missing from this plate/position
            missing = tuple(set(AAs) - set(temp_df['Residue'].unique()))

            # Loop through the missing AAs and add them to the DataFrame
            # so they still appear (empty) in the bar-plot
            for AA in missing:

                # Append row for the missing AA
                temp_df = temp_df.append(
                    pd.DataFrame(
                        [[None] * (len(temp_df.columns))],
                        columns=temp_df.columns
                    ),
                    ignore_index=True
                )

                # Set the AA value for new row to the missing AA value
                temp_df.at[len(temp_df) - 1, 'Residue'] = AA           

            # Make the sort column based on the sort dictionary made
            # outside of the loop (so axes on subplots match)
            temp_df['Sort'] = temp_df['Residue'].replace(sort)

            #### Activity Plot ####
            p = _make_SSM_activity_plot(
                df=temp_df,
                value=value,
                missing=missing,
                activity_range=activity_range,
                title=title,
                sort_counter=sort_counter,
                center_cmap=center_cmap,
                known=known,
                standard=standard,
                jitter=jitter,
                unknown_jitter=unknown_jitter,
            )

            ##### Histogram #####
            if counts:
                p_counts = _make_count_plot(
                    temp_df,
                    missing,
                    sort,
                    counts_range=counts_range,
                )

                # Since we are plotting the histogram, remove xlabel from activity plot
                p.opts(xlabel='')
                
                # Store the count histograms in a dictionary
                count_plot_dict[plate + ': ' + str(position)] = p_counts
            
            # Store the activity plots in a histogram
            activity_plot_dict[f"{plate}: {position}"] = p

    # Generate a HoloMap of the activity plots
    activity_hmap = hv.HoloMap(
        activity_plot_dict,
        kdims=['Plate: Position']
    ).opts(
        {'Bars': {'framewise': True}}
    )

    if counts:
        # Generate a HoloMap of the count histograms
        count_hmap = hv.HoloMap(
            count_plot_dict,
            kdims=['Plate: Position']
        )
        return (activity_hmap+count_hmap).cols(1)

    else:
        return activity_hmap


def check_distributions(
    df,
    bins=50,
    violin=True,
):
    """
    Generates visualizations of the distribution(s) of sequencing
    depth for a given set of evSeq data.
    
    Inputs:
    -------
    df: pandas DataFrame
        Any evSeq OutputCounts file (must contain WellSeqDepth column)
    bins: int, default 50
        Number of bins for the histogram plot.
    violin: bool, default True
        Whether or not to plot a Violin plot with distributions
        separated by Plate to check for plate-specific differences.
        
    Returns:
    --------
    holoviews plot of distribution(s)
    """

    # Grab only non-dead variants
    for column in ('VariantCombo', 'Bp', 'Aa'):
        if column in df.columns:
            clean_depths = df.loc[df[column] != '#DEAD#']['WellSeqDepth']
            continue

    # Bin the data
    binned = np.histogram(
        clean_depths,
        bins=bins,
    )

    # Create a histogram
    p = hv.Histogram(
        binned,
        kdims='Sequencing Depth',
        vdims='Counts',
    ).opts(
        color='#DCDCDC',
    )

    # Plot the median as a vertical line
    p = p*hv.VLine(
        clean_depths.median()
    ).opts(
        color='black',
        line_width=2,
    )

    if violin:

        # Look at the per-plate distributions on a Violin chart to check for outliers
        p_viol = hv.Violin(
            df.loc[df['VariantCombo'] != '#DEAD#'],
            kdims='Plate',
            vdims='WellSeqDepth',
        ).opts(
            violin_fill_color='#DCDCDC',
            xrotation=45,
        )

        p = p + p_viol

    return p
