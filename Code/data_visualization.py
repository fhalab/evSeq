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

from .logging import log_warning

hv.extension('bokeh')
hv.renderer('bokeh')


#### Heatmap ####
def generate_sequencing_heatmaps(max_combo_df):
    """Saves a heatmap html generated from from ssSeq data."""
    
    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()
    
    # dictionary for storing plots
    hm_dict = {}

    # Generate plots for each plate
    for plate in unique_plates:
        
        # Split to just the information of interest
        df = max_combo_df.loc[max_combo_df.Plate == plate].copy()
    
        # generate a holoviews plot
        hm_dict[plate] = make_heatmap(df, title=plate)

    # make logseqdepth column
    max_combo_df['logseqdepth'] = np.log(
        max_combo_df['WellSeqDepth'], 
        out=np.zeros_like(
            max_combo_df['WellSeqDepth'], 
            dtype=float
        ),
        where=max_combo_df['WellSeqDepth'] != 0
    )

    # logseqdepth heatmap
    cmap = list(reversed(cc.CET_D9))

    # Set the center
    center = np.log(10)

    # Adjust if it is greater than max of data (avoids ValueError)
    if max_combo_df['logseqdepth'].max() <= center:

        # Log a warning
        log_warning(f"All wells associated with {title} have a read depth <=10. "
                    "Be careful comparing heatmaps between this plate and others. "
                    "Be careful using this data; sequencing was not good.")

        # Adjust the center
        center = max_combo_df['logseqdepth'].median()

    # set color levels
    color_levels = stretch_color_levels(df['logseqdepth'], center, cmap)

    # Uniform color levels
    for _hm in hm_dict.values():
        _hm.opts({'HeatMap': {'color_levels': color_levels}})
    
    # plot from the dictionary
    hm_holomap = hv.HoloMap(
        hm_dict, 
        kdims=['Plate']
    )

    return hm_holomap

def save_heatmap_to_file(heatmaps, outputdir):
    
    file_path = os.path.join(outputdir, "Platemaps", "Platemaps")
    hv.renderer('bokeh').save(heatmaps, file_path)

def stretch_color_levels(data, center, cmap):
    """Stretch a color map so that its center is at `center`. Taken
    from hw4.2 solutions to 2019 bebi103a, probably with permission. 
    This is best for centering divergent color maps.
    """
    # don't allow a color map with only one color
    if len(cmap) == 1:
        raise RuntimeError("Must have `len(cmap)` > 1.")
    
    # Scale dist
    dist = max(max(data) - center, center - 0)
    dist += dist / 100

    color_levels = list(np.linspace(center-dist, center+dist, len(cmap)+1))

    # Ignore if only one value is present
    if len(np.unique(color_levels)) == 1:
        color_levels = None
    
    return color_levels


def make_heatmap(df, title):
    """Generates a heatmap from ssSeq data using Holoviews with bokeh backend."""

    # Convert SeqDepth to log for easier visualization.
    df['logseqdepth'] = np.log(df['WellSeqDepth'], out=np.zeros_like(df['WellSeqDepth'], dtype=float),
                               where=df['WellSeqDepth'] != 0)

    # Create necessary Row and Column values and sort
    df['Row'] = df.apply(lambda row: row['Well'][0], axis=1)
    df['Column'] = df.apply(lambda row: int(row['Well'][1:]), axis=1)
    df = df.sort_values(['Column', 'Row'])
    df['Column'] = df['Column'].astype('str')

    # Set some base opts
    opts = dict(invert_yaxis=True, title=title, show_legend=True)

    # logseqdepth heatmap
    cmap = list(reversed(cc.CET_D9))

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

    color_levels = stretch_color_levels(df['logseqdepth'], center, cmap)

    # Get heights
    n_rows = len(df['Row'].unique())
    n_cols = len(df['Column'].unique())
    height = int(50 * n_rows)
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
        ['Column', 'Row'],
        ['logseqdepth','VariantCombo','AlignmentFrequency','WellSeqDepth','Well']
    ).opts(
        **opts,
        colorbar=True,
        cmap=cmap,
        height=height,
        width=width,
        line_width=5,
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
    # colors = bokeh.palettes.Plasma5
    colors = ['#337D1F', '#94CD35', '#FFC300', '#C62C20']
    cmap = {bin: color for bin, color in zip(bins, colors)}

    # apply binning function to the AlignmentFrequency
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(
        bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows*1.15

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
        line_width=8,
        fill_alpha=0,
        line_alpha=1,
        legend_position='right',
        size=box_size)

    # residue labels
    def split_variant_labels(mutation_string):
        
        mutation_string = mutation_string.replace('?','')
        new_line_mutations = mutation_string.replace('_','\n')
        
        return new_line_mutations
    
    _df = df.copy()
    _df['Labels'] = _df['VariantCombo'].apply(split_variant_labels)
    
    labels = hv.Labels(
        _df,
        ['Column', 'Row'],
        'Labels',
    ).opts(text_font_size='9pt', **opts)

    # return formatted final plot
    return (hm*boxes*labels).opts(frame_height=550,
                                  frame_width=550 * 3 // 2,
                                  border=50,
                                  show_legend=True)

#### Read quality chart ####
def generate_read_qual_chart(seqpairs, output_dir):
    """Makes histograms of read qualities and saves to designated location."""
    
    # Generate forward and reverse qualities
    all_qualities = [seqpair.read_quals() for seqpair in seqpairs]
    f_qual_counts = np.unique([int(qual[0]) for qual in all_qualities if not np.isnan(qual[0])],
                              return_counts = True)
    r_qual_counts = np.unique([int(qual[1]) for qual in all_qualities if not np.isnan(qual[1])],
                              return_counts = True)
        
    # Plot counts
    p_f = plot_read_qual(f_qual_counts).opts(title='Forward Read Quality')
    p_r = plot_read_qual(r_qual_counts).opts(title='Reverse Read Quality')
    
    # Render to bokeh and combine into single chart
    p = row(
        hv.render(p_f),
        hv.render(p_r)
    )

    # Output as html
    bokeh.io.output_file(os.path.join(output_dir, "Qualities", "QualityPlot.html"))
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
    seq_output_path
    ):

    # If the well column of the data_df is not zero-padded, pad it
    data_df['Well'] = data_df['Well'].apply(ns.parsers.pad)

    # Read in the results of deSeq
    seq_df = pd.read_csv(seq_output_path+'AminoAcids_Decoupled_Max.csv')

    # Find wells that have multiple mutations since this only works for
    # single-mutant libraries
    count_df = pd.DataFrame(seq_df.groupby(
        ['Plate', 'Well'])['AaPosition'].count())
    mut_count_dict = count_df.to_dict()['AaPosition']

    # Get the number of mutations in a sequence for later filtering
    seq_df['MutCount'] = seq_df.apply(
        lambda x: mut_count_dict[x['Plate'], x['Well']], 
        axis=1)

    # Rename columns for downstream use
    seq_df = seq_df.rename(
        columns={'AaPosition': 'Position', 'Aa': 'AA'}
    )

    # Remove dead sequencing wells
    seq_df = seq_df.loc[(seq_df['Flag'] != '#DEAD#')].copy()

    # Remove sequencing wells with a flag
    seq_df = seq_df.loc[seq_df['Flag'].isna()].copy()

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
    return data_df.merge(seq_df, how='outer')

def activity_plot(
    df, 
    plate, 
    parent, 
    missing, 
    value, 
    color_levels,
    activity_range,
    title
):

    temp_df = df.copy()

    # list of AAs
    AAs = list('ACDEFGHIKLMNPQRSTVWY')

    for AA in missing:

        # Create DataFrame for these values
        temp_df = temp_df.append(
            pd.DataFrame(
                [[None] * (len(temp_df.columns))],
                columns = temp_df.columns
            ),
            ignore_index=True
        )
        temp_df.at[len(temp_df) - 1, 'AA'] = AA
    
    temp_df = temp_df.sort_values('Sort')
    
    # Placeholder plotting
    opts = dict(
        width=600,
        height=500,
        xrotation=45
    )

    # Create tooltips for plotting
    tooltips = [
        ('Well', '@Well'),
        ('Amino Acid', '@AA'),
        ('Seq. depth', '@WellSeqDepth'),
        ('Align. freq.', '@AlignmentFrequency'),
        (value, '@{'+value+'}')
    ]
    hover = HoverTool(tooltips=tooltips)

    # Plot the data for which the AA is known
    _df = temp_df[temp_df['AA'] != 'Unknown'].copy()
    
    p = ns.viz.plot_bar(
        _df,
        'AA',
        value_name=value,
        color=value,
        cmap='coolwarm'
    ).opts(
        {'Scatter': dict(tools=[hover])}
    )

    # Plot the data w/ an unknown AA
    _df = temp_df[temp_df['AA'] == 'Unknown'].copy()
    if len(_df) != 0:
        unknowns = True
    else:
        unknowns = False
        _df = pd.DataFrame({
            'AA': ['Unknown'],
            value: [np.nan],
            'Sort': 23
        })
    
    if unknowns:
        p_unknown = ns.viz.plot_bar(
            _df,
            'AA',
            value_name=value,
            color=value,
            cmap='coolwarm'
        ).opts({
            'Bars': dict(line_alpha=0, fill_alpha=0),
            'Scatter': dict(tools=[hover])
        })

        p = hv.Overlay([p, p_unknown])

    # Center colormap
    if color_levels is not None:
        p = p.opts({'Bars': {'color_levels': color_levels}})

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
    p.opts(xlabel='Residue')
    if title is not None:
        if not isinstance(title, str):
            raise ValueError(
                f'Title must be of type str, you passed {type(title)}.')
        p.opts(title=title)
    if activity_range is not None:
        p.opts(ylim=activity_range)

    return p

def count_plot(
    temp_df,
    missing,
    sort,
    hist_range
):
    
    # Get observed AA counts
    counts = temp_df['AA'].value_counts()

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
    df_counts = df_counts.sort_values('Sort')

    # Make chart
    p_counts = hv.Bars(
        df_counts,
        'Residue',
        'Counts'
    ).opts(
        width=600,
        height=200,
        xrotation=45
    )

    if hist_range is not None:
        p_counts = p_counts.opts(ylim=hist_range)
    
    return p_counts

def plot_variant_activities(
    seq_func_data,
    value,
    min_align_freq=0.8,
    min_seq_depth=10,
    hist=True,
    title=None,
    color_levels=None,
    activity_range=None,
    hist_range=None
):
    """
    Takes a DataFrame of single mutant sequence/fitness data and plots the
    activities of all sequences as well as a frequency histogram of each 
    mutant.
    """
    
    df = seq_func_data.copy()

    # list of AAs
    AAs = list('ACDEFGHIKLMNPQRSTVWY')

    ##### Activity plot

    ## TO-DO: Handle a dataframe that doesn't have positive/negative controls ##

    # Find Parent
    parent = f"Parent ({df['AA'][df['Type'] == 'Parent'].unique()[0]})"

    # Set unknowns vs controls
    def set_knowns(df):
        
        if df['Type'] == 'Negative':
            df['AA'] = 'Negative'
        elif df['Type'] == 'Sterile':
            df['AA'] = 'Sterile'
        elif df['Type'] == 'Parent':
            df['AA'] = parent
        elif pd.isna(df['AA']):
            df['AA'] = 'Unknown'
        elif df['AlignmentFrequency'] < min_align_freq:
            df['AA'] = 'Unknown'
        elif df['WellSeqDepth'] < min_seq_depth:
            df['AA'] = 'Unknown'
        elif df['MutCount'] > 1:
            df['AA'] = 'Unknown'

        return df

    working_df = df.copy()
    working_df = working_df.apply(set_knowns, axis=1)
    
    activity_plot_dict = {}
    histogram_plot_dict = {}

    for plate in working_df['Plate'].unique():
        temp_df = working_df.loc[working_df['Plate'] == plate].copy()
        
        # Sort meaningfully
        sort = {
            parent: 0,
            # Known variants : 1-19, done in update line
            'Negative': 21,
            'Sterile': 22,
            'Unknown': 23,
        }

        sort.update({AA: i for i, AA in enumerate(AAs, 1)})

        temp_df['Sort'] = temp_df['AA'].replace(sort)

        # Find missing AAs
        missing = tuple(set(AAs) - set(temp_df['AA'].unique()))
        
        #### Activity Plot ####
        p = activity_plot(
            temp_df, 
            plate, 
            parent, 
            missing, 
            value,
            color_levels=color_levels,
            activity_range=activity_range,
            title=title
        )
        
        ##### Histogram #####
        if hist:

            p_counts = count_plot(
                temp_df,
                missing,
                sort,
                hist_range=hist_range
            )

            # Since we are plotting the histogram, remove xlabel from activity plot
            p.opts(xlabel='')
            
            histogram_plot_dict[plate] = p_counts
        
        activity_plot_dict[plate] = p

    activity_hmap = hv.HoloMap(
        activity_plot_dict,
        kdims='Plate'
    ).opts(
        {'Bars': {'framewise': True}}
    )

    if hist:
        
        count_hmap = hv.HoloMap(
            histogram_plot_dict,
            kdims='Plate'
        )

        return (activity_hmap+count_hmap).cols(1)

    else:
        # print('No histogram')
        return activity_hmap
