# Import required modules
import numpy as np
import pandas as pd
import os

import holoviews as hv
import colorcet as cc
import bokeh.io
from bokeh.layouts import row

from .Logging import log_warning

hv.extension('bokeh')
hv.renderer('bokeh')


#### Heatmap ####
def generate_sequencing_heatmap(max_combo_df, output_dir):
    """Saves a heatmap html generated from from ssSeq data."""
    
    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()
    
    # Generate plots for each plate
    for plate in unique_plates:
        
        # Split to just the information of interest
        df = max_combo_df.loc[max_combo_df.Plate == plate].copy()
    
        # generate a holoviews plot
        hm = make_heatmap(df, title=plate)
        
        # render the plot using bokeh and save to html file
        hv.renderer('bokeh').save(hm, os.path.join(output_dir, "Platemaps",
                                                   f"{plate}"))

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


def generate_sequencing_heatmap(max_combo_df):
    """Saves a heatmap html generated from from ssSeq data."""

    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()
    hms = []

    # Generate plots for each plate
    for plate in unique_plates:

        # Split to just the information of interest
        df = max_combo_df.loc[max_combo_df.Plate == plate].copy()

        # generate a holoviews plot
        hm = make_heatmap(df, title=plate)
        hms.append(hm)

    layout = hv.Layout(hms).cols(2)
    return layout


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

    # generate the heatmap
    hm = hv.HeatMap(
        df,
        ['Column', 'Row'],
        'logseqdepth'
    ).opts(
        **opts,
        colorbar=True,
        cmap=cmap,
        height=height,
        width=width,
        clipping_colors={'NaN': '#DCDCDC'},
        color_levels=color_levels,
        colorbar_opts=dict(
            title='LogSeqDepth',
            background_fill_alpha=0
        )
    )

    # function to bin the alignment frequencies into more relevant groupings
    def bin_align_freq(value):
        if value > 0.99:
            bin_vals = '0.99+'
        if value <= 0.99 and value > 0.98:
            bin_vals = '0.98-0.99'
        if value <= 0.98 and value > 0.95:
            bin_vals = '0.95-0.98'
        if value <= 0.95 and value > 0.9:
            bin_vals = '0.90-0.95'

        # anything below 0.9 should really be discarded
        if value <= 0.9:
            bin_vals = '<0.90'

        return bin_vals

    # Bin alignment frequencies for easier viz
    bins = ['0.99+', '0.98-0.99', '0.95-0.98', '0.90-0.95', '<0.90']
    # colors = bokeh.palettes.Plasma5
    colors = ['#337D1F', '#94CD35', '#FFC300', '#FF5733', '#C62C20']
    cmap = {bin: color for bin, color in zip(bins, colors)}

    # apply binning function to the AlignmentFrequency
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(
        bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows*1.21

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
    labels = hv.Labels(
        df,
        ['Column', 'Row'],
        'SimpleCombo',
    ).opts(text_font_size='10pt', **opts)

    # return formatted final plot
    return (hm*boxes*labels).opts(frame_height=550,
                                  frame_width=550 * 3 // 2,
                                  border=50,
                                  show_legend=True)

    # function to bin the alignment frequencies into more relevant groupings
    def bin_align_freq(value):
        if value > 0.99:
            bin_vals = '0.99+'
        if value <= 0.99 and value > 0.98:
            bin_vals = '0.98-0.99'
        if value <= 0.98 and value > 0.95:
            bin_vals = '0.95-0.98'
        if value <= 0.95 and value > 0.9:
            bin_vals = '0.90-0.95'

        # anything below 0.9 should really be discarded
        if value <= 0.9:
            bin_vals = '<0.90'

        return bin_vals

    # Bin alignment frequencies for easier viz
    bins = ['0.99+', '0.98-0.99', '0.95-0.98', '0.90-0.95', '<0.90']
    # colors = bokeh.palettes.Plasma5
    colors = ['#337D1F', '#94CD35', '#FFC300', '#FF5733', '#C62C20']
    cmap = {bin: color for bin, color in zip(bins, colors)}

    # apply binning function to the AlignmentFrequency
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(
        bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows*1.21

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
    labels = hv.Labels(
        df,
        ['Column', 'Row'],
        'SimpleCombo'
    ).opts(text_font_size='10pt', **opts)

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
