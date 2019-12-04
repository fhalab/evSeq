import numpy as np
import pandas as pd

import holoviews as hv
import colorcet as cc
import bokeh.io
from bokeh.layouts import row

hv.extension('bokeh')
hv.renderer('bokeh')


#### Heatmap ####
def generate_sequencing_heatmap(df, plate, hm_output_file):
    """Saves a heatmap html generated from from ssSeq data."""
    hm = make_heatmap(df, title=plate)
    
    hv.renderer('bokeh').save(hm, hm_output_file)

def stretch_color_levels(data, center, cmap):
    """Stretch a color map so that its center is at `center`. Taken
    from hw4.2 solutions to 2019 bebi103a, probably with permission.
    """
    if len(cmap) == 1:
        raise RuntimeError("Must have `len(cmap)` > 1.")
        
    if not 0 < center < max(data):
        raise ValueError('Must have min(data) < center < max(data).')
    
    dist = max(max(data) - center, center - 0)
    dist += dist / 100
    
    return list(np.linspace(center-dist, center+dist, len(cmap)+1))

def make_heatmap(df, title):
    """Generates a heatmap from ssSeq data using Holoviews with bokeh backend."""
    
    # Create necessary Row and Column values and sort
    df['Row'] = df.apply(lambda row: row['Well'][0], axis=1)
    df['Column'] = df.apply(lambda row: int(row['Well'][1:]), axis=1)
    df = df.sort_values(['Column','Row'])
    
    # Set some base opts
    opts = dict(invert_yaxis=True,title=title,show_legend=True)
    
    # Convert SeqDepth to log for easier visualization.
    df['logseqdepth'] = np.log(df['WellSeqDepth'])
    
    # logseqdepth heatmap
    cmap = list(reversed(cc.CET_D9))

    hm = hv.HeatMap(
        df,
        ['Column', 'Row'],
        'logseqdepth'
        ).opts(
            **opts,
            colorbar=True,
            cmap=cmap,
            clipping_colors={'NaN': '#DCDCDC'},
            color_levels=stretch_color_levels(df['logseqdepth'], np.log(10), cmap),
            colorbar_opts=dict(
                title='LogSeqDepth',
                background_fill_alpha=0
            )
        )
    
    # Patrick did this, don't @ me
    # @Kadina 
    def bin_align_freq(value):
        if value > 0.99:
            bin_vals = '0.99+'
        if value <= 0.99 and value > 0.98:
            bin_vals = '0.98-0.99'
        if value <= 0.98 and value > 0.95:
            bin_vals = '0.95-0.98'
        if value <= 0.95 and value > 0.9:
            bin_vals = '0.90-0.95'
        if value <= 0.9:
            bin_vals = '<0.90'

        return bin_vals
    
    # Bin alignment frequencies for easier viz
    bins = ['0.99+','0.98-0.99','0.95-0.98','0.90-0.95','<0.90']
    colors = bokeh.palettes.Plasma5
    cmap = {bin: color for bin, color in zip(bins, colors)}
    
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(bin_align_freq)
    
    # alignment frequency heamap
    boxes = hv.Points(
        df.sort_values(['AlignmentFrequency'], ascending=False),
        ['Column', 'Row'],
        'AlignmentFrequencyBinned'
    ).opts(**opts)
    
    # residue labels
    labels = hv.Labels(
        df,
        ['Column', 'Row'],
        'VariantCombo'
    ).opts(**opts)
    
    return (hm*boxes*labels).opts(frame_height=550, 
                                  frame_width=550 * 3 // 2, 
                                  border=50,
                                  show_legend=True)

#### Read quality chart ####
def generate_read_qual_chart(counts, path):
    """Makes histograms of read qualities and saves to designated location."""
    # Unpack
    f_qual_counts, r_qual_counts = counts

    p_f = plot_read_qual(f_qual_counts).opts(title='Forward Read Quality')
    p_r = plot_read_qual(r_qual_counts).opts(title='Reverse Read Quality')
    
    # Render to bokeh and combine into single chart
    p = row(
        hv.render(p_f),
        hv.render(p_r)
    )

    # Output as html
    try:
        bokeh.io.output_file(path)
        bokeh.io.save(p)
    
    # Save can fail in many ways, not worried about naked exception
    except:
        pass

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