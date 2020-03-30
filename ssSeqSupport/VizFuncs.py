import numpy as np
import pandas as pd

import holoviews as hv
import colorcet as cc
import bokeh.io
from bokeh.layouts import row

from . import LogWarning

hv.extension('bokeh')
hv.renderer('bokeh')


#### Heatmap ####
def GenerateSequencingHeatmap(df, plate, hm_output_file):
    """Saves a heatmap html generated from from ssSeq data."""
    
    # generate a holoviews plot
    hm = MakeHeatmap(df, title=plate)
    
    # render the plot using bokeh and save to html file
    hv.renderer('bokeh').save(hm, hm_output_file)

def StretchColorLevels(data, center, cmap):
    """Stretch a color map so that its center is at `center`. Taken
    from hw4.2 solutions to 2019 bebi103a, probably with permission. 
    This is best for centering divergent color maps.
    """
    # don't allow a color map with only one color
    if len(cmap) == 1:
        raise RuntimeError("Must have `len(cmap)` > 1.")
        
    # check that the center passed is within the data
    if not 0 < center < max(data):
        np.save("./TestData.npy", data.values)
        raise ValueError('Must have min(data) < center < max(data).')
    
    dist = max(max(data) - center, center - 0)
    dist += dist / 100
    
    return list(np.linspace(center-dist, center+dist, len(cmap)+1))

def MakeHeatmap(df, title):
    """Generates a heatmap from ssSeq data using Holoviews with bokeh backend."""
    
    # Add 0s to wells that have no data
    well_list = [row+column for row in ['A','B','C','D','E','F','G','H'] for column in ['01','02','03','04','05','06','07','08','09','10','11','12']]

    for well in well_list:
    
    if well not in df['Well'].unique():
        
        temp_df = pd.DataFrame([[np.nan,well,np.nan,np.nan,'    ',0,0,np.nan,np.nan]],columns=df.columns)
        
        df = pd.concat([df,temp_df],sort=False)

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
    
    # Set the center
    center = np.log(10)
    
    # Adjust if it is greater than max of data (avoids ValueError)
    if df['logseqdepth'].max() <= center:
        
        # Log a warning
        LogWarning(f"All wells associated with {df.Plate.values[0]} have a read depth <=10."
                   "Be careful comparing heatmaps between this plate and others."
                   "Be careful using this data; sequencing was not good.")
        
        # Adjust the center
        center = df['logseqdepth'].median()

    # generate the heatmap
    hm = hv.HeatMap(
        df,
        ['Column', 'Row'],
        'logseqdepth'
        ).opts(
            **opts,
            colorbar=True,
            cmap=cmap,
            xmarks=100,
            ymarks=100,
            clipping_colors={'NaN': '#DCDCDC'},
            color_levels=StretchColorLevels(df['logseqdepth'], center, cmap),
            colorbar_opts=dict(
                title='LogSeqDepth',
                background_fill_alpha=0
            )
        )
    
    # Patrick did this, don't @ me
    # @Kadina 
    # @@Patrick
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
    bins = ['0.99+','0.98-0.99','0.95-0.98','0.90-0.95','<0.90']
    # colors = bokeh.palettes.Plasma5
    colors = ['#337D1F','#94CD35','#FFC300','#FF5733','#C62C20']
    cmap = {bin: color for bin, color in zip(bins, colors)}
    
    # apply binning function to the AlignmentFrequency
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(bin_align_freq)
    
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
        legend_position = 'right',
        size=62)
    
    # residue labels
    labels = hv.Labels(
        df,
        ['Column', 'Row'],
        'VariantCombo'
    ).opts(**opts)
    
    # return formatted final plot
    return (hm*boxes*labels).opts(frame_height=550, 
                                  frame_width=550 * 3 // 2, 
                                  border=50,
                                  show_legend=True)

#### Read quality chart ####
def GenerateReadQualChart(counts, path):
    """Makes histograms of read qualities and saves to designated location."""
    # Unpack
    f_qual_counts, r_qual_counts = counts

    p_f = PlotReadQual(f_qual_counts).opts(title='Forward Read Quality')
    p_r = PlotReadQual(r_qual_counts).opts(title='Reverse Read Quality')
    
    # Render to bokeh and combine into single chart
    p = row(
        hv.render(p_f),
        hv.render(p_r)
    )

    # Output as html
    bokeh.io.output_file(path)
    bokeh.io.save(p)

def PlotReadQual(counts):
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