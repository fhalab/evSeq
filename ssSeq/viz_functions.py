import numpy as np
import pandas as pd

import holoviews as hv
import colorcet as cc
import bokeh.io
hv.extension('bokeh')

def generate_sequencing_heatmap(df, plate, hm_output_file):
    
    plot = make_plot(df,title=plate)
    
    hv.renderer('bokeh').save(plot,hm_output_file)

# function taken from 2019 bebi103a hw4.2
def stretch_color_levels(data, center, cmap):
    """Stretch a color map so that its center is at `center`."""
    if len(cmap) == 1:
        raise RuntimeError("Must have `len(cmap)` > 1.")
        
    if not 0 < center < max(data):
        raise ValueError('Must have min(data) < center < max(data).')
    
    dist = max(max(data) - center, center - 0)
    dist += dist / 100
    
    return list(np.linspace(center-dist, center+dist, len(cmap)+1))

def make_plot(df,title):
    
    df['Row'] = df.apply(lambda row: row['Well'][0], axis=1)
    df['Column'] = df.apply(lambda row: int(row['Well'][1:]),axis=1)
    df = df.sort_values(['Column','Row'])
    
    opts = dict(invert_yaxis=True,title=title,show_legend=True)
    
    df['logseqdepth'] = np.log(df['WellSeqDepth'])
    
    cmap = list(reversed(cc.CET_D9))
    
    hm = hv.HeatMap(df,['Column','Row'],'logseqdepth').opts(**opts,
                                                            colorbar=True,
                                                            cmap=cmap,
                                                            clipping_colors={'NaN':'#dcdcdc'},
                                                            color_levels=stretch_color_levels(df['logseqdepth'], np.log(10), cmap),
                                                            colorbar_opts=dict(title='LogSeqDepth',
                                                                               background_fill_alpha=0))
    
    
    # Patrick did this, don't @ me    
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
    
    bins = ['0.99+','0.98-0.99','0.95-0.98','0.90-0.95','<0.90']
    colors = bokeh.palettes.Plasma5
    cmap = {bin: color for bin, color in zip(bins, colors)}
    
    df['AlignmentFrequencyBinned'] = df['AlignmentFrequency'].apply(bin_align_freq)
    
    boxes = hv.Points(df.sort_values(['AlignmentFrequency'],ascending=False),['Column','Row'],'AlignmentFrequencyBinned').opts(**opts,
                                                                                                                               marker='square',
                                                                                                                               line_color='AlignmentFrequencyBinned',
                                                                                                                               cmap=cmap,
                                                                                                                               line_width=7,
                                                                                                                               fill_alpha=0,
                                                                                                                               line_alpha=1,
                                                                                                                               legend_position = 'right',
                                                                                                                               size=60)
    
    labels = hv.Labels(df,['Column','Row'],'VariantCombo').opts(**opts)
    
    return (hm*boxes*labels).opts(frame_height=550, 
                                  frame_width=550 * 3 // 2, 
                                  border=50,
                                  show_legend=True)