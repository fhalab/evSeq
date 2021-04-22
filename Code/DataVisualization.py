# Import required modules
import numpy as np
import pandas as pd
import os

import holoviews as hv
import colorcet as cc
import bokeh.io
from bokeh.layouts import row
from bokeh.models import HoverTool

from .Logging import log_warning

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
def plot_variant_activities(
    data_df,
    seq_output_path,
    min_align_freq=0.8,
    min_seq_depth=10,
    value=None,
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
    def get_seq_data(filepath):
        seq_df = pd.read_csv(filepath+'AminoAcids_Decoupled_Max.csv')

        count_df = pd.DataFrame(seq_df.groupby(['Plate','Well'])['AaPosition'].count())
        mut_count_dict = count_df.to_dict()['AaPosition']

        def get_mutation_counts(row):
            return mut_count_dict[row['Plate'], row['Well']]
        
        seq_df['MutCount'] = seq_df.apply(get_mutation_counts, axis=1)

        seq_df = seq_df.rename(
            columns={'AaPosition': 'Position', 'Aa': 'AA'}
            )
        seq_df = seq_df.loc[(seq_df['Flag'] != '#DEAD#')].copy()
        seq_df = seq_df.astype({'Position': 'int64'})
        seq_df = seq_df.loc[seq_df['Flag'].isna()].copy()
        return seq_df

    seq_df = get_seq_data(seq_output_path)

    df = data_df.merge(seq_df, how='outer')

    # list of AAs
    AAs = list('ACDEFGHIKLMNPQRSTVWY')

    ##### Activity plot
    # Determine activity value
    if value is None:
        raise ValueError("You have not specified a column with the data value.")

    # Find missing AAs
    missing = set(AAs) - set(df['AA'].unique())

    # Find Parent
    parent = f"Parent ({df['AA'][df['Type'] == 'Parent'].unique()[0]})"

    # Set unknowns vs controls
    def set_knowns(working_df):
        if working_df['Type'] == 'Negative':
            working_df['AA'] = 'Negative'
        elif working_df['Type'] == 'Sterile':
            working_df['AA'] = 'Sterile'
        elif working_df['Type'] == 'Parent':
            working_df['AA'] = parent
        elif pd.isna(working_df['AA']):
            working_df['AA'] = 'Unknown'
        elif working_df['AlignmentFrequency'] < min_align_freq:
            working_df['AA'] = 'Unknown'
        elif working_df['WellSeqDepth'] < min_seq_depth:
            working_df['AA'] = 'Unknown'
        return working_df
    
    working_df = df.copy()
    working_df = working_df.apply(set_knowns, axis=1)

    seq_func_plot_dict = {}

    for plate in working_df['Plate'].unique():
        temp_df = working_df.loc[working_df['Plate'] == plate].copy()

        # Subset the data (might want to keep more for encodings?)
        encodings = ['AA', value, 'AlignmentFrequency', 'Well', 'WellSeqDepth']
        temp_df = temp_df[encodings]

        # Add in values for missing
        for AA in missing:
            temp_df.loc[len(temp_df)] = [AA, *[np.nan]*(len(encodings)-1)]
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
        temp_df.sort_values('Sort', inplace=True)
        # Placeholder plotting
        opts = dict(
            variable='AA',
            value=value,
            color=value,
            sort='Sort',
            cmap='coolwarm',
            show_points=True,
            width=600,
            height=500,
            xrotation=45
        )
        
        _df = temp_df[temp_df['AA'] != 'Unknown']

        p = Viz(_df).plot_bar(**opts)
        _df = temp_df[temp_df['AA'] == 'Unknown']
        if len(_df) != 0:
            unknowns = True
        else:
            unknowns = False
            _df = pd.DataFrame({
                'AA': ['Unknown'],
                value: [np.nan],
                'Sort': 23
            })
        p_unknown = Viz(_df).plot_bar(**opts).options(
            {'Bars': dict(line_alpha=0, fill_alpha=0)}
        )
        p = hv.Overlay([p, p_unknown])
        # Center colormap
        if color_levels is not None:
            p = p.options({'Bars': {'color_levels': color_levels}})
        # Place 'n.d.' for missing amino acids
        if activity_range is None:
            height = temp_df[value].max() / 20
        else:
            height = activity_range[1] / 20
        for AA in missing:
            p = p*hv.Text(AA, height, 'n.d.', rotation=90)
        if not unknowns:
            p = p*hv.Text('Unknown', height*1.25, 'none', rotation=90)

        ##### Histogram
        # Get observed AA counts
        counts = temp_df['AA'].value_counts()
        # Create DataFrame of counts, accounting for missing values
        df_counts = pd.DataFrame({
            'Residue': list(counts.index),
            'Counts' : [value if AA not in missing else 0
                        for value, AA in zip(counts.values, counts.index)],
        })
        # Add in nans
        nan_count = sum(temp_df['AA'].isna())
        df_counts.loc[len(temp_df)] = ['Unknown', nan_count]
        # Sort meaningfully, again...
        df_counts['Sort'] = df_counts['Residue'].replace(sort)
        df_counts.sort_values('Sort', inplace=True)
        # Make chart
        p_counts = hv.Bars(df_counts, 'Residue', 'Counts').opts(
            width=600, height=200, xrotation=45)
        # Relabel from AA
        p.opts(xlabel='Residue')
        if title is not None:
            if not isinstance(title, str):
                raise ValueError(
                    f'Title must be of type str, you passed {type(title)}.')
            p.opts(title=title)
        if activity_range is not None:
            p.opts(ylim=activity_range)
        if hist is True:
            p.opts(xlabel='')
            if hist_range is not None:
                p_counts = p_counts.opts(ylim=hist_range)
            p = (p+p_counts).cols(1)

        seq_func_plot_dict[plate] = p
    return hv.HoloMap(seq_func_plot_dict, kdims='Plate').collate()

# checking functions


def check_df_col(df, column, name=None):
    """Checks for the presence of a column (or columns) in a tidy
    DataFrame with an informative error message. Passes silently,
    otherwise raises error.
    """
    if column is not None:

        if not isinstance(column, (list, tuple)):
            column = [column]

        for col in column:
            if name is None:
                error_message = f"The value '{col}' is not present in any of "\
                    "the columns of your DataFrame."
            else:
                error_message = f"Your {name} value '{col}' is not present in"\
                    " any of the columns of your DataFrame."
            error_message += "\nYou may be looking for:\n  " + \
                str(list(df.columns))

            if col not in df.columns:
                raise ValueError(error_message)


def check_controls(set_type, set_controls):
    """Checks the old control nomenclature (set_type, set_controls) and
    normalizes them to provide a better infrastructure for ScreenData
    class to assign controls.
    """
    if set_type is not None:
        control_dict = set_type
        standardized = False
    if set_controls is not None:
        control_dict = set_controls
        standardized = True
    else:
        control_dict = None
        standardized = None

    controls = control_dict, standardized

    return controls


def check_replicates(df, variable, value, grouping):
    """Checks for the presence of replicates in the values of a dataset,
    given some experimental conditions. Returns True if the standard
    deviation of the values of each group (if more than one exists) is
    greater than, indicating that replicates were performed under the
    given criteria.
    Parameters
    ----------
    df : Pandas DataFrame in tidy format
        The data set to be checked for replicates
    variable : immutable object
        Name of column of data frame for the independent variable,
        indicating a specific experimental condition.
    value : immutable object
        Name of column of data frame for the dependent variable,
        indicating an experimental observation.
    group : immutable object of list of immutable objects
        Column name or list of column names that indicates how the
        data set should be split.
    
    Returns
    -------
    replicates : boolean
        True if replicates are present.
    df_out : the DataFrame containing averaged 'variable' values, if
        replicates is True. Otherwise returns the original DataFrame.
    """
    # Unpack the experimental conditions into a single list of arguments
    if not isinstance(grouping, (list, tuple)):
        grouping = [grouping]
    args = [elem for elem in [variable, *grouping] if elem is not None]

    # Get stdev of argument groups
    grouped = df.groupby(args)[value]
    group_stdevs = grouped.std().reset_index()
    group_stdev = group_stdevs[value].mean()

    # Determine if there are replicates (mean > 0)
    replicates = bool(group_stdev > 0)

    # Average the values and return
    df_mean = grouped.mean().reset_index()
    df_mean.columns = list(df_mean.columns[:-1]) + ['Mean of ' + str(value)]
    df_return = df.merge(df_mean)

    return replicates, df_return


class Viz():

    def __init__(
        self,
        cls,
        value=None,
        title=None,
        user='pja'
    ):
        """Generates common visualizations from data via different
        methods.
        Parameters:
        -----------
        cls : DataFrame or arnoldLab_utils object (Tecan, LCMS, etc.)
        value : str
            Important value in the DataFrame (i.e., not a location or
            annotation like Well or Type). If None, assumed to be the
            last (-1) column in the DataFrame.
        title : str
            Title to use on the plot (experiment, etc.).
        user : str
            The user; adjusts styling defaults. Default 'pja'.
        Returns:
        --------
        Viz object with methods for visualization.
        """
        # Passing in a DataFrame directly
        if isinstance(cls, pd.DataFrame):
            self.data = cls

        # Inheriting DataFrame from class.data attribute
        else:
            self.data = cls.data

        # Grab some default info
        self.cls = cls
        self._value = value
        self._title = title
        self._user = user
        self._user_support = ['pja', 'ejw']
        self.types = self._get_types()

        # Plotting settings
        self.variant_cmap, self.default_cmap = self.set_defaults()
        self._update_variant_cmap()

    @property
    def value(self):
        if self._value is not None:
            check_df_col(self.data, self._value, name="value")
        else:
            self._value = self.data.columns[-1]

        return self._value

    @value.setter
    def value(self, value):
        self._value = value

    @property
    def title(self):
        if self._title is not None:
            pass
        elif self.data.columns.name is not None:
            self._title = self.data.columns.name
        else:
            self._title = datetime.datetime.now().strftime('%y%m%d %H:%M%p analysis')

        return self._title

    @property
    def user(self):
        return self._user

    def set_defaults(self):
        if self._user not in self._user_support:
            raise NotImplementedError(
                f'User not specified correctly, update or choose one of: {self.user_support}'
            )

        variant_cmap, default_cmap = hv_defaults(self._user).defaults()

        return variant_cmap, default_cmap

    def _update_variant_cmap(self):
        """Adjusts color map if parent plate."""
        # For parent plate case
        if 'Variant' not in self.types:
            cmap = self.variant_cmap
            cmap.update({'Parent': cmap['Variant']})
        else:
            cmap = hv_defaults(self._user).variant_cmap

        return cmap

    def _get_types(self):
        """Adjust types for parent plate."""
        types = []
        if 'Type' in self.data.columns:
            types = self.data['Type'].unique()

            # For parent plates
            if ('Parent' not in types) and (types is not [None]):
                replace = {
                    'Variant': 'Parent',
                    'Negative': 'Negative',
                    'Sterile': 'Sterile'
                }
                self.data['Type'].replace(replace, inplace=True)
                types = self.data['Type'].unique()

        return types

# Bar plot (with scatter) #########################

    def plot_bar(
        self,
        variable,
        value,
        split=None,
        color=None,
        scheme='viridis',
        sort=None,
        cmap='default',
        show_points='default',
        show_all=True,
        legend=False,
        xrotation=0,
        height=350,
        width=300,
        additional_opts=None
    ):
        """Converts a tidy DataFrame a bar plot, taking care to show all
        the data. Bars are given as the average of each grouping of the
        variable (which can be further split by the split argument), and
        the actual data points are overlaid on top.
        Parameters:
        -----------
        variable : str
            Column in DataFrame representing the variable, plotted on
            the x-axis.
        value : str
            Column in DataFrame representing the quantitative value,
            plotted on the y-axis
        split : str
            The names of one or more columns that further specify the
            way the data is grouped. Defaults to None.
        sort : str
            Which column is used to determine the sorting of the data.
            Defaults to None, and will sort by the condition column
            (alphabetical) if present, otherwise variable.
        cmap : The colormap to use. Any Holoviews/Bokeh colormap is fine.
            Uses Holoviews default if None.
        show_all : bool
            If split is not None, whether or not to use a drop-down or
            to show all the plots (layout). Note that this can be pretty
            buggy from Holoview's layout system. There is usually a way
            to how all the info you want, in a nice way. Just play
            around.
        show_points : bool
            Shows all the data points. I don't even know why this is an
            argument. Default will show points if there are multiple
            replicates. Unless you have a really good reason, don't
            change this.
        legend : str
            First controls whether or not the legend is shown, then its
             position. Defaults to False, though 'top' would be a good
             option, or 'top_left' if using split.
        height : int
            The height of the chart.
        width : int
            The width of the chart.
        additional_opts : dictionary
            A dictionary to pass additional Holoviews options to the
            chart. Flexible; will try all options and only use the
            ones that did not raise an exception. Not verbose.
        
        Returns:
        --------
        chart : the final Holoviews chart
        """
        # Check columns
        check_df_col(self.data, variable, name="variable")
        check_df_col(self.data, value, name="value")
        check_df_col(self.data, split, name="split")
        check_df_col(self.data, sort, name="sort")

        if not isinstance(split, (list, tuple)):
            split = [split]
        replicates, data = check_replicates(self.data, variable, value, split)

        # Sort
        if sort is not None:
            check_df_col(data, sort, name="sort")
            data = data.sort_values(by=sort).reset_index(drop=True)

        # Encode color
        if color is None:
            color = variable

        # Decide colormap
        if cmap == 'default':
            number = len(data[color].unique())
            try:
                cmap = getattr(bokeh.palettes, scheme)(number+3)[1:-1]
            except:
                cmap = getattr(bokeh.palettes, 'viridis')(number+3)[1:-1]

        # Pull out available encodings (column names)
        encodings = [*list(data.columns)]

        # Set options (this is probably horribly inefficient right now)
        base_opts = dict(
            height=height,
            width=width,
            ylim=(0, 1.1*np.max(data[value])),
            xrotation=xrotation,
            color=color,
            cmap=cmap,
            show_legend=legend,
        )

        bar_opts = base_opts
        scat_opts = dict(size=6, fill_alpha=0.2, color='black')

        if additional_opts:
            if not isinstance(additional_opts, dict):
                raise ValueError(
                    'Additional options must be passed in as dict.'
                )

            if 'bar_opts' in additional_opts.keys():
                add_bar_opts = additional_opts.pop('bar_opts')
            else:
                bar_opts = None

            if 'scat_opts' in additional_opts.keys():
                add_scat_opts = additional_opts.pop('scat_opts')
            else:
                scat_opts = None

            base_opts.update(additional_opts)

            bar_opts.update(base_opts)
            bar_opts.update(add_bar_opts)

            scat_opts.update(add_scat_opts)

        # Make bar chart
        bars = hv.Bars(
            data,
            variable,
            [('Mean of ' + str(value), value), *encodings],
        ).opts(**bar_opts)

        # Determine single-point entries
        args = [elem for elem in [variable] + split if elem is not None]
        counts = (data.groupby(args).count() ==
                  1).reset_index()[[variable, value]]
        counts.columns = [variable, 'Counts']

        # Get list of singlets to drop from plotting df
        singlets = counts[counts['Counts']][variable].to_list()

        # Make scatter chart
        points = hv.Scatter(
            data[~data[variable].isin(singlets)],
            variable,
            [value, *encodings],
        ).opts(**scat_opts)

        # Make the split
        if split != [None]:
            bars = bars.groupby(split).opts(**bar_opts)
            points = points.groupby(split).opts(**scat_opts)

            # If split, show as side-by-side, or dropdown
            if show_all is True:
                bars = bars.layout()
                points = points.layout()

        # Output chart as only bars, or bars and points
        if show_points == 'default':
            if replicates:
                chart = bars * points
            else:
                chart = bars
        elif show_points:
            chart = bars * points
        else:
            chart = bars

        return chart

################################
# Defaults! Add your own!
################################


class hv_defaults():

    def __init__(self, user):

        self._user = user

        if self._user == 'pja':

            self.variant_cmap = {
                'Variant': '#3E94FA',
                'Parent': '#8DED81',
                'Negative': '#E18409',
                'Sterile': '#DEDEDE',
            }

            self.default_cmap = [
                "#4c78a8",
                "#f58518",
                "#e45756",
                "#72b7b2",
                "#54a24b",
                "#eeca3b",
                "#b279a2",
                "#ff9da6",
                "#9d755d",
                "#bab0ac",
            ]

        if self._user == 'ejw':
            self.variant_cmap = {
                'Variant': '#92c54c',
                'Sterile': '#475b5f',
                'Negative': '#f7630b',
                'Parent': '#991a4d'}

#             self.default_cmap = [
#                 "#4c78a8",
#                 "#f58518",
#                 "#e45756",
#                 "#72b7b2",
#                 "#54a24b",
#                 "#eeca3b",
#                 "#b279a2",
#                 "#ff9da6",
#                 "#9d755d",
#                 "#bab0ac",
#             ]

            self.default_cmap = ["#475b5f",
                                 "#637478",
                                 "#808e91",
                                 "#9fa9ab",
                                 "#bec5c7",
                                 "#dee2e2",
                                 "#ffffff",
                                 "#f1dadf",
                                 "#e2b5bf",
                                 "#d291a1",
                                 "#c16d84",
                                 "#ae4868",
                                 "#991a4d"]

    def defaults(self):

        if self._user == 'pja':

            hv.opts.defaults(
                hv.opts.Bars(
                    alpha=0.75,
                    cmap=bokeh.palettes.Viridis256,
                    color=hv.Cycle(self.default_cmap),
                    fontsize=dict(title=9),
                    fill_alpha=0.8,
                    line_width=1.5,
                    padding=0.05,
                    toolbar='above',
                ),

                hv.opts.Scatter(
                    size=10,
                    fill_alpha=1,
                    alpha=0.75,
                    cmap=bokeh.palettes.Viridis256,
                    color=hv.Cycle(self.default_cmap),
                    fontsize=dict(title=9),
                    line_width=2,
                    line_color='black',
                    padding=0.05,
                    toolbar='above',
                ),

                hv.opts.HeatMap(
                    color=hv.Cycle(self.default_cmap),
                    fontsize=dict(title=9),
                    line_width=1.5,
                    padding=0.05,
                    cmap='RdBu_r',
                    xaxis='top',
                    toolbar='above',
                    bgcolor='black',
                    frame_height=350,
                    frame_width=350 * 3 // 2,
                ),

                hv.opts.Points(
                    marker='square',
                    size=350 // 8 - (350 // 8)*0.05,
                    line_width=2.5,
                    fill_alpha=0,
                    legend_position='top',
                )
            )

        if self._user == 'ejw':

            hv.opts.defaults(
                hv.opts.Bars(
                    alpha=0.75,
                    cmap=bokeh.palettes.Viridis256,
                    color=hv.Cycle(self.default_cmap),
                    fontsize=dict(title=9),
                    fill_alpha=0.8,
                    line_width=1.5,
                    padding=0.05,
                    toolbar='above',
                ),

                hv.opts.Scatter(
                    line_color=None,
                    width=600,
                    height=450,
                    tools=['hover'],
                    cmap=self.variant_cmap,
                    color=hv.Cycle(self.default_cmap),
                    size=10,
                    fill_alpha=0.6,
                    padding=0.05
                ),


                hv.opts.HeatMap(
                    color=hv.Cycle(self.default_cmap),
                    fontsize=dict(title=9),
                    line_width=1.5,
                    padding=0.05,
                    #                     cmap= bokeh.palettes.diverging_palette(["#475b5f","#ffffff" ], ["#ffffff", "#991a4d"] , 256, 0.5 ),
                    cmap='viridis',
                    xaxis='top',
                    toolbar='above',
                    bgcolor='black',
                    frame_height=350,
                    frame_width=350 * 3 // 2,
                ),


                hv.opts.Points(
                    marker='square',
                    size=350 // 8 - (350 // 8)*0.05,
                    line_width=2.5,
                    fill_alpha=0,
                    legend_position='top',
                )

            )

        return self.variant_cmap, self.default_cmap
