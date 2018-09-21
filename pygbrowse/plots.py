import intervaltree
import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
import scipy
import scipy.signal
from scipy import ndimage
import seaborn

from . import utilities

DEFAULT_ARC_POINTS = 200
DEFAULT_YLABEL_PAD = 50
CHROMOSOME_DIALECT = 'ucsc'


# ToDo: Move this stuff to a separate module so we can delete wholesale when the Ellipse class gets fixed.
def compute_half_arc_points(center, a, b, theta1, theta2, num_points=DEFAULT_ARC_POINTS):
    """
    Computes the coordinates for component points of a polygonal approximation to
    an ellipse for a single quadrant.
    """
    # ToDo: Add input validation to make sure we stay within a single quadrant.
    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)

    for i in range(0, num_points):
        theta = (theta2 - theta1) * (i / max(num_points - 1, 1)) + theta1
        fi = numpy.pi / 2 - numpy.arctan(numpy.tan(theta))
        x = center[0] + a * numpy.cos(fi)
        y = center[1] + b * numpy.sin(fi)
        x_coords[i] = x
        y_coords[i] = y

    return x_coords, y_coords


def draw_arc(ax, center, height, width, theta1=0, theta2=numpy.pi, color='k', direction='down',
             num_points=DEFAULT_ARC_POINTS):
    """
    Since Matplotlib's Arc Patches are broken at the moment, we draw arcs using the ax.plot() method
    instead.
    """
    LEFT_END_THETA = numpy.pi / 2
    RIGHT_END_THETA = numpy.pi * 1.5
    MIDPOINT_THETA = numpy.pi

    vertical_baseline = center[1]

    assert LEFT_END_THETA <= theta1 <= theta2 <= RIGHT_END_THETA

    b = height
    a = width / 2

    # determine how to allocate points
    left_angle_span = min(max(MIDPOINT_THETA - theta1, 0), theta2 - theta1)
    right_angle_span = min(max(theta2 - MIDPOINT_THETA, 0), theta2 - theta1)
    total_angle_span = left_angle_span + right_angle_span
    left_points = int(num_points * left_angle_span / total_angle_span)
    right_points = num_points - left_points

    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)

    if left_points:
        # plot upper left quadrant
        left_theta2 = theta1 + left_angle_span
        x, y = compute_half_arc_points(center=(center[0], 0),
                                       a=a, b=b,
                                       theta1=theta1, theta2=left_theta2,
                                       num_points=left_points)
        x_coords[:left_points] = x[:]
        y_coords[:left_points] = y[:]
    if right_points:
        # plot upper right quadrant
        right_theta1 = theta2 - right_angle_span
        x, y = compute_half_arc_points(center=(center[0], 0),
                                       a=a, b=b,
                                       theta1=right_theta1, theta2=theta2,
                                       num_points=right_points)
        x_coords[left_points:] = x[:]
        y_coords[left_points:] = y[:]

    if direction == 'down':
        y_coords = - y_coords

    y_coords += vertical_baseline

    ax.plot(x_coords, y_coords, color=color)


def draw_visible_arc(ax, center, height, width, ws, we,
                     color='k', direction='down', num_points=DEFAULT_ARC_POINTS):
    """
    Draws a 180 degree elliptical arc truncated by an interval of x coordinates.
    Does not truncate based on y coordinates
    """
    # ToDo: Subtract 1 pi from all coordinates
    LEFT_END_THETA = numpy.pi / 2 + 0.00001
    RIGHT_END_THETA = numpy.pi * 1.5 - 0.00001

    def infer_theta_cutoff(x, arc_center, arc_width):
        a = arc_width / 2
        fi = numpy.arccos((x - arc_center[0]) / a)
        theta = numpy.arctan(1 / numpy.tan(fi)) + numpy.pi
        return theta

    if ws > center[0] - width / 2:
        theta_start = infer_theta_cutoff(x=ws, arc_center=center, arc_width=width)
    else:
        theta_start = LEFT_END_THETA

    if we < center[0] + width / 2:
        theta_end = infer_theta_cutoff(x=we, arc_center=center, arc_width=width)
    else:
        theta_end = RIGHT_END_THETA

    draw_arc(ax=ax, center=center,
             height=height, width=width,
             theta1=theta_start, theta2=theta_end,
             color=color, direction=direction,
             num_points=num_points)


def draw_arc_interaction(ax,
                         xlim, ylim,
                         left_bin_center,
                         right_bin_center,
                         color='k',
                         baseline=0.0,
                         vertical_scaling_factor=1,
                         direction='down',
                         num_points=DEFAULT_ARC_POINTS):
    """
    """
    arc_width = right_bin_center - left_bin_center

    if direction == 'down':
        vertical_span = baseline - ylim[0]
    else:
        vertical_span = ylim[1] - baseline

    arc_height = vertical_span * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor

    draw_visible_arc(ax=ax,
                     center=((left_bin_center + right_bin_center) / 2, baseline),
                     height=arc_height,
                     ws=xlim[0],
                     we=xlim[1],
                     width=arc_width,
                     color=color,
                     direction=direction,
                     num_points=num_points)


class _BrowserSubPlot:
    def __init__(self):
        # self.chrom = None
        # self.ws = None
        # self.we = None
        # self.fig_width = None
        # self.row_height = None
        pass

    # def set_globals(self, chrom, ws, we, fig_width=64, row_height=4):
    #     self.chrom = chrom
    #     self.ws = ws
    #     self.we = we
    #     self.fig_width = fig_width
    #     self.row_height = row_height

    # @property
    # def aspect_ratio(self):
    #     return self.fig_width / self.row_height

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        print('Stub method -- must be overridden by inheritors')


class InteractionPlot(_BrowserSubPlot):
    def __init__(self, interaction_df,
                 bin_size,
                 arc_color=(0.7, 0.3, 0.6),
                 direction='down',
                 baseline=None,
                 vertical_scaling_factor=1,
                 thickness_column=None,
                 show_bin_centers=True,
                 label='Plac-seq'):
        super(InteractionPlot, self).__init__()

        self.interaction_df = interaction_df
        self.bin_size = bin_size
        if baseline is None:
            if direction == 'down':
                baseline = 1
            else:
                baseline = 0
        self.baseline = baseline
        self.vertical_scaling_factor = vertical_scaling_factor
        self.arc_color = arc_color
        self.label = label
        self.direction = direction
        self.show_bin_centers = show_bin_centers
        self.thickness_column = thickness_column


    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        # Filter the interaction DataFrame to interactions with at least one anchor point within the visible window.
        visible_interactions = self.interaction_df.loc[self.interaction_df['chr1'] == chrom]
        left_bin_midpoints = (visible_interactions['end1'] + visible_interactions['start1']) / 2
        right_bin_midpoints = (visible_interactions['end2'] + visible_interactions['start2']) / 2
        left_visible = (left_bin_midpoints >= ws) & (left_bin_midpoints <= we)
        right_visible = (right_bin_midpoints >= ws) & (right_bin_midpoints <= we)
        visible_interactions = visible_interactions.loc[left_visible | right_visible]

        original_ylim = ax.get_ylim()
        ax.set_xlim(ws, we)  # ToDo: Standardize this behavior across all subplot classes

        for interaction_id in visible_interactions.index:
            draw_arc_interaction(ax,
                                 left_bin_center=left_bin_midpoints.loc[interaction_id],
                                 right_bin_center=right_bin_midpoints.loc[interaction_id],
                                 xlim=(ws, we),
                                 ylim=original_ylim,
                                 color=self.arc_color,
                                 baseline=self.baseline,
                                 vertical_scaling_factor=self.vertical_scaling_factor,
                                 direction=self.direction)

        ax.set_xlim(ws, we)
        ax.set_ylim(original_ylim)

        if self.label:
            ax.set_ylabel(self.label)

        if self.show_bin_centers:
            leftmost_tick = numpy.ceil((ws - self.bin_size / 2) / self.bin_size) * self.bin_size + self.bin_size / 2
            rightmost_tick = numpy.floor(
                (we - self.bin_size / 2) / self.bin_size + 1) * self.bin_size + self.bin_size / 2

            ax.set_xticks(numpy.arange(leftmost_tick, rightmost_tick, self.bin_size))
            ax.set_xticklabels([])

            if self.direction == 'down':
                ax.xaxis.set_ticks_position('top')

        ax.set_yticks([])

class BedPlot(_BrowserSubPlot):
    DEFAULT_PATCH_KWARGS = {'linewidth': 1, 'edgecolor': 'k'}
    CHROM_COL_NUM = 0
    START_COL_NUM = 1
    END_COL_NUM = 2

    def __init__(self, interval_data,
                 label='',
                 cmap='RdBu_r',
                 baseline=0.5,
                 color='k',
                 color_by='',
                 display_value='',
                 patch_height=0.5,
                 pad_fraction=0.1,
                 patch_kwargs=None):
        """
        Takes an iterable of tuples in the form:

        (name, DataFrame)

        to be plotted in order using the .plot() method.
        """
        super(BedPlot, self).__init__()

        self.interval_data = interval_data.data
        self.color = color
        self.pad_fraction = pad_fraction

        self.color_by = color_by
        if self.color_by:
            assert self.color_by in self.interval_data.columns, 'Color-by column {} is not in the interval data!'.format(
                self.color_by)
            extent = numpy.abs(self.interval_data[self.color_by]).max()

            self.color_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-extent,
                                                                                              vmax=extent),
                                                             cmap=cmap)
        else:
            self.color_mapper = None

        self.display_value = display_value
        if self.display_value:
            assert self.display_value in self.interval_data.columns, 'Display value column {} is not in the interval data!'.format(
                self.display_value)
        self.label = label
        self.patch_height = patch_height

        self.baseline = baseline
        self.patch_kwargs = self.DEFAULT_PATCH_KWARGS
        if patch_kwargs:
            self.patch_kwargs.update(patch_kwargs)

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        ylim = ax.get_ylim()
        vert_span = ylim[1] - ylim[0]
        
        utilities.add_label(ax=ax, tick=self.baseline, tick_label=self.label, axis='y')
        
        
        visible_intervals = self.interval_data.loc[(self.interval_data.chrom == chrom) & (
                ((ws <= self.interval_data.chromStart) & (self.interval_data.chromStart <= we)) | (
                (ws <= self.interval_data.chromEnd) & (self.interval_data.chromEnd <= we)))]

        for interval_idx in range(visible_intervals.shape[0]):
            this_interval = visible_intervals.iloc[interval_idx]
            start_loc = this_interval['chromStart']
            end_loc = this_interval['chromEnd']
            assert end_loc > start_loc, 'interval end point must be greater than interval start!'
            
            if self.color_by:
                interval_color = self.color_mapper.to_rgba(this_interval[self.color_by])
            else:
                interval_color = self.color

            rec = matplotlib.patches.Rectangle(xy=(start_loc, self.baseline - self.patch_height / 2),
                                               width=end_loc - start_loc,
                                               height=self.patch_height,
                                               facecolor=interval_color,
                                               **self.patch_kwargs)
            ax.add_patch(rec)
            if self.display_value:
                ax.text(x=(start_loc + end_loc) / 2, y=self.baseline,
                        s='{:>0.2}'.format(this_interval[self.display_value]), ha='center')
                        
        utilities.adjust_limits(ax=ax, new_position=self.baseline + self.patch_height / 2, 
                                axis='y', padding_fraction=self.pad_fraction)
        utilities.adjust_limits(ax=ax, new_position=self.baseline - self.patch_height / 2, 
                                axis='y', padding_fraction=self.pad_fraction)
                                
        # print(self.label, ax.get_yticks(), ax.get_yticklabels(), ax.get_ylim())


class WigPlot(_BrowserSubPlot):
    # ToDo: Add support for stranded data
    def __init__(self, genomic_vector_data, label=None, color=None, style='solid', alpha=1.0,  
                center_vector=False, scale_vector_to_plot=False,
                label_rotation=0,
                 # ylim=None,
                 smoothing_bandwidth=0):
        super(WigPlot, self).__init__()  # placeholder since currently the superclass constructor does nothing.
        self.data = genomic_vector_data
        self.color = color
        self.style = style
        self.alpha = alpha
        self.center = center_vector
        self.scale_vector_to_plot = scale_vector_to_plot
        self.label = label
        
        if smoothing_bandwidth:
            self.convolution_kernel = utilities.gaussian_kernel(smoothing_bandwidth)
        else:
            self.convolution_kernel = None
            
        self.label_rotation = label_rotation  
          

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        ylim = ax.get_ylim()

        vert_span = (ylim[1] - ylim[0])
        vert_center = vert_span / 2 + ylim[0]

        this_plot_vector = self.data[chrom].loc[ws:we]

        if self.convolution_kernel is not None:
            this_plot_vector = pandas.Series(
                scipy.signal.convolve(this_plot_vector, self.convolution_kernel, mode='same'),
                index=this_plot_vector.index)

        if self.scale_vector_to_plot:
            this_plot_vector /= (this_plot_vector.max() - this_plot_vector.min())
            this_plot_vector *= vert_span

        if self.center:
            this_plot_vector -= this_plot_vector.mean()
            this_plot_vector += vert_center

        this_plot_vector = this_plot_vector.loc[(this_plot_vector.index >= ws) & (this_plot_vector.index < we)]
        this_plot_vector.name = self.label
               
        if self.style == 'solid':
            ax.fill_between(x=this_plot_vector.index, y1=this_plot_vector, color=self.color, alpha=self.alpha, label=self.label)
        else:
            ax.plot(this_plot_vector.index, this_plot_vector, color=self.color, alpha=self.alpha, label=self.label)
        
        ax.autoscale(enable=True, axis='y')

        # ToDo: Allow labeling either by ylabel or by ax.legend
        if self.label:
            ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)


class GeneModelPlot(_BrowserSubPlot):
    def __init__(self,
                 gene_annotation_data,
                 label='Genes',
                 color='k',
                 feature_height=0.12,
                 chevron_height=0.05,
                 chevron_width=0.04,
                 chevron_spacing=0.10,
                 truncation_size=0.10,
                 utr_endcap_width=0.04,
                 gene_name_fontsize=8):

        super(GeneModelPlot, self).__init__()

        self.gene_annotation_data = gene_annotation_data

        self.color = color
        self.label = label
        self.feature_height = feature_height  # in inches
        self.chevron_height = chevron_height  # in inches
        self.chevron_width = chevron_width  # in inches
        self.chevron_spacing = chevron_spacing  # in inches
        self.truncation_size = truncation_size  # in inches
        self.utr_endcap_width = utr_endcap_width  # in inches
        self.gene_name_fontsize = gene_name_fontsize

    @staticmethod
    def _arrange_genes(gene_data_list):
        """
        Given an iterable of gene data dictionaries,
        returns a list of lists of gene names that
        should be displayed at various levels.
        """
        gene_data_list = sorted(gene_data_list, key=lambda x: x['end'] - x['start'], reverse=True)

        display_levels = [intervaltree.IntervalTree(), ]

        for gene_data in gene_data_list:
            found_home = False
            level_idx = 0
            while not found_home:
                if level_idx >= len(display_levels):
                    display_levels.append(intervaltree.IntervalTree())
                if display_levels[level_idx].overlaps(gene_data['start'], gene_data['end']):
                    level_idx += 1
                else:
                    display_levels[level_idx].addi(gene_data['start'], gene_data['end'], data=gene_data)
                    found_home = True

        return [[gene_interval.data['ID'] for gene_interval in this_level] for this_level in display_levels]

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        # find overlapping genes
        overlapping_genes, overlapping_transcripts, overlapping_components, ids_to_names = self.gene_annotation_data.query(
            chrom, ws, we)
        # overlapping_genes = self.genes.overlapping(chrom, ws, we)
        # overlapping_components = self.components.overlapping(chrom, ws, we)

        gene_display_levels = self._arrange_genes(overlapping_genes.values())
        ax.set_ylim((-0.5, len(gene_display_levels) - 1 + 0.5))

        # convert inches to data coordinates
        chevron_spacing_dt = (we - ws) / (fig_width / self.chevron_spacing)
        chevron_width_dt = (we - ws) / (fig_width / self.chevron_width)
        truncation_width_dt = (we - ws) / (fig_width / self.truncation_size)
        utr_endcap_width_dt = (we - ws) / (fig_width / self.utr_endcap_width)

        feature_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.feature_height)
        chevron_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.chevron_height)
        truncation_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.truncation_size)

        for gene_num, level_genes in enumerate(gene_display_levels):

            # ToDo: make this universal. Divide the gene body into non-overlapping segments, each type of which has a template.

            for gene_id in level_genes:
                gene_data = overlapping_genes[gene_id]
                #                 print(gene_id, gene_data['Name'])

                left_truncated = gene_data['start'] < ws
                right_truncated = gene_data['end'] > we

                visible_gene_start = max(gene_data['start'], ws)
                if left_truncated:
                    visible_gene_start += truncation_width_dt * 2
                visible_gene_end = min(gene_data['end'], we)
                if right_truncated:
                    visible_gene_end -= truncation_width_dt * 2

                ax.plot((visible_gene_start, visible_gene_end), (gene_num, gene_num), color=self.color)
                ax.text(x=(visible_gene_start + visible_gene_end) / 2,
                        y=gene_num + feature_height_dt * 1.5,
                        s=gene_data['Name'],
                        ha='center',
                        fontsize=self.gene_name_fontsize)

                num_chevrons = int(max((visible_gene_end - visible_gene_start) / chevron_spacing_dt, 1))
                chevron_remainder = (visible_gene_end - visible_gene_start) - (num_chevrons - 1) * chevron_spacing_dt

                if gene_data['strand'] == '+':
                    chevron_x_delta = -chevron_width_dt
                else:
                    chevron_x_delta = chevron_width_dt

                for chevron_idx in range(num_chevrons):
                    chevron_x = visible_gene_start + chevron_idx * chevron_spacing_dt + chevron_remainder / 2

                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num + chevron_height_dt),
                            color=self.color)
                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num - chevron_height_dt),
                            color=self.color)

                if left_truncated:
                    y_points = [gene_num, gene_num - truncation_height_dt, gene_num + truncation_height_dt]
                    left_x_point = ws + 1
                    right_x_point = ws + truncation_width_dt + 1

                    x_points = numpy.array([left_x_point, right_x_point, right_x_point])

                    larr1 = matplotlib.patches.Polygon(numpy.vstack([x_points, y_points]).T,
                                                       edgecolor='k',
                                                       facecolor='w',
                                                       fill=True,
                                                       transform=ax.transData,
                                                       zorder=3)
                    larr2 = matplotlib.patches.Polygon(numpy.vstack([x_points + truncation_width_dt, y_points]).T,
                                                       edgecolor='k',
                                                       facecolor='w',
                                                       fill=True,
                                                       transform=ax.transData,
                                                       zorder=3)

                    ax.add_patch(larr1)
                    ax.add_patch(larr2)

                if right_truncated:
                    y_points = [gene_num, gene_num - truncation_height_dt, gene_num + truncation_height_dt]
                    left_x_point = we - truncation_width_dt - 1
                    right_x_point = we - 1

                    x_points = numpy.array([right_x_point, left_x_point, left_x_point])

                    rarr1 = matplotlib.patches.Polygon(xy=numpy.vstack([x_points, y_points]).T,
                                                       edgecolor='k',
                                                       facecolor='w',
                                                       fill=True,
                                                       transform=ax.transData,
                                                       zorder=3)
                    rarr2 = matplotlib.patches.Polygon(numpy.vstack([x_points - truncation_width_dt, y_points]).T,
                                                       edgecolor='k',
                                                       facecolor='w',
                                                       fill=True,
                                                       transform=ax.transData,
                                                       zorder=3)
                    ax.add_patch(rarr1)
                    ax.add_patch(rarr2)

                # Identify components belonging to this gene
                # this_gene_components = set([])
                # for transcript_id in gene_data['transcripts']:
                #     for component_id in overlapping_components:
                #         this_gene_components.add(component_id)

                # plot components
                for component_id in overlapping_components:
                    component_data = overlapping_components[component_id]
                    #                     print('\t', component_id, component_data)
                    if ((component_data['start'] >= visible_gene_start) and (
                            component_data['start'] <= visible_gene_end)) or (
                            (component_data['end'] >= visible_gene_start) and (
                            component_data['end'] <= visible_gene_end)):

                        # ToDo: systematize and condense the following:
                        if component_data['type'] == 'five_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '+':
                                utr_body = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                    width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                    height=feature_height_dt,
                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=self.color)

                            else:
                                utr_body = matplotlib.patches.Rectangle(xy=(
                                    component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                    width=max(
                                        component_data['end'] - component_data[
                                            'start'] - utr_endcap_width_dt, 0),
                                    height=feature_height_dt,
                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=self.color)

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)

                        elif component_data['type'] == 'three_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '-':
                                utr_body = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                    width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                    height=feature_height_dt,
                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=self.color)

                            else:
                                utr_body = matplotlib.patches.Rectangle(xy=(
                                    component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                    width=max(
                                        component_data['end'] - component_data[
                                            'start'] - self.utr_endcap_width, 0),
                                    height=feature_height_dt,
                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=self.color)

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)

                        elif component_data['type'] == 'CDS':
                            cds = matplotlib.patches.Rectangle(
                                xy=(component_data['start'], gene_num - feature_height_dt),
                                width=component_data['end'] - component_data['start'],
                                height=feature_height_dt * 2,
                                facecolor=self.color)
                            ax.add_patch(cds)

        ax.set_yticks([])
        ax.set_ylabel(self.label)


def compute_ax_row_positions(row_heights, ax_spacing=0.1):
    """
    Given a sequence of row heights (in inches), and the size of the space to put between
    axes (in fractions of total row height), returns a list of bottom coordinates
    and heights for each row, suitable for passing to the fig.add_ax() method.
    """
    bottoms = []
    heights = []
    total_canvas_height = numpy.sum(row_heights)
    fig_height = total_canvas_height * (1 + ax_spacing * len(row_heights))
    cur_vertical_pos = 1
    for row_idx in range(len(row_heights)):
        this_row_height = row_heights[row_idx] / fig_height
        cur_vertical_pos -= this_row_height
        heights.append(this_row_height)
        bottoms.append(cur_vertical_pos)
        cur_vertical_pos -= ax_spacing
    return bottoms, heights


class GenomeBrowser:
    VECTOR_LEGEND_LOC = 0

    def __init__(self, subplot_objects):
        """
        Given a 2D nested list of genomic subplot objects it will allow the user to call the .visualize() method
        to generate plots of genomic data.

        :param plot_objects:
        """
        self.subplot_objects = subplot_objects

    def visualize(self, chrom, start, end,
                  fig_width=12,
                  row_heights=1,
                  ax_spacing=0.05,
                  num_xticks=10,
                  seaborn_style=seaborn.axes_style(style='ticks',
                                                   rc={'axes.edgecolor': 'w', 'axes.facecolor': '#EAEAF2'})):
        """
        Generate, display and return a matplotlib.Figure object comprising one or more Axes representing the genomic
        data tracks specified at initialization.

        The region to plot is specified by the parameters chrom, start, and end.

        :param:`fig_width` is specified in inches

        :param:`row_heights` can be specified as a scalar value (in inches), in which case the same row height will be
        used for all subplots, or as an iterable, in which case the row heights will be applied to the subplots
        in order.

        :param chrom:
        :param start:
        :param end:
        :param fig_width:
        :param row_heights:
        :param ax_spacing:
        :param num_xticks:
        :param seaborn_style:
        :return:
        """
        # ToDo: Add gene (or other feature) lookup instead of specifying coordinates.
        start, end = int(start), int(end)

        assert end > start, 'Window end must be greater than window start! Got: {}, {}'.format(start, end)

        # if we receive a scalar here, use it as the height for all rows
        try:
            if len(row_heights) == 1:
                row_heights = row_heights * len(self.subplot_objects)  # treat as a uniform row height
        except TypeError:
            row_heights = [row_heights] * len(self.subplot_objects)  # treat as a uniform row height

        assert len(row_heights) == len(self.subplot_objects)

        span = end - start
        xtick_increment = span / num_xticks
        rounding_increment = 5 * 10 ** numpy.round(numpy.log10(xtick_increment) - 1)
        xtick_increment = utilities.roundto(xtick_increment, rounding_increment)
        num_ticks = int(span / xtick_increment) + 1
        round_start = utilities.roundto(start, rounding_increment)

        seaborn.set_style(seaborn_style)

        fig = plt.figure(len(self.subplot_objects),
                         figsize=(fig_width, numpy.sum(row_heights) * (1 + ax_spacing * len(self.subplot_objects))))
        bottoms, heights = compute_ax_row_positions(row_heights=row_heights, ax_spacing=ax_spacing)

        for ax_idx in range(len(self.subplot_objects)):
            this_ax = fig.add_axes([0, bottoms[ax_idx], 1, heights[ax_idx]])

            if ax_idx == len(self.subplot_objects) - 1:
                this_ax.set_xticks(numpy.arange(num_ticks) * xtick_increment + round_start)
                this_ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
                this_ax.set_xlabel('{} position'.format(chrom))
                
            else:  # clear out xticks but plot objects can override this later
                this_ax.set_xlabel('')
                this_ax.set_xticks([])

            plot_object_subset = self.subplot_objects[ax_idx]

            # Set default plot limits (can be changed by client objects)
            this_ax.set_ylim(0, 1)
            this_ax.set_xlim(start, end)

            for plot_object in plot_object_subset:
                plot_object.plot(this_ax, chrom=chrom, ws=start, we=end, fig_width=fig_width,
                                 row_height=row_heights[ax_idx])

            # ToDo: Refactor legend code to get colors and names from objects not from axes handles.
            # if len(this_ax.get_legend_handles_labels()[1]):
            #     this_ax.legend(loc=self.VECTOR_LEGEND_LOC)

        return fig

class HicPlot:
    def __init__(self, data, label='', vertical_scale=1, cmap='YlOrRd', label_rotation=0, transform=lambda x: x**2, max_masked_diag=2):
        self.data = data
        self.cmap = cmap
        self.transform = transform # ToDo: Move the transform to the Data provider
        self.max_masked_diag = max_masked_diag
        self.label = label
        self.label_rotation = label_rotation
        self.vertical_scale = vertical_scale
        
    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        visible_start_bin = utilities.roundto(ws, self.data.bin_size)
        visible_end_bin = utilities.roundto(we, self.data.bin_size)
        visible_span = visible_end_bin - visible_start_bin

        data_start_bin = visible_start_bin - visible_span // 2
        data_end_bin = visible_end_bin + visible_span // 2

        plot_data = self.data.query(chrom, data_start_bin, data_end_bin)
        
        for diag in range(self.max_masked_diag):
            plot_data.values[utilities.diag_indices(plot_data.shape[0], diag)] = 0

        plot_data = pandas.DataFrame(ndimage.rotate(plot_data, 45, reshape=False),
                                     index=plot_data.index, columns=plot_data.columns)
#         print(plot_data.shape)
        
        # Trim back to visible area
        plot_data = plot_data.loc[visible_start_bin:visible_end_bin,visible_start_bin:visible_end_bin]
        
        # Only show upper diagonal
        plot_data = plot_data.iloc[:-plot_data.shape[0] // 2,:]
        
    
        plot_data = self.transform(plot_data)
        
        # Re-index plot_data to allow it to play nicely with ax limits

        ax.set_ylim(0, plot_data.shape[1])
        ax.imshow(plot_data, cmap=self.cmap, aspect='auto', extent=(ws, we, 0, plot_data.shape[1]))
        ax.set_xticks([])
        ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)
#         print('Done on {}'.format(ax))
        return plot_data
        
        
def match_ylims(fig, ax_nums):
    """
    Will make the upper ylim of each of the numbered axes of :param fig: listed in
    :param ax_nums: equal to the maximum found in any of the numbered axes.
    """
    max_extent = max([fig.get_axes()[ax_num].get_ylim()[1] for ax_num in ax_nums])
    for ax_num in ax_nums:
        fig.get_axes()[ax_num].set_ylim((0, max_extent))        