import intervaltree
import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
import scipy
import scipy.signal
import seaborn
from pgtools import toolbox

DEFAULT_ARC_POINTS = 200


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
        self.chrom = None
        self.ws = None
        self.we = None
        self.fig_width = None
        self.row_height = None

    def set_globals(self, chrom, ws, we, fig_width=64, row_height=4):
        self.chrom = chrom
        self.ws = ws
        self.we = we
        self.fig_width = fig_width
        self.row_height = row_height

    @property
    def aspect_ratio(self):
        return self.fig_width / self.row_height


class InteractionPlot(_BrowserSubPlot):
    def __init__(self, interaction_df,
                 bin_size,
                 arc_color='k',
                 direction='up',
                 baseline=0,
                 vertical_scaling_factor=1,
                 thickness_column=None,
                 show_bin_centers=True,
                 label='Plac-seq'):
        super(InteractionPlot, self).__init__()

        self.interaction_df = interaction_df
        self.bin_size = bin_size
        self.baseline = baseline
        self.vertical_scaling_factor = vertical_scaling_factor
        self.arc_color = arc_color
        self.label = label
        self.direction = direction
        self.show_bin_centers = show_bin_centers
        self.thickness_column = thickness_column

    def plot(self, ax):
        # Filter the interaction DataFrame to interactions with at least one anchor point within the visible window.
        visible_interactions = self.interaction_df.loc[self.interaction_df['chr1'] == self.chrom]
        left_bin_midpoints = (visible_interactions['end1'] + visible_interactions['start1']) / 2
        right_bin_midpoints = (visible_interactions['end2'] + visible_interactions['start2']) / 2
        left_visible = (left_bin_midpoints >= self.ws) & (left_bin_midpoints <= self.we)
        right_visible = (right_bin_midpoints >= self.ws) & (right_bin_midpoints <= self.we)
        visible_interactions = visible_interactions.loc[left_visible | right_visible]

        original_ylim = ax.get_ylim()
        ax.set_xlim(self.ws, self.we)

        for interaction_id in visible_interactions.index:
            draw_arc_interaction(ax,
                                 left_bin_center=left_bin_midpoints.loc[interaction_id],
                                 right_bin_center=right_bin_midpoints.loc[interaction_id],
                                 xlim=(self.ws, self.we),
                                 ylim=original_ylim,
                                 color=self.arc_color,
                                 baseline=self.baseline,
                                 vertical_scaling_factor=self.vertical_scaling_factor,
                                 direction=self.direction)

        ax.set_xlim(self.ws, self.we)
        ax.set_ylim(original_ylim)
        if self.label:
            ax.set_ylabel(self.label)

        if self.show_bin_centers:
            leftmost_tick = numpy.ceil((
                                               self.ws - self.bin_size / 2) / self.bin_size) * self.bin_size + self.bin_size / 2
            rightmost_tick = numpy.floor((
                                                 self.we - self.bin_size / 2) / self.bin_size + 1) * self.bin_size + self.bin_size / 2

            ax.set_xticks(numpy.arange(leftmost_tick, rightmost_tick, self.bin_size))
            ax.set_xticklabels([])

            if self.direction == 'down':
                ax.xaxis.set_ticks_position('top')


class IntervalLines(_BrowserSubPlot):
    # ToDo: Combine lines and box plot classes or superclass them
    def __init__(self, named_crds, crd_cmap='Purples', crd_outline_color=None, linewidth=5):
        """
        Takes an iterable of tuples in the form:

        (name, DataFrame)

        to be plotted in order using the .plot() method.
        """
        super(IntervalLines, self).__init__()

        self.named_crds = named_crds
        min_val = min([crd_df['score'].min() for crdset_name, crd_df in named_crds])
        max_val = max([crd_df['score'].max() for crdset_name, crd_df in named_crds])

        self.crd_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min_val,
                                                                                        vmax=max_val),
                                                       cmap=crd_cmap)
        self.crd_outline_color = crd_outline_color
        self.linewidth = linewidth

    def plot(self, ax):
        for vertical_offset, (crdset_name, crd_df) in enumerate(self.named_crds):
            for crd_name in crd_df.index:
                start_loc = crd_df.loc[crd_name, 'chromStart']
                end_loc = crd_df.loc[crd_name, 'chromEnd']
                if start_loc <= self.we or end_loc >= self.ws:
                    interval_color = self.crd_mapper.to_rgba(crd_df.loc[crd_name]['score'])
                    ax.plot((start_loc,
                             end_loc),
                            (vertical_offset,
                             vertical_offset),
                            color=interval_color,
                            linewidth=5)
        ax.set_yticks(numpy.arange(len(self.named_crds)))
        ax.set_yticklabels([crdset_name for crdset_name, crd_df in self.named_crds])
        ax.set_ylim(-1, len(self.named_crds))


class IntervalBoxes(_BrowserSubPlot):
    DEFAULT_PATCH_KWARGS = {'linewidth': 1, 'edgecolor': 'k'}

    def __init__(self, named_crds,
                 crd_cmap='Purples',
                 height=None,
                 annotation_col='',
                 patch_kwargs={}):
        """
        Takes an iterable of tuples in the form:

        (name, DataFrame)

        to be plotted in order using the .plot() method.
        """
        super(IntervalBoxes, self).__init__()

        self.named_crds = named_crds[::-1]
        extent = max([numpy.abs(crd_df['score']).max() for crdset_name, crd_df in named_crds])

        self.color_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-extent,
                                                                                          vmax=extent),
                                                         cmap=crd_cmap)

        self.height = height
        self.annotation_column = annotation_col
        self.patch_kwargs = self.DEFAULT_PATCH_KWARGS
        self.patch_kwargs.update(patch_kwargs)

    def plot(self, ax):
        ylim = ax.get_ylim()
        vert_span = ylim[1] - ylim[0]
        num_tracks = len(self.named_crds)
        if self.height is None:
            self.height = vert_span / (len(self.named_crds) * 2 - 1)

        vert_pad = (vert_span - num_tracks * self.height) / (num_tracks + 1)
        canvas_span = vert_span - (2 * vert_pad) - self.height

        yticks = []
        yticklabels = []
        for crdset_idx, (crdset_name, crd_df) in enumerate(self.named_crds):
            if num_tracks == 1:
                bottom = vert_pad
            else:
                bottom = vert_pad + (crdset_idx / (num_tracks - 1) * canvas_span)

            vert_midpoint = bottom + self.height / 2

            yticks.append(vert_midpoint)
            yticklabels.append(crdset_name)

            this_crds = crd_df.loc[(crd_df.chrom == self.chrom) & (
                    ((self.ws <= crd_df.chromStart) & (crd_df.chromStart <= self.we)) | (
                    (self.ws <= crd_df.chromEnd) & (crd_df.chromEnd <= self.we)))]

            for crd_name in this_crds.index:
                start_loc = this_crds.loc[crd_name, 'chromStart']
                end_loc = this_crds.loc[crd_name, 'chromEnd']
                interval_color = self.color_mapper.to_rgba(crd_df.loc[crd_name]['score'])
                #                     print('plotting {} in {}'.format(crd_name, interval_color))

                rec = matplotlib.patches.Rectangle(xy=(start_loc, bottom),
                                                   width=end_loc - start_loc,
                                                   height=self.height,
                                                   facecolor=interval_color,
                                                   **self.patch_kwargs)
                ax.add_patch(rec)
                if self.annotation_column and self.annotation_column in crd_df.columns:
                    #                         print((start_loc+end_loc)/2, midpoint)
                    ax.text(x=(start_loc + end_loc) / 2, y=vert_midpoint,
                            s='{:>0.2}'.format(crd_df.loc[crd_name, self.annotation_column]), ha='center')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)


#         ax.set_ylim(-1, len(self.named_crds))


class VectorPlot(_BrowserSubPlot):
    def __init__(self, vector_series, label=None, color=None, center=False, scale=False, ylim=None,
                 convolution_kernel=None):
        super(VectorPlot, self).__init__()
        self.color = color
        self.vector = vector_series
        self.center = center
        self.scale = scale
        self.label = label
        self.convolution_kernel = convolution_kernel
        self.ylim = ylim

    def plot(self, ax):
        if self.ylim:
            ylim = self.ylim
        else:
            ylim = ax.get_ylim()

        vert_span = (ylim[1] - ylim[0])
        vert_center = vert_span / 2 + ylim[0]

        this_plot_vector = self.vector.copy()

        if self.convolution_kernel is not None:
            this_plot_vector = pandas.Series(
                scipy.signal.convolve(this_plot_vector, self.convolution_kernel, mode='same'),
                index=this_plot_vector.index)

        if self.scale:
            this_plot_vector /= (this_plot_vector.max() - this_plot_vector.min())
            this_plot_vector *= vert_span

        if self.center:
            this_plot_vector -= this_plot_vector.mean()
            this_plot_vector += vert_center
        #         print(self.center, self.scale, vert_center, vert_span, this_plot_vector.min(), this_plot_vector.max())
        this_plot_vector = this_plot_vector.loc[
            (this_plot_vector.index >= self.ws) & (this_plot_vector.index < self.we)]
        ax.plot(this_plot_vector.index, this_plot_vector, color=self.color, label=self.label)
        ax.set_ylim(ylim)


# class FeatureStats(_BrowserSubPlot):
#     def __init__(self, features_df, annotated_regions_df,
#                  plot_mean=True, plot_var=False, plot_sd=False,
#                  plot_relvar=False, plot_cv=True, znorm_tracks=False):
#         super(FeatureStats, self).__init__()
#
#         self.displayed_features = features_df.T.copy()
#         self.displayed_features.columns = (annotated_regions_df.loc[self.displayed_features.columns].Start + \
#                                            annotated_regions_df.loc[self.displayed_features.columns].End) / 2
#
#         self.znorm_tracks = znorm_tracks
#         self.plot_mean = plot_mean
#         self.plot_var = plot_var
#         self.plot_sd = plot_sd
#         self.plot_relvar = plot_relvar
#         self.plot_cv = plot_cv
#
#     def plot(self, ax):
#         def plot_possibly_znormed(track_data):
#             if self.znorm_tracks: track_data = toolbox.znorm(track_data)
#             ax.plot(track_data)
#
#         if self.plot_mean or self.plot_cv:
#             feature_mean = self.displayed_features.mean(axis=0)
#             if self.plot_mean:
#                 plot_possibly_znormed(feature_mean)
#         if self.plot_var or self.plot_relvar:
#             feature_var = self.displayed_features.var(axis=0)
#             if self.plot_var:
#                 plot_possibly_znormed(feature_var)
#         if self.plot_sd or self.plot_cv:
#             feature_sd = self.displayed_features.std(axis=0)
#             if self.plot_sd:
#                 plot_possibly_znormed(feature_sd)
#         if self.plot_relvar:
#             feature_relvar = feature_var / feature_mean
#             plot_possibly_znormed(feature_relvar)
#         if self.plot_cv:
#             feature_cv = feature_sd / feature_mean
#             plot_possibly_znormed(feature_cv)

class GeneModels(_BrowserSubPlot):
    def __init__(self, genome, label=None, color='k',
                 feature_height=0.12,
                 chevron_height=0.05,
                 chevron_width=0.04,
                 chevron_spacing=0.10,
                 truncation_size=0.10,
                 utr_endcap_width=0.04,
                 gene_name_fontsize=8):
        super(GeneModels, self).__init__()

        self.color = color
        self.genome = genome
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

    def plot(self, ax):
        # find overlapping genes
        overlapping_genes = self.genome.genes.overlapping(self.chrom, self.ws, self.we)
        overlapping_components = self.genome.components.overlapping(self.chrom, self.ws, self.we)

        gene_display_levels = self._arrange_genes(overlapping_genes.values())
        ax.set_ylim((-0.5, len(gene_display_levels) - 1 + 0.5))

        # convert inches to data coordinates
        chevron_spacing_dt = (self.we - self.ws) / (self.fig_width / self.chevron_spacing)
        chevron_width_dt = (self.we - self.ws) / (self.fig_width / self.chevron_width)
        truncation_width_dt = (self.we - self.ws) / (self.fig_width / self.truncation_size)
        utr_endcap_width_dt = (self.we - self.ws) / (self.fig_width / self.utr_endcap_width)

        feature_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.feature_height)
        chevron_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.chevron_height)
        truncation_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.truncation_size)

        for gene_num, level_genes in enumerate(gene_display_levels):

            # ToDo: make this universal. Divide the gene body into non-overlapping segments, each type of which has a template.

            for gene_id in level_genes:
                gene_data = overlapping_genes[gene_id]
                #                 print(gene_id, gene_data['Name'])

                left_truncated = gene_data['start'] < self.ws
                right_truncated = gene_data['end'] > self.we

                visible_gene_start = max(gene_data['start'], self.ws)
                if left_truncated:
                    visible_gene_start += truncation_width_dt * 2
                visible_gene_end = min(gene_data['end'], self.we)
                if right_truncated:
                    visible_gene_end -= truncation_width_dt * 2

                ax.plot((visible_gene_start, visible_gene_end), (gene_num, gene_num), color=self.color)
                ax.text(x=(visible_gene_start + visible_gene_end) / 2,
                        y=gene_num + feature_height_dt * 1.5,
                        s=gene_data['Name'],
                        ha='center',
                        #                         rotation=10,
                        fontsize=self.gene_name_fontsize)

                num_chevrons = int(max((visible_gene_end - visible_gene_start) / chevron_spacing_dt, 1))
                chevron_remainder = (visible_gene_end - visible_gene_start) - (num_chevrons - 1) * chevron_spacing_dt

                if gene_data['strand'] == '+':
                    chevron_x_delta = -chevron_width_dt
                else:
                    chevron_x_delta = chevron_width_dt

                for chevron_idx in range(num_chevrons):
                    #                     if num_chevrons =
                    chevron_x = visible_gene_start + chevron_idx * chevron_spacing_dt + chevron_remainder / 2

                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num + chevron_height_dt),
                            color=self.color)
                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num - chevron_height_dt),
                            color=self.color)

                if left_truncated:
                    y_points = [gene_num, gene_num - truncation_height_dt, gene_num + truncation_height_dt]
                    left_x_point = self.ws + 1
                    right_x_point = self.ws + truncation_width_dt + 1

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
                    left_x_point = self.we - truncation_width_dt - 1
                    right_x_point = self.we - 1

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
                this_gene_components = set([])
                for transcript_id in gene_data['transcripts']:
                    for component_id in self.genome.transcripts[transcript_id]['components']:
                        if component_id in overlapping_components:
                            #                         print('\t' + component_id)
                            this_gene_components.add(component_id)

                # plot components
                for component_id in this_gene_components:
                    component_data = self.genome.components[component_id]
                    #                     print('\t', component_id, component_data)
                    if ((component_data['start'] >= visible_gene_start) and (
                            component_data['start'] <= visible_gene_end)) or (
                            (component_data['end'] >= visible_gene_start) and (
                            component_data['end'] <= visible_gene_end)):

                        # ToDo: systematize and condense the following
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

            ax.set_yticks(range(len(gene_display_levels)))


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


def visualize(plot_objects,
              chrom='',
              ws=0, we=0,
              genome=None,
              fig_width=12,
              row_heights=1,
              ax_spacing=0.1,
              show_vector_legend=False,
              vector_legend_loc=0,
              #                           show_x_span=False,
              #                           x_span_kwargs={},
              num_ticks=10,
              seaborn_style=seaborn.axes_style(style='ticks', rc={
                  'axes.edgecolor': 'w',
                  'axes.facecolor': '#EAEAF2',
              })):
    if we == 0:
        we = genome.contig_lengths[chrom]

    ws, we = int(ws), int(we)

    assert we > ws, 'Window end must be greater than window start! Got: {}, {}'.format(ws, we)
    try:
        if len(row_heights) == 1:
            row_heights = row_heights * len(plot_objects)  # treat as a uniform row height
    except TypeError:
        row_heights = [row_heights] * len(plot_objects)  # treat as a uniform row height

    assert len(row_heights) == len(plot_objects)

    span = we - ws
    xtick_increment = span / num_ticks
    rounding_increment = 5 * 10 ** numpy.round(numpy.log10(xtick_increment) - 1)
    xtick_increment = toolbox.roundto(xtick_increment, rounding_increment)
    num_ticks = int(span / xtick_increment) + 1
    round_start = toolbox.roundto(ws, rounding_increment)

    seaborn.set_style(seaborn_style)

    fig = plt.figure(len(plot_objects),
                     figsize=(fig_width, numpy.sum(row_heights) * (1 + ax_spacing * len(plot_objects))))
    bottoms, heights = compute_ax_row_positions(row_heights=row_heights, ax_spacing=ax_spacing)

    for ax_idx in range(len(plot_objects)):
        this_ax = fig.add_axes([0, bottoms[ax_idx], 1, heights[ax_idx]])

        if ax_idx == len(plot_objects) - 1:
            this_ax.set_xticks(numpy.arange(num_ticks) * xtick_increment + round_start)
            this_ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
            this_ax.set_xlabel('{} position'.format(chrom))
        else:  # clear out xticks but plot objects can override this later
            this_ax.set_xlabel('')
            this_ax.set_xticks([])

        plot_object_subset = plot_objects[ax_idx]

        # Set default plot limits (can be changed by client objects)
        this_ax.set_ylim(0, 1)
        this_ax.set_xlim(ws, we)

        for plot_object in plot_object_subset:
            plot_object.set_globals(chrom=chrom, ws=ws, we=we, fig_width=fig_width, row_height=row_heights[ax_idx])
            plot_object.plot(this_ax)

        # ToDo: Refactor legend code to get colors and names from objects not from axes handles.
        if show_vector_legend and len(this_ax.get_legend_handles_labels()[1]):
            this_ax.legend(loc=vector_legend_loc)

    return fig
