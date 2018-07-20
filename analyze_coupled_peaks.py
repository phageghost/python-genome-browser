import collections
import itertools

import intervaltree
import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
import scipy
import scipy.signal
import seaborn
from pgtools import myplots
from pgtools import toolbox
from pgtools.genomicwrappers import GENE_TYPES, COMPONENT_TYPES
from pgtools.toolbox import log_print

DEFAULT_GENOME = object()
NUM_ARC_POINTS = 200
TRAIN_CHROM = 'chr18'


def clean_col_names(columns, name_mapper):
    return [name_mapper[col.split(' ')[0].split('/')[-1]] for col in columns]


def load_annotations(data_type, condition, chrom, remove_clonal=False, rip_norm=True, znorm=False, log_transform=True):
    annotated_peaks_fname = call_coupled_peaks.generate_matching_peak_annotation_fname(data_type, condition,
                                                                                           annotation_basepath=ANNOTATED_PEAK_FOLDER)
    log_print('Loading annotated peak file from {}'.format(annotated_peaks_fname))
    anno_peaks = pandas.read_csv(annotated_peaks_fname, sep='\t', index_col=0).sort_values(['Chr', 'Start'])
    anno_peaks.index = [anno_peaks.loc[this_peak, 'Chr'] + '-' +  str(anno_peaks.loc[this_peak, 'Start']) for this_peak in anno_peaks.index]
    anno_peaks = call_coupled_peaks.normalize_peak_data(anno_peaks, 
                                                        remove_clonal=remove_clonal, 
                                                        rip_norm=rip_norm, 
                                                        znorm=znorm,
                                                        log_transform=log_transform)
    renamer = toolbox.invert_dict(construct_data_renamer(data_type, condition))
#     print(renamer, anno_peaks.columns[6:])
    anno_peaks.columns = list(anno_peaks.columns[:6]) + clean_col_names(anno_peaks.columns[6:], renamer)
    anno_peaks = anno_peaks.loc[anno_peaks.Chr == chrom]
    return anno_peaks


def extract_features(annotated_peaks):
    chrom_features = annotated_peaks.iloc[:,6:]
    chrom_features.index.name = 'peak'
    chrom_features.columns.name = 'strain'
    return chrom_features


def construct_data_renamer(data_type, condition, input_dna=False, strainlist=('C57', 'BALB', 'NOD', 'PWK', 'SPRET')):
    if data_type.startswith('H3'):
        if condition == 'KLA':
            condition = 'KLA_1h'
        data_type = 'ChIP_' + data_type
        
        if input_dna :
            data_type = 'Input_FA'
        
        return {strain:'{}_BMDM_{}_{}_BC_merged'.format(strain, data_type, condition) for strain in strainlist}
    else:
        if input_dna :
            template = 'pooled_Input_{}_{}_{}'
        else:
            template = 'pooled_IDR_{}_{}_{}'
        return {strain:template.format(strain, data_type, condition) for strain in strainlist}


def dd():
    return collections.defaultdict(lambda:{})


def compute_feature_deltas(feature_df, reference_column='C57'):
    query_cols = [col for col in feature_df.columns if col != reference_column]
    return (feature_df.loc[:, reference_column].T - feature_df.loc[:,query_cols].T).T


def compute_feature_deltas_from_mean(feature_df,):
    return(feature_df.T - feature_df.T.mean()).T


def znorm_rows(feature_df):
    return ((feature_df.T - feature_df.T.mean()) / feature_df.T.std()).T


def znorm_cols(feature_df):
    return ((feature_df - feature_df.mean()) / feature_df.std())


def extract_crd_features(crd_df, feature_df):
    crd_features = {}
    for crd_name in crd_df.index:
        component_peaks = crd_name.split('_')
        crd_features[crd_name] = feature_df.loc[component_peaks]
    return crd_features


def extract_crd_components(crd_df):
    """
    Returns a set of all component peak identifiers comprising the CRDs in :param:`crd_df`
    """
    all_crd_components = set([])
    for crd_name in crd_df.index:
        all_crd_components.update(crd_name.split('_'))
    return all_crd_components


def melt_meta_profiles(feature_df, crd_df, value_name):
    """
    Given a DataFrame of features and CRDs,
    returns a long-form pandas Dataframe of the 
    strain profiles of the component peaks
    suitable for plotting with seaborn barplot
    or violinplot
    """
    crd_components = extract_crd_components(crd_df)
    non_crd_components = set(feature_df.index).difference(crd_components)  
    
    crd_pos_features = feature_df.loc[crd_components]
    crd_neg_features = feature_df.loc[non_crd_components]

    crd_pos_features_long = pandas.melt(crd_pos_features, var_name='strain', value_name=value_name)
    crd_pos_features_long['CRD status'] = 'inside CRD'
    crd_neg_features_long = pandas.melt(crd_neg_features, var_name='strain', value_name=value_name)
    crd_neg_features_long['CRD status'] = 'outside CRDs'

    return pandas.concat((crd_pos_features_long, crd_neg_features_long), axis=0)


def analyze_meta_profiles(feature_df, crd_df, figsize=(8,3)):
    
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    
    standard_feature_profiles = melt_meta_profiles(feature_df, crd_df, value_name='log2 tag count')
    seaborn.violinplot(data=standard_feature_profiles, x='strain', y='log2 tag count', hue='CRD status', ax=ax[0])    
    
    rowznormed_features = znorm_rows(feature_df)
    rowznormed_feature_profiles = melt_meta_profiles(rowznormed_features, crd_df, value_name='log2 tag z-score')
    seaborn.violinplot(data=rowznormed_feature_profiles, x='strain', y='log2 tag z-score', hue='CRD status', ax=ax[1])     
    
    ax[0].get_legend().set_visible(False)
    ax[1].legend(bbox_to_anchor=(1.0,0.5), )
    
    fig.tight_layout()
    return fig

# def analyze_meta_profiles(feature_df, crd_df, figsize=(8,3)):
    
#     fig, ax = plt.subplots(1, 1, figsize=figsize)
    
#     standard_feature_profiles = melt_meta_profiles(feature_df, crd_df, value_name='log2 tag count')
#     seaborn.violinplot(data=standard_feature_profiles, x='strain', y='log2 tag count', hue='CRD status', ax=ax)    
    
# #     rowznormed_features = znorm_rows(feature_df)
# #     rowznormed_feature_profiles = melt_meta_profiles(rowznormed_features, crd_df, value_name='tag Z')
# #     seaborn.violinplot(data=rowznormed_feature_profiles, x='strain', y='tag Z', hue='CRD status', ax=ax[1])     
    
#     ax.get_legend().set_visible(False)
#     ax.legend(bbox_to_anchor=(1.0,0.5), )
    
#     fig.tight_layout()
#     return fig


def process_features(annotated_peak_df, rip_norm, znorm, log_transform, chrom='chr18'):
    this_features = call_coupled_peaks.normalize_peak_data(annotated_peak_df.loc[annotated_peak_df.Chr == chrom], 
                                                           rip_norm=rip_norm,
                                                           znorm=znorm,
                                                           log_transform=log_transform).iloc[:,6:]
    this_features.columns = clean_col_names(this_features.columns, name_mapper=toolbox.invert_dict(atac_strain_renamer))
    return this_features

def count_tags_in_regions(bed_fname, output_fname, data_tag_folders, genome='mm10'):
    cmd_line = ['annotatePeaks.pl', bed_fname, genome, '-size given', '-d', ' '.join(data_tag_folders), '>', output_fname]
    cmd_line = ' '.join(cmd_line)
    print(cmd_line)

    
def convert_chrom_coords_to_peak_coords(annotated_peaks_df, query_chrom, query_start, query_end):
    this_chrom_peaks = annotated_peaks_df.loc[annotated_peaks_df.Chr == 'chr18']
    this_chrom_peaks = this_chrom_peaks.sort_values('Start')

    peak_start = numpy.max(numpy.nonzero(numpy.less_equal(this_chrom_peaks.Start, query_start)))
    peak_end = numpy.min(numpy.nonzero(numpy.greater_equal(this_chrom_peaks.End, query_end)))
    
    return peak_start, peak_end


def determine_size_distribution(sizes):
    counts, values = numpy.histogram(sizes, bins=numpy.arange(sizes.min(), sizes.max()+2) - 0.5)

    values = numpy.arange(sizes.min(), sizes.max()+1)
    frequencies = counts / counts.sum()
    return values, frequencies

def find_uncoupled_blocks(annotated_peaks, coupled_peak_set):
    """
    Find contiguous blocks of uncoupled peaks
    """
    
    uncoupled_blocks = []

    for chrom, chrom_data in annotated_peaks.groupby('Chr'):
        chrom_data = chrom_data.sort_values('Start')
        n = chrom_data.shape[0]
        i = 0
        while i < n:
            while i < n and chrom_data.index[i] in coupled_peak_set:
                i += 1
#                 print('\tin', i,)
                
            if i < n:
                start = i

                while i < n and chrom_data.index[i] not in coupled_peak_set:
#                     print('\tout', i,)
                    i += 1
                end = i
    #             print(start, end, chrom_data.index[start:end])
#                 print(start,end, n)
        
                assert sum([peak_name in coupled_peak_set for peak_name in chrom_data.index[start:end]]) == 0
                size = end - start - 1
                uncoupled_blocks.append((chrom, start, end, size))
    return uncoupled_blocks


def sample_uncoupled_peaks(annotated_peaks, coupled_peak_df, num_samples=1000):
    print('Finding {} size-matched uncoupled regions'.format(num_samples))
    coupled_peak_set = set([])
    
    coupled_size_values = numpy.zeros(coupled_peak_df.shape[0])
    for i, name in enumerate(coupled_peak_df.index):
        split_names = name.split('_')
        coupled_size_values[i] = len(split_names)
        coupled_peak_set.update(set(split_names))
    
    coupled_sizes, coupled_size_frequencies = determine_size_distribution(coupled_size_values)
#     print(coupled_sizes, coupled_size_frequencies)

    uncoupled_blocks = find_uncoupled_blocks(annotated_peaks, coupled_peak_set)
#     print(uncoupled_blocks)
    uncoupled_sizes = numpy.array([region[3] for region in uncoupled_blocks], dtype=int)
    sampled_sizes = numpy.random.choice(a=coupled_sizes, p=coupled_size_frequencies, size=num_samples).astype(int)
    
#     print(uncoupled_sizes, sampled_sizes)
    
    annotations_grouped = {chrom:chrom_data.sort_values('Start') for chrom, chrom_data in annotated_peaks.groupby('Chr')}
    
    sampled_peak_data = {}

    for size in sampled_sizes:
        # sample the uncoupled blocks proportionally
        ways = numpy.maximum((uncoupled_sizes - size) + 1,0).astype(int) # number of ways this size can fit within the sampled block
        weights = ways / ways.sum()
#         print(ways, weights)
        sampled_block_idx = numpy.random.choice(a=numpy.arange(len(uncoupled_sizes)), p=weights, size=1)[0]
#         print(sampled_block_idx)
        sampled_block_chrom, sampled_block_start, sampled_block_end, sampled_block_size  = uncoupled_blocks[sampled_block_idx]
        
        starting_peak_idx = numpy.random.randint(0, ways[sampled_block_idx]) + sampled_block_start
        ending_peak_idx = starting_peak_idx + size 
        
#         print(sampled_block_idx, uncoupled_blocks[sampled_block_idx], starting_peak_idx, ending_peak_idx)

        this_chrom_annotated_peaks = annotations_grouped[sampled_block_chrom]
#         print(this_chrom_annotated_peaks.shape)
        
        compound_name = '_'.join(this_chrom_annotated_peaks.index[starting_peak_idx:ending_peak_idx])
        
        assert sum([peak_name in coupled_peak_set for peak_name in this_chrom_annotated_peaks.index[starting_peak_idx:ending_peak_idx]]) == 0, list(zip(this_chrom_annotated_peaks.index[starting_peak_idx:ending_peak_idx], 
                                                            [peak_name in coupled_peak_set for peak_name in this_chrom_annotated_peaks.index[starting_peak_idx:ending_peak_idx]]) )
        starting_chrom_coord = this_chrom_annotated_peaks.iloc[starting_peak_idx].loc['Start']
        ending_chrom_coord = this_chrom_annotated_peaks.iloc[ending_peak_idx].loc['End']
        
#         print(sampled_block_chrom, starting_chrom_coord, ending_chrom_coord, compound_name)
        
        sampled_peak_data[compound_name] = {'chrom':sampled_block_chrom,
                                            'chromStart':starting_chrom_coord,
                                            'chromEnd':ending_chrom_coord,
                                            'score:':0,
                                            'strand':'+'}
        
    sampled_peak_data = pandas.DataFrame(sampled_peak_data).T
    sampled_peak_data = sampled_peak_data.loc[:, ['chrom', 'chromStart', 'chromEnd', 'score', 'strand']]
    sampled_peak_data = sampled_peak_data.sort_values(['chrom', 'chromStart'])
    return sampled_peak_data


def tuples_from_df(df):
    results = {}
    for chrom, chrom_data in df.groupby('chrom'):
        results[chrom] = list(zip(chrom_data.index, chrom_data.chromStart, chrom_data.chromEnd))
    return results


def find_overlapping_genes(crd_df, genome):
    gene_tuples = genome.features.as_tupledict()
    crd_tuples = tuples_from_df(crd_df)
    overlapping_genes = set([])
    for chrom in crd_tuples:
        if chrom in gene_tuples:
            this_overlaps = compute_interval_overlaps(crd_tuples[chrom], gene_tuples[chrom])
            overlapping_genes.update([tup[1] for tup in this_overlaps])
    return overlapping_genes


def count_tags_in_regions(bed_fname, output_fname, data_tag_folders, genome='mm10', cmdlineonly=True):
    cmd_line = ['annotatePeaks.pl', bed_fname, genome, '-size given', '-d', ' '.join(data_tag_folders), '>', output_fname]
    cmd_line = ' '.join(cmd_line)
    print(cmd_line)

    if cmdlineonly:
        print(cmd_line)
    else:
        print(subprocess.check_output(cmd_line, shell=True))

        
def analyze_peak_sizes(crd_df, stats, genome, fig_saver):
    crd_sizes = {}
    for chrom in genome.contig_names:
        crd_sizes.update({peak_name:len(peak_name.split('_')) for peak_name in crd_df.loc[crd_df['chrom'] == chrom].index})
    crd_sizes = pandas.Series(crd_sizes).sort_index()
    crd_sizes.name = 'CRD sizes'
    
    fig, ax = plt.subplots(1, figsize=(4,3))
    seaborn.kdeplot(crd_sizes, ax=ax)
    ax.set_xlabel('# Member peaks')
    ax.set_ylabel('relative frequency')
    ax.set_title('CRD size distribution')
    save_fig(fig, 'peak_size_distribution')
    
    stats['median_peaksize'] = crd_sizes.median()
    stats['max_peaksize'] = crd_sizes.max()

    
def analyze_nucleotide_coverage(crd_df, neg_df, stats, genome_size, fig_saver):
    crd_coverages = compute_nuc_coverage(crd_df)
    crd_coverages.name = 'CRDs'
    stats['CRD_coverage_sum_bp'] = crd_coverages.sum()
    stats['CRD_coverage_sum_genomefraction'] = stats['CRD_coverage_sum_bp'] / genome_size
    stats['CRD_coverage_median_bp'] = crd_coverages.median()
    
    neg_coverages = compute_nuc_coverage(neg_df)
    stats['shuffled_coverage_sum_bp'] = neg_coverages.sum()
    stats['shuffled_coverage_median_bp'] = neg_coverages.median()
    neg_coverages.name = 'shuffled'
    
    if fig_saver:
        fig, ax = plt.subplots(1, figsize=(4,3))
        seaborn.kdeplot(numpy.log10(crd_coverages.astype(float)), ax=ax)
        seaborn.kdeplot(numpy.log10(neg_coverages.astype(float)), ax=ax)
        ax.set_ylabel('relative frequency')
        ax.set_xlabel('Log10 bp coverage')
        ax.set_title('Covered bp')
        save_fig(fig, 'coverage_distributions')
        fig.show()

        
def shuffle_crd_members(crds_df, all_peaks_df, random_seed=None, preserve_gap_sizes=False):
    """
    Given a dataframe of CCDs :param:`crds_df` with compound peak names as an index,
    and a dataframe of individual peak information :param:`all_peaks_df`,
    return a dictionary, keyed by chromosome, of lists of individual peak IDs representing 
    CCDs formed by randomly joining adjacent peaks into contiguous blocks with the same size 
    distribution as the original CCDs (per chromosome).
    """
    numpy.random.seed(random_seed)
    # Just wrap "shuffle_regions" using peak indices instead of genomic coordinates
    region_tuples = []
    chrom_sizes = {}
    peak_nums_to_names_by_chrom = {}
    
    for this_chrom, this_chrom_crds in crds_df.groupby('chrom'):
        this_chrom_crds = this_chrom_crds.sort_values('chromStart')
#         print(this_chrom, this_chrom_crds)
        these_annotated_peaks = all_peaks_df.loc[all_peaks_df['Chr'] == this_chrom].sort_values('Start')
        these_peaknames = list(these_annotated_peaks.index)
        
        chrom_sizes[this_chrom] = len(these_peaknames)
        
        peak_name_to_number = {peak_name:num for num, peak_name in enumerate(these_peaknames)}
        peak_nums_to_names_by_chrom[this_chrom] = toolbox.invert_dict(peak_name_to_number)
    
        for crd_name in this_chrom_crds.index:
            member_peaks = crd_name.split('_')
            member_numbers = [peak_name_to_number[peak_name] for peak_name in member_peaks]
            region_tuples.append(('', this_chrom, min(member_numbers), max(member_numbers)))
            

    shuffled_region_tuples = shuffle_regions(region_tuples=region_tuples, 
                                             chrom_lengths=chrom_sizes,
                                             preserve_gap_sizes=False)
    
    shuffled_crd_members_by_chrom = {}
    for chrom in shuffled_region_tuples:
        shuffled_crd_members_by_chrom[chrom] = []
        for name, first_peak, end_peak in shuffled_region_tuples[chrom]:
            member_peaknames = [peak_nums_to_names_by_chrom[chrom][this_peak] for this_peak in range(first_peak, end_peak)]
            shuffled_crd_members_by_chrom[chrom].append(member_peaknames)
        
    return construct_crd_bed_from_member_peaks(shuffled_crd_members_by_chrom, all_peaks_df)
        
def construct_crd_bed_from_member_peaks(member_peaks_by_chrom, annotated_peaks_df):
    """
    Given a dictionary (keyed by chromosome) of lists of lists of peak names that make
    up each CCD, return a dataframe of CCDs made by merging those peaks together
    """
    data_dict = {}
    for chrom in member_peaks_by_chrom:
        for peakset in member_peaks_by_chrom[chrom]:
            compound_name = '_'.join(peakset)
            start_peak = peakset[0]
            end_peak = peakset[-1]
            this_data = {'chrom':chrom, 
                         'chromStart':annotated_peaks_df.loc[start_peak, 'Start'],
                         'chromEnd':annotated_peaks_df.loc[end_peak, 'End'],
                         'strand':'+',
                         'score':0}
            data_dict[compound_name] = this_data
    data_df = pandas.DataFrame(data_dict).T
    data_df = data_df.loc[:, ['chrom', 'chromStart', 'chromEnd', 'strand', 'score']]
    data_df.index.name = 'name'
    return data_df

def compute_nuc_coverage(region_df):
    return region_df.loc[:,'chromEnd'] - region_df.loc[:,'chromStart']


def generate_chrom_vector_from_pcs(chrom_name, pc_series, ws=0, we=0, genome=DEFAULT_GENOME):
    """
    Creates an interpolated basepair-resolution vector from a pandas.Series of Hi-C PC1 data.
    """
    if we == 0:
        we = genome.contig_lengths[chrom_name]
    assert we > ws
    assert we <= genome.contig_lengths[chrom_name]
        
    chrom_vector = numpy.zeros(shape=we-ws, dtype=float)
    for pos_idx in range(pc_series.shape[0]):
        this_pos, this_val = pc_series.index[pos_idx], pc_series.values[pos_idx]
        if ws <= this_pos < we:
            if pos_idx < pc_series.shape[0] -1:
                next_pos, next_val = pc_series.index[pos_idx+1], pc_series.values[pos_idx+1]

            chrom_vector[this_pos:next_pos] = numpy.linspace(this_val, next_val, num=next_pos-this_pos)
    return chrom_vector

    
def analyze_crds(peakset_name, annotated_peaks, tag_folders, genome=DEFAULT_GENOME,
                chroms_to_include=('chr18',)):
    stats = pandas.Series()
    output_dir = os.path.join(ANALYSIS_BASEPATH, peakset_name)
    os.makedirs(output_dir, exist_ok=True)
    
    save_fig = myplots.generate_figure_saver(output_dir)
    
    # Load full CRD set
    crd_bed_fname = os.path.join(PROJECT_PATH, '{}.homer'.format(peakset_name))
    crd_df = pandas.read_csv(crd_bed_fname, index_col=0, sep='\t')
    
    # Filter to specified chroms
    log_print('Filtering to chromosomes {}'.format(chroms_to_include))
    crd_df = crd_df.loc[numpy.in1d(crd_df.chrom, chroms_to_include)]
    
    # Save filtered data to output dir
    crd_bed_fname = os.path.join(output_dir, '{}.homer'.format(peakset_name))
    crd_df[['chrom', 'chromStart', 'chromEnd', 'strand', 'score']].to_csv(crd_bed_fname, sep='\t')
    
    neg_bed_fname = os.path.join(output_dir, '{}_shuffled.homer'.format(peakset_name))
    neg_df = shuffle_crd_members(crd_df, atac_notx_anno,1)
    neg_df.to_csv(neg_bed_fname, sep='\t')    
    
    # Peak sizes
    log_print('Computing peak sizes ...')
    analyze_peak_sizes(crd_df=crd_df, stats=stats, genome=genome, fig_saver=save_fig)
    # Nucleotide coverage
    log_print('Computing nucleotide coverage ...')
    analyze_nucleotide_coverage(crd_df, neg_df, stats=stats, genome_size=genome.size, fig_saver=save_fig)
        
    pos_matrix_fname = os.path.join(output_dir, '{}_tag_matrix.tsv'.format(peakset_name))
    try:
        pos_tag_matrix = pandas.read_csv(pos_matrix_fname,
                                 index_col=0,
                                 sep='\t')
    except FILE_EXCEPTIONS:
        log_print('Annotating CRDs with tag folders {} ...'.format(', '.join(tag_folders.values())))
        count_tags_in_regions(crd_bed_fname, pos_matrix_fname, [tag_folders[key] for key in sorted(tag_folders)], cmdlineonly=False)

        pos_tag_matrix = pandas.read_csv(pos_matrix_fname,
                                     index_col=0,
                                     sep='\t')
        
        log_print('Done',2)
    else:
        log_print('Found pre-generated CRD tag matrix at {}'.format(pos_matrix_fname))
        
    pos_tag_matrix.index.name = 'PeakID'
    pos_sizes = pos_tag_matrix['End'] - pos_tag_matrix['Start'] 
    clean_pos_matrix = pos_tag_matrix.iloc[:, 18:]
    clean_pos_matrix.columns = clean_col_names(clean_pos_matrix.columns, name_mapper=toolbox.invert_dict(tag_folders))
    clean_pos_matrix = clean_pos_matrix.divide(pos_sizes/10000, axis=0) # normalize by basepair size of crd -- units of tags per 10 kbp
    log_crd_matrix = numpy.log2(clean_pos_matrix + 1)        
    
    neg_matrix_fname = os.path.join(output_dir, '{}_shuffled_tag_matrix.tsv'.format(peakset_name))
    try:
        neg_tag_matrix = pandas.read_csv(neg_matrix_fname,
                                     index_col=0,
                                     sep='\t')
    except FILE_EXCEPTIONS:       
        log_print('Annotating shuffled control with tag folders {} ...'.format(', '.join(tag_folders.values())))
        count_tags_in_regions(neg_bed_fname, neg_matrix_fname, sorted(tag_folders.values()), cmdlineonly=False) 
        neg_tag_matrix = pandas.read_csv(neg_matrix_fname,
                                         index_col=0,
                                         sep='\t')
        log_print('Done',2)

    else:
        log_print('Found pre-generated shuffled tag matrix at {}'.format(neg_matrix_fname))

    neg_tag_matrix.index.name = 'PeakID'
    neg_sizes = neg_tag_matrix['End'] - neg_tag_matrix['Start'] 
    clean_neg_matrix = neg_tag_matrix.iloc[:, 18:]
    clean_neg_matrix.columns = clean_col_names(clean_neg_matrix.columns, name_mapper=toolbox.invert_dict(tag_folders))
    clean_neg_matrix = clean_neg_matrix.divide(neg_sizes/10000, axis=0) # normalize by basepair size of crd -- units of tags per 10 kbp
    log_neg_matrix = numpy.log2(clean_neg_matrix + 1)               
        
    log_print('Computing enrichments ...')
    for col in clean_pos_matrix.columns:
        pos_series = log_crd_matrix[col]
        pos_series.name = 'CRDs'
        
        neg_series = numpy.log2(clean_neg_matrix.loc[:,col] + 1)
        neg_series.name = 'neg_control'
        neg_series = log_neg_matrix[col]
        
        pos_mean, pos_std = pos_series.mean(), pos_series.std()
        neg_mean, neg_std = neg_series.mean(), neg_series.std()
        enrichment = pos_mean - neg_mean
        _, pval=scipy.stats.ttest_ind_from_stats(pos_mean, pos_std, len(pos_series), neg_mean, neg_std, len(neg_series))

        stats['{}_enrichment'.format(col)] = enrichment
        stats['{}_pval'.format(col)] = pval
        
    print(stats)
    
    pos_melted = pandas.melt(log_crd_matrix, var_name='Data type', value_name='Log2 tags per 10 kbp')
    pos_melted['dataset'] = 'CRDs'
    neg_melted = pandas.melt(log_neg_matrix, var_name='Data type', value_name='Log2 tags per 10 kbp')
    neg_melted['dataset'] = 'shuffled'

    merged_melted = pandas.concat((pos_melted, neg_melted), axis=0)
    
    fig, ax = plt.subplots(1, figsize=(1 * len(tag_folders),4))
    seaborn.violinplot(x='Data type', y='Log2 tags per 10 kbp', hue='dataset', data=merged_melted, ax=ax)
    save_fig(fig, 'tagcount_violinplot')
#     for fig_ext in myplots.FIG_EXTS:
#         plt.savefig(os.path.join(output_dir, 
#                                  'tagcount_violin_plot.{}'.format(fig_ext)),
#                                  dpi=myplots.FIG_DPI, bbox_inches='tight')
    
    return stats
    
    
def generate_placseq_tuple_dicts(placseq_df):
    starts = {}
    for chrom, chrom_data in placseq_df.groupby('chr1'):
        this_starts = set(zip(chrom_data.start1, chrom_data.start1, chrom_data.end1))
        this_starts.update(zip(chrom_data.start2, chrom_data.start2, chrom_data.end2))
        starts[chrom] = sorted(this_starts)
    return starts

def generate_crd_tuple_dicts(crd_df):
    return {chrom:list(zip(chrom_data.index, chrom_data.chromStart, chrom_data.chromEnd)) for chrom, chrom_data in crd_df.groupby('chrom')}


def convert_overlap_tuples_to_dicts(overlap_tuples):
    """
    Takes the output of an intervaloverlaps run (list of tuples containing names of overlapping regions)
    and returns a pair of dictionaries giving the interaction partners of each side.
    """
    a_dict = {}
    b_dict = {}
    for a_name, b_name, _, _ in overlap_tuples:
        if a_name not in a_dict:
            a_dict[a_name] = set([])
        a_dict[a_name].add(b_name)
        
        if b_name not in b_dict:
            b_dict[b_name] = set([])
        b_dict[b_name].add(a_name)
    
    return a_dict, b_dict


def generate_placseq_interaction_dict(placseq_df):
    """
    Given a DataFrame of placseq interactions, return a nested dictionary,
    keyed by chromosome, then by bin start point, of all the bin starts that interact
    with it.
    """
    interactions_by_chrom = {}
    for chrom, chrom_data in placseq_df.groupby('chr1'):        
#         interactions_by_chrom[chrom] = {}
        this_chrom_dict = {}
        for interaction_id in chrom_data.index:            
            left_start = chrom_data.loc[interaction_id, 'start1']
            right_start = chrom_data.loc[interaction_id, 'start2']
            if left_start not in this_chrom_dict:
                this_chrom_dict[left_start] = set([right_start])
            else:
                this_chrom_dict[left_start].add(right_start)
            
            if right_start not in this_chrom_dict:
                this_chrom_dict[right_start] = set([left_start])
            else:
                this_chrom_dict[right_start].add(left_start)
                
        interactions_by_chrom[chrom] = this_chrom_dict
    return interactions_by_chrom

def annotate_interaction_region_overlaps(interaction_tuples, interaction_dict, interactions_to_regions):
    """
    For a single chromosome, use pre-processed data to annotate overlaps between plac-seq interactions
    and CRDs in 4 categories:
    
    1. Number of interactions where both endpoints overlap the same region.
    2. Number of interactions where both endpoints overlap different regions.
    3. Number of interactions where one endpoint overlaps a region and the other does not overlap any region.
    4. Number of interactions where neither endpoint overlaps any region.
    
    :param:`interaction_dict`    an adjacency dictionary of interaction endpoints. Keys are window start positions and
        values are sets of start positions of windows that are linked by plac-seq interactions.
    :param:`interactions_to_region`    a dictionary keyed by window start positions containing sets of identifiers for
        CRDs that overlap that window.
    """
    interaction_status = {}

    for tup in interaction_tuples:
        interaction_a = tup[1]
        interaction_partners = [interaction_b for interaction_b in interaction_dict[interaction_a] if interaction_b > interaction_a]
        for interaction_b in interaction_partners:
            if interaction_a in interactions_to_regions:
                interaction_a_overlapping_regions = interactions_to_regions[interaction_a]
            else:
                interaction_a_overlapping_regions = {}
            if interaction_b in interactions_to_regions:
                interaction_b_overlapping_regions = interactions_to_regions[interaction_b]
            else:
                interaction_b_overlapping_regions = {}           
                
            if len(interaction_a_overlapping_regions) or len(interaction_b_overlapping_regions):            
                # we now have eliminated cat 4
                if len(interaction_a_overlapping_regions) and len(interaction_b_overlapping_regions):            
                    # we have now eliminated cat 3
                    # and just need to check if a and b share a CRD
                    if len(interaction_a_overlapping_regions.intersection(interaction_b_overlapping_regions)) > 0:
                        # it's cat 1
                        interaction_status[(interaction_a, interaction_b)] = 1
                    else:
                        interaction_status[(interaction_a, interaction_b)] = 2
                else:
                    interaction_status[(interaction_a, interaction_b)] = 3
            else:
                interaction_status[(interaction_a, interaction_b)] = 4

    return interaction_status

def count_interaction_region_overlaps_one_chrom(interaction_tuples, interaction_dict, interactions_to_regions):
    """
    For a single chromosome, use pre-processed data to count overlaps between plac-seq interactions
    and CRDs in 4 categories:
    
    1. Number of interactions where both endpoints overlap the same region.
    2. Number of interactions where both endpoints overlap different regions.
    3. Number of interactions where one endpoint overlaps a region and the other does not overlap any region.
    4. Number of interactions where neither endpoint overlaps any region.
    
    :param:`interaction_dict`    an adjacency dictionary of interaction endpoints. Keys are window start positions and
        values are sets of start positions of windows that are linked by plac-seq interactions.
    :param:`interactions_to_region`    a dictionary keyed by window start positions containing sets of identifiers for
        CRDs that overlap that window.
    """
    overall_count = 0
    cat1_count, cat2_count, cat3_count, cat4_count = 0, 0, 0, 0

    for tup in interaction_tuples:
        interaction_a = tup[1]
        interaction_partners = interaction_dict[interaction_a]
        for interaction_b in interaction_partners: # this loop will execute once per interaction (unique pair of interacting start points)
            overall_count += 1
            if interaction_a in interactions_to_regions:
                # we now have eliminated cat 4
                if interaction_b in interactions_to_regions:
                    # we have now eliminated cat 3
                    # and just need to check if a and b share a CRD
                    if len(interactions_to_regions[interaction_a].intersection(interactions_to_regions[interaction_b])) > 0:
                        # it's cat 1
                        cat1_count += 1
#                         print(interaction_a, interaction_b, interactions_to_regions[interaction_a].intersection(interactions_to_regions[interaction_b]))
                    else:
                        cat2_count += 1
                else:
                    cat3_count += 1
            else:
                cat4_count += 1
    assert overall_count == cat1_count + cat2_count + cat3_count + cat4_count # sanity check
    return cat1_count, cat2_count, cat3_count, cat4_count
    
def count_interaction_region_overlaps(placseq_df, crd_df):
    overall_count = 0
    cat1_count, cat2_count, cat3_count, cat4_count = 0, 0, 0, 0
    
    placseq_edges = generate_placseq_interaction_dict(placseq_df)
    placseq_regions = generate_placseq_tuple_dicts(placseq_df)
    crd_regions = generate_crd_tuple_dicts(crd_df)
    
    for chrom in placseq_regions:
        if chrom in crd_regions:
            overlaps = genomicoverlaps2.compute_interval_overlaps(placseq_regions[chrom], crd_regions[chrom])
            this_chrom_placseq_to_regions, this_chrom_regions_to_placseq = convert_overlap_tuples_to_dicts(overlaps)
            
            c1, c2, c3, c4 = count_interaction_region_overlaps_one_chrom(placseq_regions[chrom], placseq_edges[chrom], this_chrom_placseq_to_regions)
            print(chrom, c1, c2, c3, c4)
            cat1_count += c1
            cat2_count += c2
            cat3_count += c3
            cat4_count += c4
            
    return {'intra_region_interactions':cat1_count, 'inter_region_interactions':cat2_count, 'xenic_interactions':cat3_count, 'nonoverlapping_interactions':cat4_count }


def annotate_interaction_region_overlaps(interaction_tuples, interaction_dict, interactions_to_regions):
    """
    :param:`interaction_dict`    an adjacency dictionary of interaction endpoints. Keys are window start positions and
        values are sets of start positions of windows that are linked by plac-seq interactions.
    :param:`interactions_to_region`    a dictionary keyed by window start positions containing sets of identifiers for
        CRDs that overlap that window.
    """
    interaction_status = {}

    for tup in interaction_tuples:
        interaction_a = tup[1]
        interaction_partners = [interaction_b for interaction_b in interaction_dict[interaction_a] if interaction_b > interaction_a]
        for interaction_b in interaction_partners:
            if interaction_a in interactions_to_regions:
                interaction_a_overlapping_regions = interactions_to_regions[interaction_a]
            else:
                interaction_a_overlapping_regions = {}
            if interaction_b in interactions_to_regions:
                interaction_b_overlapping_regions = interactions_to_regions[interaction_b]
            else:
                interaction_b_overlapping_regions = {}           
                
            if len(interaction_a_overlapping_regions) or len(interaction_b_overlapping_regions):            
                # we now have eliminated cat 4
                if len(interaction_a_overlapping_regions) and len(interaction_b_overlapping_regions):            
                    # we have now eliminated cat 3
                    # and just need to check if a and b share a CRD
                    if len(interaction_a_overlapping_regions.intersection(interaction_b_overlapping_regions)) > 0:
                        # it's cat 1
                        interaction_status[(interaction_a, interaction_b)] = 1
                    else:
                        interaction_status[(interaction_a, interaction_b)] = 2
                else:
                    interaction_status[(interaction_a, interaction_b)] = 3
            else:
                interaction_status[(interaction_a, interaction_b)] = 4

    return interaction_status

def construct_interaction_patch(xlim, ylim,
                                left_bin_center, right_bin_center, 
                                color='k',
                                baseline=0.0,
                                vertical_scaling_factor=1,
                                direction='down'):
    arc_width = right_bin_center - left_bin_center
    
    arc_height = (ylim[1] - ylim[0]) * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor # at a vetical_scaling_factor of 1 an arc that oompletely spans the plot horizontally will also span it vertically.
    xy = ((left_bin_center + right_bin_center)/2, baseline)
    
    if direction == 'down':
        theta1, theta2 = 180, 0
    else:
        theta1, theta2 = 0, 180
        
#     theta1, theta2 = 0, 360
        
#     print('thetas: {} {}'.format(theta1, theta2))
    
#     print('constructing arc with xy {}, height {}, width {}'.format(xy, arc_height, arc_width))
    return matplotlib.patches.Arc(xy=xy, height=arc_height, width=arc_width,
                                  theta1=theta1, theta2=theta2, color=color)
def draw_square_interaction(ax,
                            left_bin_center, 
                            right_bin_center, 
                            color='k',
                            baseline=0.0,
                            vertical_scaling_factor=1,
                            direction='down'):
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    arc_width = right_bin_center - left_bin_center
    
    arc_height = (ylim[1] - ylim[0]) * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor # at a vetical_scaling_factor of 1 an arc that oompletely spans the plot horizontally will also span it vertically.
    xy = ((left_bin_center + right_bin_center)/2, baseline)
    
    if direction == 'down':
        y_offset = -arc_height
    else:
        y_offset = arc_height

    ax.plot((left_bin_center, left_bin_center, right_bin_center, right_bin_center), (baseline, y_offset, y_offset, baseline), color=color)

    
def draw_triangle_interaction(ax, xlim, ylim,
                            left_bin_center, 
                            right_bin_center, 
                            color='k',
                            baseline=0.0,
                            vertical_scaling_factor=1,
                            direction='down'):
    
    arc_width = right_bin_center - left_bin_center
    
    arc_height = (ylim[1] - ylim[0]) * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor # at a vetical_scaling_factor of 1 an arc that oompletely spans the plot horizontally will also span it vertically.
    
    if direction == 'down':
        y_offset = -arc_height
    else:
        y_offset = arc_height
        
#     print('Drawing triangular interaction between {} and {}, height {}'.format(left_bin_center, right_bin_center, arc_height))
#     print('\t', ylim[1] - ylim[0], arc_width, xlim[1] - xlim[0], vertical_scaling_factor)
    ax.plot((left_bin_center, (left_bin_center + right_bin_center) / 2, right_bin_center), (baseline, y_offset, baseline), color=color)    
    
    
def compute_half_arc_points(center, a, b, theta1, theta2, direction='down', num_points=NUM_ARC_POINTS):
    """
    Computes the coordinates for component points of a polygonal approximation to
    an ellipse for a single quadrant.
    """
    # ToDo: Add input validation to make sure we stay within a single quadrant.
    all_points = []
    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)
#     print('coordinates for half-arc. center: {}, a: {}, b: {}'.format(center,a,b))
        
    for i in range(0, num_points):
        theta = (theta2 - theta1) * (i / max(num_points - 1,1)) + theta1    
        fi = numpy.pi / 2 - numpy.arctan(numpy.tan(theta))
        x = center[0] + a * numpy.cos(fi)
        y = center[1] + b * numpy.sin(fi)
        x_coords[i] = x
        y_coords[i] = y
#         print(i, theta, fi, x, y)

    return x_coords, y_coords


def draw_arc(ax, center, height, width, theta1=0, theta2=numpy.pi, color='k', direction='down', num_points=NUM_ARC_POINTS):
    """
    Since Matplotlib's Arc Patches are broken at the moment, we draw arcs using the .plot() method
    instead.
    """
    LEFT_END_THETA = numpy.pi / 2
    RIGHT_END_THETA = numpy.pi * 1.5
    MIDPOINT_THETA = numpy.pi
    
    vertical_baseline=center[1]
    
    assert LEFT_END_THETA <= theta1 <= theta2 <= RIGHT_END_THETA
    
    b = height
    a = width / 2

    all_points = []
    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)

    # determine how to allocate points
#     print('overall thetas: {}-{}'.format(theta1, theta2))
    left_angle_span =  min(max(MIDPOINT_THETA - theta1, 0), theta2 - theta1)
    right_angle_span = min(max(theta2 - MIDPOINT_THETA, 0), theta2 - theta1)
    total_angle_span = left_angle_span + right_angle_span
#     print('left theta span: {}, right theta span: {}, total span: {}'.format(left_angle_span, right_angle_span, total_angle_span))
    left_points = int(num_points * left_angle_span / total_angle_span)
    right_points = num_points - left_points
#     print(left_points, right_points)
    
    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)
    
    if left_points:
        # plot upper left quadrant
        left_theta2 = theta1 + left_angle_span
#         print('Drawing left arc from {} to {}'.format(theta1, left_theta2))
        x,y = compute_half_arc_points(center=(center[0], 0), 
                      a=a, b=b, 
                      theta1=theta1, theta2=left_theta2,
                      num_points=left_points)
        x_coords[:left_points] = x[:]
        y_coords[:left_points] = y[:]
    if right_points:
        # plot upper right quadrant
        right_theta1 = theta2 - right_angle_span
#         print('Drawing right arc from {} to {}'.format(right_theta1, theta2))
        x,y = compute_half_arc_points(center=(center[0], 0), 
                      a=a, b=b, 
                      theta1=right_theta1, theta2=theta2,
                      num_points=right_points)
        x_coords[left_points:] = x[:]
        y_coords[left_points:] = y[:]
        
    if direction == 'down':
        y_coords = - y_coords
        
#     print('y-coords range: {}, {}'.format(y_coords.min(), y_coords.max()))   
#     print('addeding baseline {}'.format(vertical_baseline))

    y_coords += vertical_baseline
          
#     print('y-coords range: {}, {}'.format(y_coords.min(), y_coords.max()))

    ax.plot(x_coords, y_coords, color=color) 
    
    
def draw_visible_arc(ax, center, height, width, ws, we,
                     color='k', direction='down', num_points=NUM_ARC_POINTS):
    """
    Draws a 180 degree arc truncated by an interval of x coordinates.
    Does not truncate based on y coordinates
    """
    LEFT_END_THETA = numpy.pi / 2 + 0.00001
    RIGHT_END_THETA = numpy.pi * 1.5 - 0.00001
    MIDPOINT_THETA = numpy.pi
    
    def infer_theta_cutoff(x, center, width, height):
        a = width / 2
        fi = numpy.arccos((x - center[0]) / a)
        theta = numpy.arctan(1/ numpy.tan(fi)) + numpy.pi
        return theta
    
    if ws > center[0] - width / 2 :
        theta_start = infer_theta_cutoff(x=ws, center=center, width=width, height=height)
#         print('\tTruncating left. Original arc left: {}, new arc left: {}, theta: {}, old_left_theta: {}'.format(center[0] - width / 2, ws, theta_start, LEFT_END_THETA))
    else:
        theta_start = LEFT_END_THETA
    
    if we < center[0] + width / 2:
#         print('truncating right')
        theta_end = infer_theta_cutoff(x=we, center=center, width=width, height=height)
    else:
        theta_end = RIGHT_END_THETA 
    
#     print(theta_start, theta_end)
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
                         num_points=NUM_ARC_POINTS):
    """
    """
    arc_width = right_bin_center - left_bin_center
    horizontal_span = xlim[1] - xlim[0]
    
    if direction=='down':
        vertical_span = baseline - ylim[0]
    else:
        vertical_span = ylim[1] - baseline

#     vertical_span = 1    
#     print('target axes has span ({},{})'.format(horizontal_span, vertical_span))

    arc_height = vertical_span * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor
#     print('Drawing interaction (arc) between {} and {}, height {}'.format(left_bin_center, right_bin_center, arc_height))

    draw_visible_arc(ax=ax, 
                     center=((left_bin_center + right_bin_center) / 2, baseline),
                     height=arc_height, 
                     ws=xlim[0],
                     we=xlim[1],
                     width=arc_width, 
                     color=color, 
                     direction=direction,
                     num_points=num_points)
    
    
def generate_rotation_matrix(angle):
    angle = numpy.deg2rad(angle)
    return numpy.array([[numpy.cos(angle), -numpy.sin(angle)],
                        [numpy.sin(angle), numpy.cos(angle)]])


def rotate_vector(vec, angle):
    return generate_rotation_matrix(angle).dot(vec.T)

def scatter_heat(ax, data, marker='s', marker_size=20, cmap='YlOrRd', alpha=1):
    """
    Assuming :param:`data` is a DataFrame containing series of data in rows and whose
    column index corresponds to the x position of these data, scatters the data
    while coloring according :param:`cmap`
    """
#     print(data.shape)
    color_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(), cmap=cmap)
    for row_idx in range(data.shape[0]):
        x = data.columns
        y = numpy.full(data.shape[1], fill_value=row_idx)
#         print(len(x), len(y))
        ax.scatter(data.columns, numpy.full(data.shape[1], fill_value=row_idx),
                   s=marker_size,
                   marker=marker,
                   c=color_mapper.to_rgba(data.iloc[row_idx],
                                          alpha=alpha))
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index)

    
def diagonal_triangle_heat(ax, data, annotated_peaks, 
                          cmap='RdBu_r',
                          max_distance=100000,
                          alpha=1,
                          linewidth=3,
                          vmin=-1, vmax=1,
                          marker_s=15,    
                          scatter_kwargs={}):
    """
    Assuming :param:`data` is a DataFrame containing a square of matrix of pairwise
    values, scatter the upper triangle of these data in a diagonal/triangular format
    while coloring according :param:`cmap`.
    """
    assert data.shape[0] == data.shape[1]
    n = data.shape[0]
    starts = []
    mids = []
    ends = []
    ys = []
    values = []

    color_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
    for row_idx, row_name in enumerate(data.index):
        row_peak_midpoint = (annotated_peaks.loc[row_name, 'End'] + annotated_peaks.loc[row_name, 'Start']) / 2
        for col_idx, col_name in enumerate(list(data.columns)[row_idx+1:]):
            col_peak_midpoint = (annotated_peaks.loc[col_name, 'End'] + annotated_peaks.loc[col_name, 'Start']) / 2
            distance = col_peak_midpoint - row_peak_midpoint

            if max_distance == 0 or distance <= max_distance:
                intersection_x = (col_peak_midpoint + row_peak_midpoint) / 2
                intersection_y = distance
                starts.append(col_peak_midpoint)
                mids.append(intersection_x)
                ends.append(row_peak_midpoint)
                ys.append(intersection_y)
                normed_value = (data.loc[row_name, col_name] - vmin) / (vmax - vmin) 
                values.append(normed_value)
    
    if linewidth > 0:
        for start, mid, end, intersection_y in zip(starts, mids, ends, ys):
            ax.plot((start, mid, end),
                    (0, intersection_y, 0),
                    c='k', linewidth=0.5, zorder=1)

    sort_order = numpy.argsort(values)
    
    ax.scatter(numpy.array(mids)[sort_order], numpy.array(ys)[sort_order],
            c=color_mapper.to_rgba(numpy.array(values)[sort_order], alpha=alpha),
               alpha=alpha, linewidth=0,
           marker='D', s=marker_s, zorder=2, **scatter_kwargs)

    
class CRDSubPlot():
    def set_globals(self, chrom, ws, we, fig_width=64, row_height=4):
        self.chrom = chrom
        self.ws = ws
        self.we = we
        self.fig_width = fig_width
        self.row_height = row_height
        self.aspect_ratio = fig_width / row_height
        


class PlacseqPlot(CRDSubPlot):
    def __init__(self, placseq_df,
                 placseq_bin_size, 
                 crd_df, 
                 arc_color='k',
                 direction='up',
                 baseline=0,
                 vertical_scaling_factor=1,
                 thickness_column=None,
                 show_bin_centers=True,
                 label='Plac-seq'):
         
        self.placseq_df = placseq_df
        self.placseq_bin_size = placseq_bin_size
        self.baseline = baseline
        self.vertical_scaling_factor = vertical_scaling_factor
        self.arc_color = arc_color
        self.label = label
        self.direction = direction
        self.show_bin_centers = show_bin_centers
        self.thickness_column = thickness_column
     
    def plot(self, ax):
        # Filter the interacton DataFrame to eliminate anything without one anchor point within the visible window.
        visible_interactions = self.placseq_df.loc[self.placseq_df['chr1'] == self.chrom]
        left_bin_midpoints = (visible_interactions['end1'] + visible_interactions['start1']) / 2
        right_bin_midpoints = (visible_interactions['end2'] + visible_interactions['start2']) / 2
        left_visible = (left_bin_midpoints >= self.ws) & (left_bin_midpoints <= self.we)
        right_visible = (right_bin_midpoints >= self.ws) & (right_bin_midpoints <= self.we)
        visible_interactions = visible_interactions.loc[left_visible | right_visible]
        
        original_ylim = ax.get_ylim()
#         print('original ylim: {}'.format(original_ylim))
        ax.set_xlim(self.ws, self.we)
#         if direction == 'down':
#             ax.set_ylim(
        
        for interaction_id in visible_interactions.index:
#             print('Drawing interaction from {} to {}'.format(left_bin_midpoints.loc[interaction_id],
#                                                      right_bin_midpoints.loc[interaction_id]))
#             ax.add_patch(construct_interaction_patch(left_bin_center=left_bin_midpoints.loc[interaction_id],
#                                                      right_bin_center=right_bin_midpoints.loc[interaction_id],
#                                                      xlim=(self.ws, self.we),                                                                            
#                                                      ylim=ax.get_ylim(),
#                                                      vertical_scaling_factor=self.vertical_scaling_factor,
#                                                      color=self.arc_color, 
#                                                      baseline=self.baseline,
#                                                      direction=self.direction))
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
            
            leftmost_tick = numpy.ceil((self.ws - self.placseq_bin_size/2) / self.placseq_bin_size) * self.placseq_bin_size + self.placseq_bin_size/2
            rightmost_tick = numpy.floor((self.we - self.placseq_bin_size/2) / self.placseq_bin_size + 1) * self.placseq_bin_size + self.placseq_bin_size/2

            ax.set_xticks(numpy.arange(leftmost_tick, rightmost_tick, self.placseq_bin_size))
            ax.set_xticklabels([])
            
            if self.direction == 'down':
                ax.xaxis.set_ticks_position('top')

        
class DiagScatter(CRDSubPlot):
    def __init__(self, features_df, annotated_features_df, corr_power=2,
                 max_triangle_distance='auto', diagonal_triangle_heat_kwargs={'cmap':'Greys', 'alpha':0.01}):
#         self.features_df = features_df
        features_df = features_df.T
        self.triangle_data = features_df.corr()**corr_power
        self.annotated_features_df = annotated_features_df

        self.max_triangle_distance = max_triangle_distance
        self.diagonal_triangle_heat_kwargs = diagonal_triangle_heat_kwargs
    
    def plot(self, ax):
        if self.max_triangle_distance == 'auto':
            self.max_triangle_distance = (self.we-self.ws) / self.aspect_ratio
        diagonal_triangle_heat(ax=ax, 
                               data=self.triangle_data,
                               annotated_peaks=self.annotated_features_df,
                               max_distance=self.max_triangle_distance,
                               linewidth=0, **self.diagonal_triangle_heat_kwargs)

        if self.max_triangle_distance > 0:
            ax.set_ylim(0, self.max_triangle_distance)        
        ax.set_yticks(numpy.arange(toolbox.roundto(self.ws, 5000), toolbox.roundto(self.we, 5000), 5000))

        
class CRDLines(CRDSubPlot):
    def __init__(self, named_crds, crd_cmap='Purples', crd_outline_color=None, linewidth=5):
        """
        Takes an iterable of tuples in the form:

        (name, DataFrame)
        
        to be plotted in order using the .plot() method.
        """
#         try: # flexible input of either single DataFrame or iterable of DataFrames
#             assert len(crd_dfs[0].shape) == 2
#         except (AssertionError, AttributeError):
#             crd_dfs = [crd_dfs]
                
        self.named_crds = named_crds
        min_val = min([crd_df['score'].min() for crdset_name, crd_df in named_crds])
        max_val = max([crd_df['score'].max() for crdset_name, crd_df in named_crds])

        self.crd_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min_val,
                                                                                        vmax=max_val), 
                                                       cmap=crd_cmap)
        self.crd_outline_color = crd_outline_color
        self.linewidth = 5

    def plot(self, ax):
        for vertical_offset, (crdset_name, crd_df) in enumerate(self.named_crds):
            for crd_name in crd_df.index:            
                start_loc = crd_df.loc[crd_name, 'chromStart']
                end_loc = crd_df.loc[crd_name, 'chromEnd']
                if start_loc <= self.we or end_loc >= self.ws:
                    interval_color = self.crd_mapper.to_rgba(crd_df.loc[crd_name]['score'])
                    ax.plot((start_loc,
                             end_loc),
                             (vertical_offset,\
                              vertical_offset),
                            color=interval_color, 
                            linewidth=5)
#             print(vertical_offset)
        ax.set_yticks(numpy.arange(len(self.named_crds)))
        ax.set_yticklabels([crdset_name for crdset_name, crd_df in self.named_crds])
        ax.set_ylim(-1, len(self.named_crds))        
#         print(ax.get_ylim())


class CRDBoxes(CRDSubPlot):
    DEFAULT_PATCH_KWARGS={#'boxstyle':matplotlib.patches.BoxStyle('round', pad=0),
                          'linewidth':1,
                          'edgecolor':'k'}
    
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
        self.named_crds = named_crds[::-1]
        extent = max([numpy.abs(crd_df['score']).max() for crdset_name, crd_df in named_crds])

        self.crd_mapper = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-extent,
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
            self.height = vert_span / (num_crds * 2 -1)
            
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
            
            this_crds = crd_df.loc[(crd_df.chrom == self.chrom) & (((self.ws <= crd_df.chromStart) & (crd_df.chromStart <= self.we)) | ((self.ws <= crd_df.chromEnd) & (crd_df.chromEnd <= self.we)))]
#             print('looking at {}:{}-{}, found {} visible items out of {} in set {}'.format(self.chrom, self.ws, self.we,
#                 this_crds.shape[0], crd_df.shape[0], crdset_name))
            for crd_name in this_crds.index:            
                start_loc = this_crds.loc[crd_name, 'chromStart']
                end_loc = this_crds.loc[crd_name, 'chromEnd']
                interval_color = self.crd_mapper.to_rgba(crd_df.loc[crd_name]['score'])
#                     print('plotting {} in {}'.format(crd_name, interval_color))

                rec = matplotlib.patches.Rectangle(xy=(start_loc, bottom),
                                                        width=end_loc-start_loc,
                                                        height=self.height,
                                                        facecolor=interval_color,
                                                        **self.patch_kwargs)
                ax.add_patch(rec)
                if self.annotation_column and self.annotation_column in crd_df.columns:
#                         print((start_loc+end_loc)/2, midpoint)
                    ax.text(x=(start_loc+end_loc)/2, y=vert_midpoint, s='{:>0.2}'.format(crd_df.loc[crd_name, self.annotation_column]), ha='center')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
#         ax.set_ylim(-1, len(self.named_crds))  


class VectorOverlay(CRDSubPlot):
    def __init__(self, vector_series, label=None, color=None, center=False, scale=False, ylim=None, convolution_kernel=None):
        self.color = color
        # Check if vector has string labels -- if so, reindex by integer start positions
        if type(vector_series.index[0]) == str:
            vector_series = vector_series.copy()
            vector_series.index = [int(peak_id.split('-')[1]) for peak_id in vector_series.index]
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
            this_plot_vector = pandas.Series(scipy.signal.convolve(this_plot_vector, self.convolution_kernel, mode='same'), index=this_plot_vector.index)
        
        if self.scale:
            this_plot_vector /= (this_plot_vector.max() - this_plot_vector.min())
            this_plot_vector *= vert_span 
            
        if self.center:
            this_plot_vector -= this_plot_vector.mean()
            this_plot_vector += vert_center
#         print(self.center, self.scale, vert_center, vert_span, this_plot_vector.min(), this_plot_vector.max())
        this_plot_vector = this_plot_vector.loc[(this_plot_vector.index >= self.ws) & (this_plot_vector.index < self.we)]
        ax.plot(this_plot_vector.index, this_plot_vector, color=self.color, label=self.label)
        ax.set_ylim(ylim)


class FeatureStats(CRDSubPlot):
    def __init__(self, features_df, annotated_regions_df, 
                 plot_mean=True, plot_var=False, plot_sd=False, 
                 plot_relvar=False, plot_cv=True, znorm_tracks=False):
        self.displayed_features = features_df.T.copy()
        self.displayed_features.columns = (annotated_regions_df.loc[self.displayed_features.columns].Start + \
                                           annotated_regions_df.loc[self.displayed_features.columns].End) / 2
        
        self.znorm_tracks = znorm_tracks
        self.plot_mean = plot_mean
        self.plot_var = plot_var
        self.plot_sd = plot_sd
        self.plot_relvar = plot_relvar
        self.plot_cv = plot_cv
    
    def plot(self, ax):
        def plot_possibly_znormed(track_data):
            if self.znorm_tracks: track_data = toolbox.znorm(track_data)
            ax.plot(track_data)
        
        if self.plot_mean or self.plot_cv:
            feature_mean = self.displayed_features.mean(axis=0)
            if self.plot_mean:
                plot_possibly_znormed(feature_mean)
        if self.plot_var or self.plot_relvar:
            feature_var = self.displayed_features.var(axis=0)
            if self.plot_var:
                plot_possibly_znormed(feature_var)
        if self.plot_sd or self.plot_cv:
            feature_sd = self.displayed_features.std(axis=0)
            if self.plot_sd:
                plot_possibly_znormed(feature_sd)
        if self.plot_relvar:
            feature_relvar = feature_var / feature_mean
            plot_possibly_znormed(feature_relvar)                
        if self.plot_cv:
            feature_cv = feature_sd / feature_mean
            plot_possibly_znormed(feature_cv)
            
def arrange_genes(gene_data_list):
    """
    Given an iterable of gene data dictionaries,
    returns a list of lists of gene names to display at various levels
    """
    gene_data_list = sorted(gene_data_list, key=lambda x: x['end'] - x['start'], reverse=True)

    display_levels = [intervaltree.IntervalTree(),]

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
                found_home=True
                
    return [[gene_interval.data['ID'] for gene_interval in this_level] for this_level in display_levels]


class GeneModels(CRDSubPlot):
    def __init__(self, genome, label=None, color='k', 
                 feature_height=0.12, 
                 chevron_height=0.05, 
                 chevron_width=0.04, 
                 chevron_spacing=0.10,
                 truncation_size=0.10,
                 utr_endcap_width=0.04,
                gene_name_fontsize=8,
                valid_gene_types=GENE_TYPES,
                valid_component_types=COMPONENT_TYPES):
        self.color = color
        self.genome = genome
        self.label = label
        self.feature_height = feature_height # in inches
        self.chevron_height = chevron_height # in inches 
        self.chevron_width = chevron_width # in inches 
        self.chevron_spacing = chevron_spacing # in inches
        self.truncation_size = truncation_size # in inches
        self.utr_endcap_width = utr_endcap_width # in inches
        self.gene_name_fontsize = gene_name_fontsize
        
    def plot(self, ax):
        # find overlapping genes
        overlapping_genes = self.genome.genes.overlapping(self.chrom, self.ws, self.we)
        overlapping_components = self.genome.components.overlapping(self.chrom, self.ws, self.we)
        
        gene_display_levels = arrange_genes(overlapping_genes.values())
        ax.set_ylim((-0.5, len(gene_display_levels) -1 + 0.5))
               
        # convert inches to data coordinates
        chevron_spacing_dt = (self.we - self.ws) / (self.fig_width / self.chevron_spacing) 
        chevron_width_dt = (self.we - self.ws) / (self.fig_width / self.chevron_width) 
        truncation_width_dt = (self.we - self.ws) / (self.fig_width / self.truncation_size) 
        utr_endcap_width_dt = (self.we - self.ws) / (self.fig_width / self.utr_endcap_width) 
        
        feature_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.feature_height)
        chevron_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.chevron_height)
        feature_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.feature_height)
        truncation_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (self.row_height / self.truncation_size)
  
        for gene_num, level_genes in enumerate(gene_display_levels):
        
            # ToDo: make this universal. Divide the gene body into non-overlapping segments, each type of which has a template.
        
            for gene_id in level_genes:
                gene_data = overlapping_genes[gene_id]
#                 print(gene_id, gene_data['Name'])
                
                left_truncated = gene_data['start'] < self.ws
                right_truncated = gene_data['end'] > self.we
                
                visible_gene_start = max(gene_data['start'], self.ws)
                if left_truncated: visible_gene_start += truncation_width_dt * 2
                visible_gene_end = min(gene_data['end'], self.we)
                if right_truncated: visible_gene_end -= truncation_width_dt * 2

                ax.plot((visible_gene_start, visible_gene_end), (gene_num, gene_num), color=self.color)
                ax.text(x=(visible_gene_start + visible_gene_end)/2,
                        y=gene_num + feature_height_dt * 1.5, 
                        s=gene_data['Name'], 
                        ha='center',
#                         rotation=10,
                        fontsize=self.gene_name_fontsize)

                num_chevrons = int(max((visible_gene_end - visible_gene_start) / chevron_spacing_dt,1))
                chevron_remainder = (visible_gene_end - visible_gene_start) - (num_chevrons-1) * chevron_spacing_dt

                if gene_data['strand'] == '+':
                    chevron_x_delta = -chevron_width_dt
                else:
                    chevron_x_delta = chevron_width_dt
                
                for chevron_idx in range(num_chevrons):
#                     if num_chevrons =
                    chevron_x = visible_gene_start + chevron_idx * chevron_spacing_dt + chevron_remainder / 2
                    
                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num + chevron_height_dt), color=self.color)
                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num - chevron_height_dt), color=self.color)

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
                    if ((component_data['start'] >= visible_gene_start) and (component_data['start'] <= visible_gene_end)) or ((component_data['end'] >= visible_gene_start) and (component_data['end'] <= visible_gene_end)):
#                         visible_component_start = max(component_data['start'], visible_gene_start)
#                         visible_component_end = min(component_data['end'], visible_gene_end)
                        
                        # ToDo: systematize and condense the following
                        if component_data['type'] == 'five_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '+':
                                utr_body = matplotlib.patches.Rectangle(xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                                                        width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                                                    height=feature_height_dt,
                                                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                                                        width=utr_endcap_width_dt,
                                                                    height=feature_height_dt * 2,
                                                                    facecolor=self.color)

                            else:
                                utr_body = matplotlib.patches.Rectangle(xy=(component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                                                        width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                                                    height=feature_height_dt,
                                                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(xy=(component_data['start'], gene_num - feature_height_dt),
                                                                        width=utr_endcap_width_dt,
                                                                    height=feature_height_dt * 2,
                                                                    facecolor=self.color)                                

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)  
                            
                        elif component_data['type'] == 'three_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '-':
                                utr_body = matplotlib.patches.Rectangle(xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                                                        width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                                                    height=feature_height_dt,
                                                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                                                        width=utr_endcap_width_dt,
                                                                    height=feature_height_dt * 2,
                                                                    facecolor=self.color)

                            else:
                                utr_body = matplotlib.patches.Rectangle(xy=(component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                                                        width=max(component_data['end'] - component_data['start'] - self.utr_endcap_width, 0),
                                                                    height=feature_height_dt,
                                                                    facecolor=self.color)
                                utr_endcap = matplotlib.patches.Rectangle(xy=(component_data['start'], gene_num - feature_height_dt),
                                                                        width=utr_endcap_width_dt,
                                                                    height=feature_height_dt * 2,
                                                                    facecolor=self.color)                                
#                             print('{} utr body: xy: {}, width: {}, height: {}'.format(gene_data['Name'], (component_data['start'], gene_num - feature_height_dt / 2), component_data['end'] - component_data['start'] - utr_endcap_width_dt, feature_height_dt))
#                             print('utr endcap: xy: {}, width: {}, height: {}'.format((component_data['end'] - self.utr_endcap_width, gene_num - feature_height_dt), utr_endcap_width_dt, feature_height_dt * 2))

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)      
                        elif component_data['type'] == 'CDS':
                            cds = matplotlib.patches.Rectangle(xy=(component_data['start'], gene_num - feature_height_dt),
                                                                    width=component_data['end'] -  component_data['start'],
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
                  
    
def visualize_strain_data(plot_objects, 
                          chrom='', 
                          ws=0, we=0, 
                          genome=DEFAULT_GENOME,
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
        
    assert we > ws, 'Window end must be greater than window start! Got: {}, {}'.format(ws,we)
    try:
        if len(row_heights) == 1:
            row_heights = row_heights * len(plot_objects) # treat as a uniform row height
    except TypeError:
        row_heights = [row_heights] * len(plot_objects) # treat as a uniform row height
            
    assert len(row_heights) == len(plot_objects)
        
    span = we - ws
    xtick_increment = span / num_ticks
    rounding_increment = 5*10**numpy.round(numpy.log10(xtick_increment)-1)
    xtick_increment = toolbox.roundto(xtick_increment, rounding_increment)
    num_ticks = int(span / xtick_increment) + 1
#     print(span, xtick_increment, num_ticks)
    round_start = toolbox.roundto(ws, rounding_increment)

    seaborn.set_style(seaborn_style)
    
#     fig, axes = plt.subplots(len(plot_objects), figsize=(fig_width, numpy.sum(row_height)), sharex=True)
    fig = plt.figure(len(plot_objects), figsize=(fig_width, numpy.sum(row_heights) * (1 + ax_spacing * len(plot_objects))))
#     print(row_heights, ax_spacing)
    bottoms, heights = compute_ax_row_positions(row_heights=row_heights, ax_spacing=ax_spacing)

#     if len(plot_objects) == 1: axes = [axes]
    
    for ax_idx in range(len(plot_objects)):
#         print(ax_idx, len(plot_objects), bottoms, heights)
        this_ax = fig.add_axes([0, bottoms[ax_idx], 1, heights[ax_idx]])

        if ax_idx == len(plot_objects) - 1:
            this_ax.set_xticks(numpy.arange(num_ticks)*xtick_increment + round_start)
            this_ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
            this_ax.set_xlabel('{} position'.format(chrom))
        else: # clear out xticks but plot objects can override this later
            this_ax.set_xlabel('')
            this_ax.set_xticks([])
            
        plot_object_subset = plot_objects[ax_idx]
        
        # Set default plot limits (can be changed by client objects)
        this_ax.set_ylim((0,1)) 
        this_ax.set_xlim((ws, we))
        
        for plot_object in plot_object_subset:
            plot_object.set_globals(chrom=chrom, ws=ws, we=we, fig_width=fig_width, row_height=row_heights[ax_idx])            
            plot_object.plot(this_ax)
             
        # Redo legend code to grab colors and names from objects not from axes handles.
        if show_vector_legend and len(this_ax.get_legend_handles_labels()[1]):
            this_ax.legend(loc=vector_legend_loc)
        
#         if show_x_span and ax_idx == 0:
#             myplots.annotate_ax_span(this_ax, **x_span_kwargs)
            
    return fig    
    
    
def generate_coverage_vector(crd_df, genome=DEFAULT_GENOME):
    """
    Given a BED dataframe of CRDs, returns an int numpy.Array
    of size :param:`chrom_size` where each covered nucleotide is set to 1
    and uncovered nucleotides are set to 0.
    """
    chrom_size=genome.contig_lengths[TRAIN_CHROM]
    covered_nucs = numpy.zeros(chrom_size, dtype=int)
    for region_name in crd_df.index:
        covered_nucs[crd_df.loc[region_name, 'chromStart']:crd_df.loc[region_name, 'chromEnd']] = 1
    return covered_nucs

def pairwise_jaccard(matrix):
    """
    Return a matrix giving all pairwise jaccard indices of binary vectors
    contained as the columns of :param:`matrix`
    """
    try:
        this_index = matrix.columns
    except AttributeError:
        this_index = None
    else:
        matrix = numpy.array(matrix)
        
    def jaccard_binary_vec(vec_a, vec_b):
        return numpy.logical_and(vec_a, vec_b).sum() / numpy.logical_or(vec_a, vec_b).sum()
        
    n = matrix.shape[1]
    result = numpy.zeros((n, n))
    for col_idx_a, col_idx_b in itertools.combinations(numpy.arange(n), 2):
#         print(col_idx_a, col_idx_b)
        this_jaccard = jaccard_binary_vec(matrix[:,col_idx_a], matrix[:,col_idx_b])
        result[col_idx_a, col_idx_b] = this_jaccard
        result[col_idx_b, col_idx_a] = this_jaccard
    if len(this_index):
        return pandas.DataFrame(result, index=this_index, columns=this_index)
    else:
        return result                         
        
        
def generate_vis_vector_from_peaks(annotations,
                                     data_series,
                                     chrom,
                                     ws=0, we=0,
                                     data_transform=None, 
                                     peak_center_only=True,
                                     interpolate=False,
                                     convolution_kernel=toolbox.gaussian_kernel(66),
                                     genome=DEFAULT_GENOME):
    """
    Can be used on either feature vectors or feature deltas
    """
#     toolbox.validate_param('method', method, ['convolve', 'interpolate'])
    if we == 0:
        we = genome.contig_lengths[chrom_name]
    assert we > ws
    
    chrom_vector = numpy.zeros(shape=we-ws, dtype=float)
    
    visible_peaks = annotations.loc[(annotations.Chr == chrom) & (annotations.Start >= ws) & (annotations.End <= we)].index
    for peak_name in visible_peaks:
        this_peak_start, this_peak_end = annotations['Start'].loc[peak_name], annotations['End'].loc[peak_name]
        
        this_peak_val = data_series.loc[peak_name]
            
        if data_transform is not None:
            this_peak_val = data_transform(this_peak_val)
            
        if peak_center_only:
            this_peak_center = int((this_peak_start + this_peak_end) / 2)
#                 print('peak center: {}'.format(this_peak_center))
            chrom_vector[this_peak_center - ws] = this_peak_val   
        else:
            chrom_vector[this_peak_start - ws:this_peak_end - ws] = this_peak_val
      
    if convolution_kernel is not None:
        chrom_vector = scipy.signal.convolve(chrom_vector, convolution_kernel, mode='same')
 
    return pandas.Series(chrom_vector, index=numpy.arange(ws,we))


from scipy import ndimage
def plot_diagonal_corr_matrix(corr_matrix, figsize=(10,3), cmap='RdBu_r', max_distance=None, vmin=-1, vmax=1):
    assert corr_matrix.shape[0] == corr_matrix.shape[1]
    n = corr_matrix.shape[0]
    if not max_distance:
        max_distance = n
    rotated_n = int(n * numpy.sqrt(2))
    plot_data = numpy.array(corr_matrix)
    vert_crop_start = int(rotated_n/2 - (max_distance * numpy.sqrt(2))/2)
    vert_crop_end = rotated_n//2
    plot_data = ndimage.rotate(plot_data, 45)[vert_crop_start:vert_crop_end]
    fig, ax = plt.subplots(1, figsize=figsize)
    ax.imshow(plot_data, cmap='RdBu_r', aspect='auto', vmin=vmin, vmax=vmax)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.axis('off')
    return fig


def generate_crd_peak_subsets(crd_df, annotated_peaks):
    """
    Given a pandas.DataFrame of CRDs and an annotated peak DataFrame, return a tuple of subsets of component
    peak IDs as follows:
    
    peaks_in_crds, peaks_outside_crds, edge_peaks, interior_peaks
    """
    in_crd_peaks = set([])
    interior_peaks = set([])
    edge_peaks = set([])

    for crd_name in crd_df.index:
        crd_peaks = crd_name.split('_')
        in_crd_peaks.update(crd_peaks)
        edge_peaks.update([crd_peaks[0], crd_peaks[-1]])
        interior_peaks.update(crd_peaks[1:-1])
    outside_crd_peaks = set(annotated_peaks.index).difference(in_crd_peaks)

    return in_crd_peaks, outside_crd_peaks, edge_peaks, interior_peaks   


def generate_crd_peak_mapping(crd_df):
    """
    Given a pandas.DataFrame of CRDs, return a dictionary of lists of component peak_ids, keyed by CRD name.
    """
    peaks_by_crd = {crd_name: crd_name.split('_') for crd_name in crd_df.index} 
    crds_by_peak = {}
    for crd_name, peak_list in peaks_by_crd.items():
        for peak_id in peak_list:
            crds_by_peak[peak_id] = crd_name
            
    return peaks_by_crd, crds_by_peak
    
    