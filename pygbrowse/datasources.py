import os

import numpy
import pandas
import pysam
from scipy.signal import convolve

from . import utilities
from .utilities import log_print

DEFAULT_TAG_COUNT_NORMALIZATION_TARGET = 10000000
DEFAULT_FEATURE_SOURCES = ('ensembl', 'havana', 'ensembl_havana')
DEFAULT_GENE_TYPES = (
'gene', 'RNA', 'mt_gene', 'lincRNA_gene', 'miRNA_gene', 'ncRNA_gene', 'rRNA_gene', 'snRNA_gene', 'snoRNA_gene',
'processed_transcript')
DEFAULT_TRANSCRIPT_TYPES = ('mRNA', 'transcript', 'lincRNA', 'lnc_RNA', 'miRNA', 'ncRNA', 'snRNA', 'snoRNA')
DEFAULT_COMPONENT_TYPES = ('CDS', 'three_prime_UTR', 'five_prime_UTR')


# DEFAULT_MAXIMUM_TRANSCRIPT_SUPPORT = 5

# ToDo: For each class, allow option of loading into memory or leaving on disk (where applicable)
# ToDo: Add a transform function and smoothing.
# ToDo: Add indexing of on-disk csv-like files


class _ChromWrapper:
    def __init__(self, chrom, parent_data_source):
        self.chrom = chrom
        self.parent_data_source = parent_data_source

    def __getitem__(self, key):
        # print(key)
        # ToDo: Add support for step argument
        try:  # See if key is a slice
            query_start = key.start
            query_end = key.stop
        except TypeError:  # if not, treat as a scalar index
            query_start = key
            query_end = key + 1

        return self.parent_data_source.query(query_chrom=self.chrom, query_start=query_start, query_end=query_end)


class _DataVector:
    def __init__(self, chrom, parent_data_source):
        self.loc = _ChromWrapper(chrom=chrom, parent_data_source=parent_data_source)


class _VectorDataSource:
    # ToDo: Add methods for arithmetic and such, as done for old Pileups class
    def __init__(self, transform=None, smoothing_bandwidth=0):
        self.transform = transform
        if smoothing_bandwidth:
            self.convolution_kernel = utilities.gaussian_kernel(smoothing_bandwidth)
        else:
            self.convolution_kernel = None

    def _query(self, query_chrom, query_start, query_end):
        print('Stub method -- must be overridden by inheritors')

    def query(self, query_chrom, query_start, query_end):
        result = self._query(query_chrom=query_chrom, query_start=query_start, query_end=query_end)
        if self.convolution_kernel is not None:
            result = pandas.Series(convolve(result, self.convolution_kernel, mode='same'), index=result.index)
        if self.transform:
            result = self.transform(result)
        return result

    def __getitem__(self, key):
        return _DataVector(chrom=key, parent_data_source=self)


class TagDirectory(_VectorDataSource):
    tag_strand_translator = {0: '+', 1: '-'}

    def __init__(self, tag_directory_path, normalize_to=DEFAULT_TAG_COUNT_NORMALIZATION_TARGET, transform=None,
                 smoothing_bandwidth=0):
        super(TagDirectory, self).__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)

        self.tag_directory_path = tag_directory_path

        if normalize_to:
            # extract total tag count from tagInfo.txt
            tag_info_fname = os.path.join(tag_directory_path, 'tagInfo.txt')
            with open(tag_info_fname, 'rt') as tag_info_file:
                sizeline = tag_info_file.readlines()[1].strip().split('\t')
            num_tags = int(float(sizeline[2]))

            self.normalization_factor = normalize_to / num_tags

    def _query(self, query_chrom, query_start, query_end, read_handling='reads'):
        # ToDo: Add argument validation to all functions and methods with string parameters
        # ToDo: Add verbosity-based logging output
        # ToDo; Compare performance with memory-mapped pandas DataFrames
        query_result = pandas.Series(numpy.zeros(query_end - query_start), index=numpy.arange(query_start, query_end))

        tag_filename = os.path.join(self.tag_directory_path, '{}.tags.tsv'.format(query_chrom))
        start_offset = utilities.binary_search_tag_file(tag_filename=tag_filename, search_target=query_start + 1)

        done = False
        with open(tag_filename, 'rt') as tag_file:
            tag_file.seek(start_offset)
            # print(start_offset)
            while not done:
                line_fields = tag_file.readline().strip().split('\t')
                # print(line_fields)
                if len(line_fields) > 1:
                    # chrom = line_fields[0]
                    read_start = int(line_fields[1]) - 1
                    # strand = self.tag_strand_translator[int(line_fields[2])]
                    depth = float(line_fields[3])

                    if read_handling == 'starts':
                        assert read_start > query_start
                        if read_start < query_end:
                            query_result.loc[read_start] += depth
                        else:
                            done = True

                    elif read_handling == 'reads':
                        # ToDo: Hard to do this in a streaming fashion because we don't know how far upstream to seek to capture left-overhanging reads.
                        read_len = int(line_fields[4])
                        if query_start < read_start <= query_end or query_start < read_start + read_len <= query_end:
                            # print(max(read_start, query_start), min(read_start + read_len,
                            #                                         query_end))
                            query_result.loc[max(read_start, query_start):min(read_start + read_len,
                                                                              query_end)] += depth  # trim to visible vector
                        else:
                            done = True

        query_result *= self.normalization_factor

        return query_result


class IntervalData:
    HOMER_PEAKFILE_HEADER_ROW = 39
    HOMER_PEAKFILE_COLUMN_RENAMER = {'chr': 'chrom', 'start': 'chromStart', 'end': 'chromEnd'}
    HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER = {'Chr': 'chrom', 'Start': 'chromStart', 'End': 'chromEnd', 'Strand': 'strand'}

    def __init__(self, interval_data, format='bed'):
        """
        Loads genomic interval information in various formats and stores them in a standardized form as a
        pandas.DataFrame in self.data.

        :param:`interval_data` should be a pandas.DataFrame representing BED-formatted genomic data, or,
        alternatively, a filename pointing to one of the following file formats:

        * A BED file
        * A HOMER peak file
        * A HOMER annotated peak file.

        If a filename is passed instead of a DataFrame, :param:`format` should be specified. Allowed values are:
            'bed', 'homer', 'homer_annotated'

        :param interval_data:
        :param format:
        """
        try:
            _ = interval_data.loc[:, ['chrom', 'chromStart', 'chromEnd', 'strand']]

        except KeyError:  # maybe it's a BED DataFrame without column names?
            log_print('Guessing this is a BED-style DataFrame without column names')

            assert interval_data.shape[1] >= 3, 'Not enough columns (got {})!'.format(interval_data.shape[1])

            if interval_data.shape[1] >= 6:  # assume name is still separate column
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'name',
                                     'score', 'strand'] + list(self.data.columns)[6:]
                self.data.index = self.data['name']

            elif interval_data.shape[1] == 5:  # assume name has been made the index and deleted from the columns
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'score',
                                     'strand']
            else:
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'score',
                                     'strand'][:interval_data.shape[1] - 5]

            self.data.index.name = 'IntervalID'

        except (AttributeError,):  # guessing it's a filename string
            log_print('Guessing {} is a filename'.format(interval_data))
            # if format == 'auto':
            #     extension = filename.split('.')[-1]
            #     if extension.lower() == 'bed':
            #         format = 'bed'
            #     elif extension.lower() == 'homer':
            #         # ToDo: Add more sophisticated methods of detecting formats since, e.g. .txt can refer to many.
            #         format = 'homer'
            # ToDo: allow filename parameters to be file handles -- autodetect

            if format == 'bed':
                self.data = pandas.read_csv(interval_data, sep='\t', index_col=3, comment='#', header=None,
                                            names=['chrom', 'chromStart', 'chromEnd', 'score', 'strand'])
            elif format == 'homer':
                self.data = pandas.read_csv(interval_data, sep='\t', index_col=0, header=self.HOMER_PEAKFILE_HEADER_ROW)
                self.data.index.name = self.data.index.name[1:]
                self.data = self.data.rename(columns=self.HOMER_PEAKFILE_COLUMN_RENAMER)

            elif format == 'homer_annotated':
                self.data = pandas.read_csv(interval_data, index_col=0, sep='\t')
                self.data.index.name = self.data.index.name.split(' ')[0]
                self.data = self.data.rename(columns=self.HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER)

        else:  # seems to be a properly-formatted DataFrame so just  store it
            self.data = interval_data

        self.data = self.data.sort_values(['chrom', 'chromStart'])


class _GeneModels():
    def __init__(self):
        pass

    def _query(self, query_chromosome, query_start, query_end):
        print('Must be overridden by inheritors!')

    def query(self, chromosome, start, end):
        return self._query(query_chromosome=chromosome, query_start=start, query_end=end)


class Gff3Annotations(_GeneModels):
    def __init__(self,
                 gff3_filename,
                 incoming_chromosome_name_converter=lambda x: utilities.convert_chromosome_name(x, dialect='ensembl'),
                 outgoing_chromosome_name_converter=lambda x: utilities.convert_chromosome_name(x, dialect='ucsc'),
                 feature_sources=DEFAULT_FEATURE_SOURCES,
                 gene_types=DEFAULT_GENE_TYPES,
                 transcript_types=DEFAULT_TRANSCRIPT_TYPES,
                 component_types=DEFAULT_COMPONENT_TYPES,
                 # maximum_transcript_support=DEFAULT_MAXIMUM_TRANSCRIPT_SUPPORT
                 ):

        super(Gff3Annotations, self).__init__()

        self.tabix_file = pysam.TabixFile(gff3_filename)
        self.incoming_chromosome_name_converter = incoming_chromosome_name_converter
        self.outgoing_chromosome_name_converter = outgoing_chromosome_name_converter
        self.feature_sources = feature_sources
        self.gene_types = gene_types
        self.transcript_types = transcript_types
        self.component_types = component_types
        # self.maximum_transcript_support = maximum_transcript_support

    def _query(self, query_chromosome, query_start, query_end):
        gene_names_to_ensembl_ids = {}
        genes = {}
        transcripts = {}
        components = {}
        component_num = 0  # serial index for components without IDs

        query_rows = self.tabix_file.fetch(self.incoming_chromosome_name_converter(query_chromosome), query_start,
                                           query_end)

        for line in query_rows:
            split_line = line.strip('\n').split('\t')
            source, feature_type = split_line[1], split_line[2]

            if source in self.feature_sources:
                contig = split_line[0]
                start = int(split_line[3])
                end = int(split_line[4])
                strand = split_line[6]

                fields = dict(field_value_pair.split('=') for field_value_pair in split_line[8].split(';'))
                # print(line_num, line)
                if feature_type in self.gene_types:
                    ensembl_id = fields['ID']
                    gene_name = fields['Name']
                    #                     assert ensembl_id not in genes, 'Duplicate entry for gene {} on line {}'.format(ensembl_id,
                    #                                                                                                     line_num)

                    genes[ensembl_id] = {'contig': contig,
                                         'start': start - 1,  # convert 1-based to 0-based
                                         'end': end,
                                         'strand': strand,
                                         'transcripts': []}

                    genes[ensembl_id].update(fields)
                    if gene_name not in gene_names_to_ensembl_ids:
                        gene_names_to_ensembl_ids[gene_name] = []
                    gene_names_to_ensembl_ids[gene_name].append(ensembl_id)
                    # print('\t added gene {}'.format(ensembl_id))

                elif feature_type in self.transcript_types:
                    parent = fields['Parent']
                    # print('\ttranscript parent {}'.format(parent))
                    # try:
                    #     transcript_support_level = int(fields['transcript_support_level'].split(' ')[0])
                    # except ValueError:
                    #     passed_support_filter = False
                    # else:
                    #     passed_support_filter = transcript_support_level < self.maximum_transcript_support

                    if parent in genes:
                        ensembl_id = fields['ID']
                        transcripts[ensembl_id] = {'contig': contig,
                                                   'start': start - 1,  # convert 1-based to 0-based
                                                   'end': end,
                                                   'strand': strand,
                                                   'components': []}
                        transcripts[ensembl_id].update(fields)

                        genes[parent]['transcripts'].append(ensembl_id)
                        # print('\t added transcript {} with parent {}'.format(ensembl_id, parent))


                elif feature_type in self.component_types:
                    parent = fields['Parent']
                    # print('\thas parent {}. {}'.format(parent, parent in transcripts))
                    if parent in transcripts:
                        if 'exon_id' in fields:
                            ensembl_id = fields['exon_id']
                        else:
                            ensembl_id = str(component_num)
                            component_num += 1

                        components[ensembl_id] = {'contig': contig,
                                                  'start': start - 1,  # convert 1-based to 0-based
                                                  'end': end,
                                                  'strand': strand,
                                                  'type': feature_type}
                        components[ensembl_id].update(fields)
                        transcripts[parent]['components'].append(ensembl_id)

        return genes, transcripts, components, gene_names_to_ensembl_ids
