import os

import numpy
import pandas
from scipy.signal import convolve

from . import utilities
from .utilities import log_print

DEFAULT_TAG_COUNT_NORMALIZATION_TARGET = 10000000


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
        try:
            query_start = key.start
            query_end = key.stop
        except TypeError:  # Handle scalar indices
            query_start = key
            query_end = key + 1

        return self.parent_data_source.query(query_chrom=self.chrom, query_start=query_start, query_end=query_end)


class _DataVector:
    def __init__(self, chrom, parent_data_source):
        self.loc = _ChromWrapper(chrom=chrom, parent_data_source=parent_data_source)


class _DataSource:
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
        result: pandas.Series = self._query(query_chrom=query_chrom, query_start=query_start, query_end=query_end)
        if self.convolution_kernel is not None:
            result = pandas.Series(convolve(result, self.convolution_kernel, mode='same'), index=result.index)
        if self.transform:
            result = self.transform(result)
        return result

    def __getitem__(self, key):
        return _DataVector(chrom=key, parent_data_source=self)


class TagDirectory(_DataSource):
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

            self.normalization_factor = num_tags / normalize_to

    def _query(self, query_chrom, query_start, query_end, read_handling='starts'):
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
                            query_result.loc[read_start] = depth
                        else:
                            done = True

                    elif read_handling == 'reads':
                        # ToDo: Hard to do this in a streaming fashion because we don't know how far upstream to seek to capture left-overhanging reads.
                        read_len = int(line_fields[4])
                        if query_start < read_start <= query_end or query_start < read_start + read_len <= query_end:
                            print(max(read_start, query_start), min(read_start + read_len,
                                                                    query_end))
                            query_result.loc[max(read_start, query_start):min(read_start + read_len,
                                                                              query_end)] = depth  # trim to visible vector
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
