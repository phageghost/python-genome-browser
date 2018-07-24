import os

import numpy
import pandas

from . import utilities

DEFAULT_TAG_COUNT_NORMALIZATION_TARGET = 10000000


# ToDo: For each class, allow option of loading into memory or leaving on disk (where applicable)
# ToDo: Add a transform function and smoothing.

class GenomicData:
    def query(self, chrom, start, end):
        return None


class SeriesDict(GenomicData):
    def __init__(self, series_dict):
        self.series_by_chrom = series_dict

    def query(self, query_chrom, query_start, query_end):
        return self.series_by_chrom[query_chrom].loc[query_start:query_end]


class TagDirectory(GenomicData):
    tag_strand_translator = {0: '+', 1: '-'}

    def __init__(self, tag_directory_path, normalize_to=DEFAULT_TAG_COUNT_NORMALIZATION_TARGET):
        self.tag_directory_path = tag_directory_path

        if normalize_to:
            # extract total tag count from tagInfo.txt
            tag_info_fname = os.path.join(tag_directory_path, 'tagInfo.txt')
            with open(tag_info_fname, 'rt') as tag_info_file:
                sizeline = tag_info_file.readlines()[1].strip().split('\t')
            num_tags = int(float(sizeline[2]))

            self.normalization_factor = num_tags / normalize_to

    def __getitem__(self, key):
        return

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

