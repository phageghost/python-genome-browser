import collections

import intervaltree


class IntervalDict(object):
    """
    This class stores genomic regions and their locations in a way that facilitates indexing by region_id or location.
    It supports setting and deleting, and querying by position with the .overlapping() method.
    """

    def __init__(self, region_dict=None):
        """
        Create a new IntervalDict, optionally populating with the regions in <region_dict>,
        a dictionary, keyed by region_identifier, of sub-dictionaries defining the 
        regions. These subdictionaries contain, at a minimum, fields called
        `contig`, `start`, and `end` to define their genomic location.
        """
        self._regions = collections.OrderedDict()
        self._locations = {}
        if region_dict is None: region_dict = {}

        # Update with a given region_dict
        for region_id, region in list(region_dict.items()):
            self.__setitem__(region_id, region)

    def __getitem__(self, region_id):
        """
        Getter.
        :param region_id: str
        :return: region dictionary {chrom: X, start: Y, stop: Z} query result
        """
        return self._regions[region_id]

    def __setitem__(self, new_region_id, new_region):
        """
        Setter.
        :param new_region_id: str
        :param new_region: dict
        :return: None
        """
        assert 'contig' in new_region
        assert 'start' in new_region
        assert 'end' in new_region

        if new_region['end'] <= new_region['start']:
            raise InvalidCoordinates(start=new_region['start'], end=new_region['end'],
                                     message='Start coordinate {} is greater than end coordinate {} when trying to create region {}'.format(
                                         new_region['start'], new_region['end'], new_region_id))

        # If rewriting this region with different coordinates, need to first delete it from the region tree
        if new_region_id in self._regions and (
                new_region['start'] != self._regions[new_region_id]['start'] or new_region['end'] !=
                self._regions[new_region_id]['end']):
            self.__delitem__(new_region_id)
        self._regions[new_region_id] = new_region

        # Now add the updated region to the location tree
        if new_region['contig'] not in self._locations:
            # Add as a new region entry if region ID is a new ID
            self._locations[new_region['contig']] = intervaltree.IntervalTree()
        self._locations[new_region['contig']].addi(new_region['start'], new_region['end'], new_region_id)

    def __delitem__(self, key):
        """
        Delete :param:`key` from the IntervalDict.
        :param region_id: str
        :return: None
        """
        if key in self._regions:
            self._locations[self._regions[key]['contig']].removei(self._regions[key]['start'],
                                                                  self._regions[key]['end'],
                                                                  key)
            del (self._regions[key])
        else:
            raise KeyError('Region {} not found.'.format(key))

    def __contains__(self, item):
        return item in self._regions

    def __len__(self):
        """
        Get length.
        :return: int
        """
        return len(self._regions)

    def __repr__(self):
        repr_text = ''
        for region_id, region in list(self._regions.items()):
            repr_text += '{}: {}\n'.format(region_id, ', '.join(
                ['='.join([str(pair) for pair in item]) for item in list(region.items())]))
        return repr_text

    def __iter__(self):
        return iter(list(self._regions.keys()))

    def keys(self):
        return list(self._regions.keys())

    def values(self):
        return list(self._regions.values())

    def items(self):
        return list(self._regions.items())

    def iterkeys(self):
        return iter(list(self._regions.keys()))

    def itervalues(self):
        return iter(list(self._regions.values()))

    def iteritems(self):
        return iter(list(self._regions.items()))

    def as_tupledict(self, usenames='gene_name'):
        """
        Returns a dictionary, keyed by chromosome, of lists of region tuples 
        in the form: (name, start, end), suitable for computing overlaps 
        using the two pointer algorithm.
        
        :param:`usenames` dictates which field will be used for the name element, 
        and can be either 'gene_name' or 'ensembl_id'
        """
        assert usenames in ('gene_name', 'ensembl_id'), 'invalid value for parameter \' usenames\': {}'.format(usenames)
        results = {}
        for region in self._regions.values():
            if region['contig'] not in results:
                results[region['contig']] = []
            results[region['contig']].append((region[usenames], region['start'], region['end']))
        return results

    def overlapping(self, contig, start=0, end=0, strict=False):
        """
        If <strict> is False (default), return a new RegionDict of all regions overlapping the query region by at least one base pair
        If <strict> is True, return a new RegionDict of all regions completely enclosed by the query region
        """
        if end < start:
            raise InvalidCoordinates(start=start, end=end)

        results = IntervalDict()
        if end > 0:
            for overlap in self._locations[contig].search(start, end, strict=strict):
                results[overlap.data] = self._regions[overlap.data]
        else:
            for overlap in self._locations[contig]:
                results[overlap.data] = self._regions[overlap.data]
        return results

    def intersection(self, other):
        """
        Return the regions in self that overlap with <other>.
        """
        result = IntervalDict()
        for chrom in self._locations:
            for location in self._locations[chrom]:
                if chrom in other._locations:
                    if other._locations[chrom].search(location):
                        result[location.data] = self._regions[location.data]
        return result

    def __and__(self, other):
        return self.intersection(other)

    def difference(self, other):
        """
        Return the regions in self that do not overlap with <other>
        """
        result = IntervalDict()
        for chrom in self._locations:
            for location in self._locations[chrom]:
                if not other._locations[chrom].search(location):
                    result[location.data] = self._regions[location.data]
        return result

    def __sub__(self, other):
        return self.difference(other)

    def copy(self):
        """
        Returns a deep copy of the object
        """
        return self.deep_copy()

    def shallow_copy(self):
        """
        """
        new_copy = IntervalDict()
        new_copy._regions = self._regions
        new_copy._locations = self._locations
        return new_copy

    def deep_copy(self):
        """
        """
        new_copy = IntervalDict()

        new_copy._regions.update(self._regions)

        new_copy._locations = {chrom: intervaltree.IntervalTree(self._locations[chrom]) for chrom in self._locations}
        return new_copy

    def trim_contigs(keep_sex=True):
        """
        Modifies the IntervalDict in-place to remove any non-canonical contigs
        """
        contigs_to_remove = set(['_' in contig_name or (
                not keep_sex and sum([contig_name[-1] == char for char in ('Z', 'X', 'Y', 'W')]) > 1) for contig in
                                 self.contig_names])

        for region_name, region in self._regions.items():
            if region['contig'] in contigs_to_remove:
                del (self._regions[region_name])

        for chrom in self._locations:
            if chrom in contigs_to_remove:
                del (self._locations[chrom])

    @property
    def contig_names(self):
        return list(self._locations.keys())

    def to_bed(self, bed_filename, ucsc_header=''):
        """
        Outputs the regions to a BED6 file at :param:`bed_filename`
        
        :param:`ucsc_header` defines an optional string to include at the top of the file , typically used for
        defining a UCSC genome browser track (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        
        """
        with open(bed_filename, 'wt') as bed_file:
            if ucsc_header:
                bed_file.write(ucsc_header.strip() + '\n')
            for contig in toolbox.numerical_string_sort(list(self._locations.keys())):
                for interval in sorted(list(self._locations[contig].items()), key=lambda x: x.begin):
                    region_name = interval.data
                    region_info = self._regions[region_name]

                    if 'strand' in region_info:
                        this_strand = region_info['strand']
                    else:
                        this_strand = 0

                    bed_file.write(
                        '{chrom}\t{chromStart}\t{chromEnd}\t{uniqueID}\t{dummy}\t{Strand}\n'.format(chrom=contig,
                                                                                                    chromStart=interval.begin,
                                                                                                    chromEnd=interval.end,
                                                                                                    uniqueID=region_name,
                                                                                                    dummy=' ',
                                                                                                    Strand=this_strand))

    @classmethod
    def from_bed(cls, bed_filename):
        """
        Generates a regions object from the intervals in :param:`bed_filename`. The first column of the BED is treated as an index and
        the next 4 fields of are assumed to be contig, start, end, strand. Any remaining fields will be added as data fields.
        """

        peak_regions = cls()

        bed_data = pandas.read_csv(bed_filename, index_col=0, comment='#', sep='\t')

        if len(bed_data.columns) > 4:
            additional_fields = bed_data.columns[4:]

        for region_num, region_name in enumerate(bed_data.index):
            contig = bed_data.iloc[region_num, 0]
            start = bed_data.iloc[region_num, 1]
            end = bed_data.iloc[region_num, 2]
            strand = bed_data.iloc[region_num, 3]

            if strand == '0' or strand == '+':
                strand = 1
            else:
                strand = -1

            this_region = {'contig': contig, 'start': start, 'end': end, 'strand': strand}

            for field_number, field_name in zip(range(4, len(bed_data.columns)), additional_fields):
                this_region[field_name] = bed_data.iloc[region_num, field_number]

            peak_regions[region_name] = this_region

        return peak_regions

    @classmethod
    def from_homer_motif_bed(cls, bed_filename):
        """
        Generates a regions object from the intervals in :param:`bed_filename` which should point to a file
        containing the output of a HOMER motif scan.
        """

        peak_regions = cls()

        bed_data = pandas.read_csv(bed_filename, index_col=None, header=None, comment='#', sep='\t')

        for region_num, region_name in enumerate(bed_data.index):
            if region_num % 100000 == 0:
                print('Loading region {} of {} ({:>0.2} %)'.format(region_num + 1, bed_data.shape[0],
                                                                   (region_num + 1) / bed_data.shape[0] * 100))
            contig = bed_data.iloc[region_num, 0]
            start = bed_data.iloc[region_num, 1]
            end = bed_data.iloc[region_num, 2]
            motif_name = bed_data.iloc[region_num, 3]
            motif_score = bed_data.iloc[region_num, 4]
            strand = bed_data.iloc[region_num, 5]

            if strand == '0' or strand == '+':
                strand = 1
            else:
                strand = -1

            this_region = {'contig': contig, 'start': start, 'end': end, 'strand': strand, 'motif_score': motif_score,
                           'motif_name': motif_name}

            peak_regions[region_name] = this_region

        return peak_regions

    @classmethod
    def from_homer_peaks(cls, peak_filename):
        """
        Generates a regions object from the intervals in :param:`peak_filename`
        """

        peak_regions = cls()

        with open(peak_filename, 'rt') as peak_file:
            for line in peak_file:
                if not line.startswith('#'):
                    split_line = line.split('\t')
                    region_name = split_line[0]
                    contig = split_line[1]
                    start = int(split_line[2])
                    end = int(split_line[3])
                    if split_line[4] == '0' or split_line[4] == '+':
                        strand = 1
                    else:
                        strand = -1
                    peak_regions[region_name] = {'contig': contig, 'start': start, 'end': end, 'strand': strand}

        return peak_regions
