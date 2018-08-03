import datetime
import os
import subprocess

import numpy
from scipy.stats import norm

from . import romannumerals


# ToDo: Bring back scale bar
# ToDo: Add option for solid fill of vectors


def roundto(num, nearest):
    """
    Rounds :param:`num` to the nearest increment of :param:`nearest`
    """
    return int((num + (nearest / 2)) // nearest * nearest)


def convert_chromosome_name(chrom_string, dialect='ucsc'):
    """
    Try to auto-detect chromosome number and convert it to the specified "dialect".

    Valid dialects are "ucsc", "ensembl" and "yeast".

    :param chrom_string:
    :param source:
    :param dest:
    :return:
    """
    try:
        chrom_string = str(romannumerals.roman_to_int(chrom_string))
    except ValueError:
        pass

    if dialect == 'ensembl':
        if chrom_string == 'chrM':
            return 'dmel_mitochonrdion_genome'
        elif chrom_string[:3].lower() == 'chr':
            return chrom_string[3:]
        else:
            return chrom_string
    elif dialect == 'ucsc':
        if chrom_string == 'dmel_mitochondrion_genome':
            return 'chrM'
        elif chrom_string[:3].lower() == 'chr':
            return chrom_string
        else:
            return 'chr{}'.format(chrom_string)
    elif dialect == 'yeast':
        if chrom_string[:3].lower() == 'chr':
            chrom_string = chrom_string[3:]
        try:
            return romannumerals.int_to_roman(int(chrom_string))
        except ValueError:
            return chrom_string
    else:
        raise ValueError('Unknown dialect {}'.format(dialect))


def binary_search_tag_file(tag_filename, search_target):
    """
    Find the offset (in bytes) in :param:`tag_filename` that corresponds
    to the start of the first tag that is equal to or greater than :param:`search_target`.

    If none of the reads have a start position greater than :param:`search_target`,
    return None.

    Note that positions in tag files have a 1-based index.
    """

    def get_read_start(file_offset):
        tag_file.seek(file_offset)
        if file_offset > 0:
            _ = tag_file.readline()  # read forward to get to a line start
        this_line = tag_file.readline().strip()
        if tag_file.tell() >= filesize:
            # We've reached the end of the file and the reads are still upstream of the target
            return None
        else:
            return int(this_line.split('\t')[1])

    filesize = os.path.getsize(tag_filename)
    search_window_start = 0
    search_window_end = filesize - 1
    guess_genomic_start = -1
    guess = int((search_window_start + search_window_end) / 2)

    with open(tag_filename, 'rt') as tag_file:
        first_genomic_start = get_read_start(search_window_start)
        # last_genomic_start = get_read_position(search_window_end)

        if search_target < first_genomic_start:
            return search_window_start

        while search_window_end - search_window_start > 1:
            guess = int((search_window_start + search_window_end) / 2)
            guess_genomic_start = get_read_start(guess)

            if guess_genomic_start == None:
                return None

            # print(search_window_start, guess, search_window_end, guess_genomic_start)

            if guess_genomic_start < search_target:
                # print('\ttoo low!')
                search_window_start = guess

            elif guess_genomic_start > search_target:
                search_window_end = guess

                # print('\ttoo high!')
            else:
                # print('\tjust right!')
                break

        if guess_genomic_start == -1:
            return None

        if guess_genomic_start < search_target:
            guess += 1

        tag_file.seek(guess)
        _ = tag_file.readline()
        guess = tag_file.tell()

        return guess


def bgzip_gff(gff3_fname, bgzipped_fname):
    """
    Compress a GFF3 file in block-gzip format (requires that bgzip be accessible on the current path).

    If :param gff3_fname: ends with '.gz' assumes that the file is gzipped, otherwise assumes it is uncompressed.

    :param gzipped_fname:
    :param bgzipped_fname:
    :return:
    """
    if bgzipped_fname == gff3_fname:
        log_print('Destination and source file cannot have the same name!')

    cmd_line = '{} {} | sort -k1,1 -k4,4n | bgzip > {}'.format(('cat', 'zcat')[gff3_fname.endswith('.gz')], gff3_fname,
                                                               bgzipped_fname)
    try:
        assert os.path.isfile(gff3_fname)  # needed since no error occurs otherwise
        subprocess.check_call(cmd_line, shell=True)

    except subprocess.CalledProcessError as cpe:
        log_print('Unsuccessful. Got return code {}'.format(cpe.returncode))

    except AssertionError:
        log_print('{} not found!'.format(gff3_fname))

    else:
        log_print('Successfully generated block-gzipped file {} from {}'.format(bgzipped_fname, gff3_fname))


def generate_tabix_index(target_fname):
    """
    Index :param target_fname: with tabix. Requires that the directory in which :param:target_fname: resides is
    writeable.

    :param target_fname:
    :return:
    """
    cmd_line = 'tabix -f -p gff {}'.format(target_fname)
    try:
        return_code = subprocess.check_call(cmd_line, shell=True)
    except subprocess.CalledProcessError as cpe:
        log_print('Unsuccessful. Got return code {}'.format(cpe.returncode))
    else:
        log_print('Successfully indexed block-gzipped file {}'.format(target_fname))


def pretty_now():
    """
    Returns the current date/time in a nicely formatted string (without decimal seconds)
    """
    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%b-%d %H:%M:%S')


def log_print(message, tabs=1):
    """
    Print a chunk of text preceded by a timestamp and an optional number of tabs (default 1).

    :param message:
    :param tabs:
    :return:
    """
    print('{}{}{}'.format(pretty_now(), '\t' * tabs, message))


def gaussian_kernel(sd, sd_cutoff=3, normalize=False):
    """
    Generate and return a numpy.Array whose elements are proportional to the PDF of a normal distribution
    having standard deviation :param:`sd`.

    :param sd:
    :param sd_cutoff:
    :param normalize:
    :return:
    """
    bw = sd_cutoff * sd * 2 + 1
    midpoint = sd_cutoff * sd
    kern = numpy.zeros(bw)
    frozen_rv = norm(scale=sd)
    for i in range(bw):
        kern[i] = frozen_rv.pdf(i - midpoint)
    if normalize:
        kern = kern / kern.max()
    return kern
