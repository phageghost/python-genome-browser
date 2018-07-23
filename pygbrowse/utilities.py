import datetime
import os

from . import romannumerals


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
    Find the offset (in bytes) in :param:`tagf_filename` that corresponds
    to the start of the first tag that is equal to or greater than :param:`search_target`.

    If none of the reads have a start position greater than :param:`search_target`,
    return None.

    Note that positions in tag files have a 1-based index.
    """
    filesize = os.path.getsize(tag_filename)
    search_window_start = 0
    search_window_end = filesize - 1
    guess = int((search_window_start + search_window_end) / 2)

    with open(tag_filename, 'rt') as tag_file:
        while search_window_end - search_window_start > 1:
            tag_file.seek(guess)
            if guess > 0:
                _ = tag_file.readline()  # read forward to get to a line start

            thisline = tag_file.readline().strip()

            if thisline == '':
                # We've reached the end of the file and the reads are still upstream of the target
                return None

            genomic_start = int(thisline.split('\t')[1])

            print(search_window_start, guess, search_window_end, genomic_start, )

            if genomic_start < search_target:
                search_window_start = guess

            elif genomic_start >= search_target:
                search_window_end = guess

            guess = int((search_window_start + search_window_end) / 2)

    return guess

def pretty_now():
    """
    Returns the current date/time in a nicely formatted string (without decimal seconds)
    """
    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%b-%d %H:%M:%S')


def log_print(message, tabs=1):
    """
    Print a chunk of text preceded by the timestamp and an optional number of tabs.

    :param message:
    :param tabs:
    :return:
    """
    print('{}{}{}'.format(pretty_now(), '\t' * tabs, message))
