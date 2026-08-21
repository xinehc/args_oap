import gzip
import sys
import re
from .settings import logger


def buffer_count(file):
    '''
    Count number of lines of a fa/fq file, then divide by 2/4 to get reads.
    '''
    nlines = 0
    if re.search('.gz$', file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file)

    char = f.read(1)  # get first character > or @
    f.seek(0)  # roll back

    buf_size = 1024 * 1024
    read_f = f.read

    buf = read_f(buf_size)
    while buf:
        nlines += buf.count('\n')
        buf = read_f(buf_size)

    f.close()
    if char == '>':
        return int(nlines / 2)
    elif char == '@':
        return int(nlines / 4)
    else:
        logger.critical(f'Unrecognized file format <{file}>. Only fasta or fastq supported.')
        sys.exit(2)


def nucleotide_count(file):
    '''
    Count A/T/C/G nucleotides in a fasta/fastq file.
    '''
    if re.search('.gz$', file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file)

    char = f.read(1)
    f.seek(0)

    if char == '>':
        sequences = (line.rstrip().upper() for line in f if not line.startswith('>'))
    elif char == '@':
        sequences = (line.rstrip().upper() for i, line in enumerate(f) if i % 4 == 1)
    else:
        f.close()
        logger.critical(f'Unrecognized file format <{file}>. Only fasta or fastq supported.')
        sys.exit(2)

    counts = {base: 0 for base in 'ATCG'}
    for sequence in sequences:
        for base in counts:
            counts[base] += sequence.count(base)

    f.close()
    return counts


def simple_count(file):
    '''
    Simple count bps of a file.
    '''
    nbps, nlines = 0, 0
    f = open(file)

    char = '>'
    record = False
    for line in f:
        if line[0] == char:
            nlines += 1
            record = True  # record next row

        elif record:
            nbps += len(line) - 1
            record = False

    f.close()
    return nbps, nlines
