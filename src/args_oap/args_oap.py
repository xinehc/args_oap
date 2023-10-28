import sys

from argparse import ArgumentParser
from .stage_one import run_stage_one
from .stage_two import run_stage_two
from .make_db import run_make_db
from . import __version__


## stage_one parameters
def parse_stage_one(parser):
    parser_stage_one = parser.add_parser('stage_one', help='run stage_one pipeline')

    required = parser_stage_one.add_argument_group('required arguments')
    optional = parser_stage_one.add_argument_group('optional arguments')
    database = parser_stage_one.add_argument_group('database arguments')

    required.add_argument(
        '-i',
        '--indir',
        required=True,
        metavar='DIR',
        help='Input folder.')

    required.add_argument(
        '-o',
        '--outdir',
        required=True,
        metavar='DIR',
        help='Output folder.')

    optional.add_argument(
        '-t',
        '--thread',
        metavar='INT',
        default=8,
        type=int,
        help='Number of threads. [8]')

    optional.add_argument(
        '-f',
        '--format',
        metavar='EXT',
        default='fq',
        help='File format in input folder (--indir), gzipped (*.gz) optional. [fq]')

    optional.add_argument(
        '--e1',
        metavar='FLOAT',
        default=1e-10,
        type=float,
        help='E-value cutoff for GreenGenes. [1e-10]')

    optional.add_argument(
        '--e2',
        metavar='FLOAT',
        default=3,
        type=float,
        help='E-value cutoff for Essential Single Copy Marker Genes. [3]')

    optional.add_argument(
        '--id',
        metavar='FLOAT',
        default=45,
        type=float,
        help='Identity cutoff (in percentage) for Essential Single Copy Marker Genes. [45]')

    optional.add_argument(
        '--qcov',
        metavar='FLOAT',
        default=0,
        type=float,
        help='Query cover cutoff (in percentage) for Essential Single Copy Marker Genes. [0]')

    optional.add_argument(
        '--skip',
        action='store_true',
        help='Skip cell number and 16S copy number estimation.')

    optional.add_argument(
        '--keep',
        action='store_true',
        help='Keep all temporary files (*.tmp) in output folder (--outdir).')

    database.add_argument(
        '--database',
        metavar='FILE',
        default=None,
        help='Customized database (indexed) other than SARG. [None]')

    parser_stage_one.set_defaults(func=run_stage_one)


## stage_two parameters
def parse_stage_two(parser):
    parser_stage_two = parser.add_parser('stage_two', help='run stage_two pipeline')

    required = parser_stage_two.add_argument_group('required arguments')
    optional = parser_stage_two.add_argument_group('optional arguments')
    database = parser_stage_two.add_argument_group('database arguments')

    required.add_argument(
        '-i',
        '--indir',
        required=True,
        metavar='DIR',
        help='Input folder. Should be the output folder of stage_one (containing <metadata.txt> and <extracted.fa>).')

    optional.add_argument(
        '-t',
        '--thread',
        metavar='INT',
        default=8,
        type=int,
        help='Number of threads. [8]')

    optional.add_argument(
        '-o',
        '--outdir',
        metavar='DIR',
        default=None,
        help='Output folder, if not given then same as input folder (--indir). [None]')

    optional.add_argument(
        '--e',
        metavar='FLOAT',
        default=1e-7,
        type=float,
        help='E-value cutoff for target sequences. [1e-7]')

    optional.add_argument(
        '--id',
        metavar='FLOAT',
        default=80,
        type=float,
        help='Identity cutoff (in percentage) for target sequences. [80]')

    optional.add_argument(
        '--qcov',
        metavar='FLOAT',
        default=75,
        type=float,
        help='Query cover cutoff (in percentage) for target sequences. [75]')

    optional.add_argument(
        '--length',
        metavar='INT',
        default=25,
        type=int,
        help='Aligned length cutoff (in amino acid) for target sequences. [25]')

    optional.add_argument(
        '--blastout',
        metavar='FILE',
        default=None,
        help='BLAST output (-outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"), if given then skip BLAST. [None]')

    database.add_argument(
        '--database',
        metavar='FILE',
        default=None,
        help='Customized database (indexed) other than SARG. [None]')

    database.add_argument(
        '--structure1',
        metavar='FILE',
        default=None,
        help='Customized structure file (weight 1, single-component). [None]')

    database.add_argument(
        '--structure2',
        metavar='FILE',
        default=None,
        help='Customized structure file (weight 1/2, two-component). [None]')

    database.add_argument(
        '--structure3',
        metavar='FILE',
        default=None,
        help='Customized structure file (weight 1/3, multi-component). [None]')

    parser_stage_two.set_defaults(func=run_stage_two)


## makedb parameters
def parse_make_db(parser):
    parser_make_db = parser.add_parser('make_db', help='build customized database')
    required = parser_make_db.add_argument_group('required arguments')

    required.add_argument(
        '-i',
        '--infile',
        required=True,
        metavar='FILE',
        help='Database FASTA file. Can be either nucleotide or protein.')

    parser_make_db.set_defaults(func=run_make_db)


def main(argv=sys.argv):
    '''entry point'''

    parser = ArgumentParser(description=f'ARGs-OAP v{__version__}: online analysis pipeline for antibiotic resistance genes detection.')
    subparsers = parser.add_subparsers(help='descriptions', metavar='{stage_one, stage_two, make_db}')

    ## attach parsers
    parse_stage_one(subparsers)
    parse_stage_two(subparsers)
    parse_make_db(subparsers)

    ## print help
    if len(argv) < 2 or argv[1] not in {'stage_one', 'stage_two', 'make_db', '-v', '--version'}:
        parser.print_help()
        sys.exit(0)

    if argv[1] in {'-v', '--version'}:
        print(__version__)
        sys.exit(0)

    if len(argv) < 3:
        if argv[1] == 'stage_one':
            subparsers.choices['stage_one'].print_help()
        elif argv[1] == 'stage_two':
            subparsers.choices['stage_two'].print_help()
        elif argv[1] == 'make_db':
            subparsers.choices['make_db'].print_help()
        else:
            parser.print_help()
        sys.exit(0)

    options = parser.parse_args(argv[1:])
    options.func(options)


if __name__ == '__main__':
    main(sys.argv)
