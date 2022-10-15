import sys

from argparse import ArgumentParser
from .stage_one import run_stage_one
from .stage_two import run_stage_two
from .utils import *
from . import __version__

def parse_options(argv):

    parser = ArgumentParser(description='ARGs-OAP v{}:  online analysis pipeline for antibiotic resistance genes detection.'.format(__version__))
    subparsers = parser.add_subparsers(help='descriptions', metavar='{stage_one, stage_two, make_db}')

    ## stage_one parameters
    parser_stage_one = subparsers.add_parser('stage_one', help='run stage_one pipeline')

    required = parser_stage_one.add_argument_group('required arguments')
    optional = parser_stage_one.add_argument_group('optional arguments')
    database = parser_stage_one.add_argument_group('database arguments')

    required.add_argument('-i', '--indir', required=True, metavar='DIR', help='Input folder.')
    required.add_argument('-o', '--outdir', required=True, metavar='DIR', help='Output folder.')

    optional.add_argument('-t', '--thread', metavar='INT', help='Threads. [8]', default=8, type=int)
    optional.add_argument('-f', '--format', metavar='FORMAT', help='Format (extension) of files in input folder (--indir). [fq]', default='fq')
    optional.add_argument('--e1', metavar='FLOAT', help='E-value cutoff for Greengenes. [1e-10]', default=1e-10, type=float)
    optional.add_argument('--e2', metavar='FLOAT', help='E-value cutoff for Universal Essential Single Copy Marker Genes. [3]', default=3, type=float)
    optional.add_argument('--id', metavar='FLOAT', help='Identity cutoff (in percentage) for Universal Essential Single Copy Marker Genes. [40]', default=40, type=float)
    optional.add_argument('--qcov', metavar='FLOAT', help='Query cover cutoff (in percentage) for Universal Essential Single Copy Marker Genes. [65]', default=65, type=float)
    optional.add_argument('--skip', help='Skip cell number and 16s calculation.', action='store_true')
    optional.add_argument('--keep', help='Keep all temporary files (*.tmp) in the output folder (--outdir).', action='store_true')

    database.add_argument('--database', metavar='FILE', help='Customized database (indexed) other than SARG. [None]', default=None)
    parser_stage_one.set_defaults(func=run_stage_one)

    ## stage_two parameters
    parser_stage_two = subparsers.add_parser('stage_two', help='run stage_two pipeline')
    required = parser_stage_two.add_argument_group('required arguments')
    optional = parser_stage_two.add_argument_group('optional arguments')
    database = parser_stage_two.add_argument_group('database arguments')

    required.add_argument('-i', '--indir', required=True, metavar='DIR', help='Input folder. Should be the output folder of stageone (--outdir), and contain files "meta_data.txt" and "extracted.fa".')

    optional.add_argument('-t', '--thread', metavar='INT', help='Threads. [8]', default=8, type=int)
    optional.add_argument('--e', metavar='FLOAT', help='E-value cutoff for target sequences. [1e-7]', default=1e-7, type=float)
    optional.add_argument('--id', metavar='FLOAT', help='Identity cutoff (in percentage) for target sequences. [80]', default=80, type=float)
    optional.add_argument('--qcov', metavar='FLOAT', help='Query cover cutoff (in percentage) for target sequences. [75]', default=75, type=float)

    database.add_argument('--database', metavar='FILE', help='Customized database (indexed) other than SARG. [None]', default=None)
    database.add_argument('--structure1', metavar='FILE', help='Customized structure file (weight 1, single-component). [None]', default=None)
    database.add_argument('--structure2', metavar='FILE', help='Customized structure file (weight 1/2, two-component). [None]', default=None)
    database.add_argument('--structure3', metavar='FILE', help='Customized structure file (weight 1/3, multi-component). [None]', default=None)
    
    parser_stage_two.set_defaults(func=run_stage_two)

    ## makedb parameters
    parser_make_db = subparsers.add_parser('make_db', help='build customized database')
    required = parser_make_db.add_argument_group('required arguments')

    required.add_argument('-i', '--infile', required=True, metavar='FILE', help='Database FASTA file. Can be either "nucleotide" or "protein".')

    parser_make_db.set_defaults(func=run_make_db)

    ## print help
    if len(argv) < 2:
        parser.print_help()
        sys.exit(0)
    elif argv[1] not in {'stage_one', 'stage_two', 'make_db', '-h'}:
        parser.print_help()
        sys.exit(0)

    if len(argv) < 3:
        if argv[1] == 'stage_one':
            parser_stage_one.print_help()
        elif argv[1] == 'stage_two':
            parser_stage_two.print_help()
        elif argv[1] == 'make_db':
            parser_make_db.print_help()
        else:
            parser.print_help()
        sys.exit(0)
    return(parser.parse_args(argv[1:]))

def main(argv=sys.argv):

    options = parse_options(argv)
    options.func(options)

if __name__ == '__main__':
    main(sys.argv)
