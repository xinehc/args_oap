
from argparse import ArgumentParser
from .args_oap import run_stage_one
from .stage_two import run_stage_two
from .make_db import run_make_db
from . import __version__


def parse_options():

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
    optional.add_argument('-f', '--format', metavar='EXT', help='Format (extension) of files in input folder (--indir), gzipped (*.gz) optional. [fq]', default='fq')
    optional.add_argument('--e1', metavar='FLOAT', help='E-value cutoff for GreenGenes. [1e-10]', default=1e-10, type=float)
    optional.add_argument('--e2', metavar='FLOAT', help='E-value cutoff for Essential Single Copy Marker Genes. [3]', default=3, type=float)
    optional.add_argument('--id', metavar='FLOAT', help='Identity cutoff (in percentage) for Essential Single Copy Marker Genes. [45]', default=45, type=float)
    optional.add_argument('--qcov', metavar='FLOAT', help='Query cover cutoff (in percentage) for Essential Single Copy Marker Genes. [0]', default=0, type=float)
    optional.add_argument('--skip', help='Skip cell number and 16S calculation.', action='store_true')
    optional.add_argument('--keep', help='Keep all temporary files (*.tmp) in the output folder (--outdir).', action='store_true')

    database.add_argument('--database', metavar='FILE', help='Customized database (indexed) other than SARG. [None]', default=None)
    parser_stage_one.set_defaults(func=run_stage_one)

    ## stage_two parameters
    parser_stage_two = subparsers.add_parser('stage_two', help='run stage_two pipeline')
    required = parser_stage_two.add_argument_group('required arguments')
    optional = parser_stage_two.add_argument_group('optional arguments')
    database = parser_stage_two.add_argument_group('database arguments')

    required.add_argument('-i', '--indir', required=True, metavar='DIR', help='Input folder. Should be the output folder of stage_one (contain "meta_data.txt" and "extracted.fa").')

    optional.add_argument('-t', '--thread', metavar='INT', help='Threads. [8]', default=8, type=int)
    optional.add_argument('-o', '--outdir', metavar='DIR', help='Output folder, if not given then samve as input folder (--indir). [None]', default=None)
    optional.add_argument('--e', metavar='FLOAT', help='E-value cutoff for target sequences. [1e-7]', default=1e-7, type=float)
    optional.add_argument('--id', metavar='FLOAT', help='Identity cutoff (in percentage) for target sequences. [80]', default=80, type=float)
    optional.add_argument('--qcov', metavar='FLOAT', help='Query cover cutoff (in percentage) for target sequences. [75]', default=75, type=float)
    optional.add_argument('--length', metavar='INT', help='Aligned length cutoff (in amino acid) for target sequences. [25]', default=25, type=int)
    optional.add_argument('--blastout', metavar='FILE', help='BLAST output (-outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"), if given then skip BLAST. [None]', default=None)

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

    return parser

def call_stage_one(args,logger):
    parser = parse_options()
    parser = parser.parse_args(args)
    parser.logger = logger
    print(vars(parser))
    run_stage_one(parser)


def call_stage_two(args,logger):
    parser = parse_options()
    parser = parser.parse_args(args)
    parser.logger = logger
    print(vars(parser))
    run_stage_two(parser)
