import sys
import os
import subprocess

from argparse import ArgumentParser
from .stage_one import stage_one
from .stage_two import stage_two

workingdir = os.path.abspath(os.path.dirname(__file__))
def parse_options(argv):

    parser = ArgumentParser(prog='args_oap')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{stage_one, stage_two}')

    ## stage_one parameters
    parser_stage_one = subparsers.add_parser('stage_one', help='run stage_one')
    mandatory_one = parser_stage_one.add_argument_group('MANDATORY')
    mandatory_one.add_argument('-i', help='input files directory', metavar='DIR', required=True)
    mandatory_one.add_argument('-o', help='output files directory', metavar='DIR', required=True)
    mandatory_one.add_argument('-m', help='meta data file', metavar='FILE', required=True)

    optional_one = parser_stage_one.add_argument_group('OPTIONAL')
    optional_one.add_argument('-n', help='number of threads, default 1', default=1)
    optional_one.add_argument('-f', help='the format of processed files, default fq', default='fq')
    optional_one.add_argument('-q', help='quality control of fastq sequences defualt not take effect, set to 1, then will do QC with fastp', default=0)
    optional_one.add_argument('-z', help='whether the fq files were .gz format, if 1, then firstly gzip -d, default 0', default=0)
    optional_one.add_argument('-x', help='evalue for searching 16S in usearch default 1e-10', default=1e-10)
    optional_one.add_argument('-y', help='evalue for searching universal single copy marker gene default 3', default=3)
    optional_one.add_argument('-v', help='the identity cutoff for diamond to search the USCMGs default 40', default=40)
    optional_one.add_argument('-w', help='the query cover cutoff for diamond to search the USCMGs default 65', default=65)
    optional_one.add_argument('-db', help='database, default sarg', default='')
    optional_one.add_argument('-skip', help='skip cell number and 16s calculation', action='store_true')
    parser_stage_one.set_defaults(func=stage_one)

    ## stage_one parameters
    parser_stage_two = subparsers.add_parser('stage_two', help='run stage_two')
    mandatory_two = parser_stage_two.add_argument_group('MANDATORY')
    mandatory_two.add_argument('-i', help='the potential arg reads from stage one', required=True)
    mandatory_two.add_argument('-o', help='Output prefix', required=True)
    mandatory_two.add_argument('-m', help='meta data online from stage one', required=True)

    optional_two = parser_stage_two.add_argument_group('OPTIONAL')
    optional_two.add_argument('-n', help='number of threads used for blastx, default 1', default=1)
    optional_two.add_argument('-b', help='if set then process the blastx results directly [default off], useful if user want to accelerate the stage two by running blastx paralell', default=False)
    optional_two.add_argument('-l', help='length filtering default 75 percent', default=75)
    optional_two.add_argument('-e', help='evalue cutoff, default 1e-7', default=1e-7)
    optional_two.add_argument('-d', help='identity cutoff, default 80', default=80)

    optional_two.add_argument('-db', help='reference ARG database which is blastxable after makeblastdb, default is SARG database', default=workingdir+"/DB/SARG.fasta")
    optional_two.add_argument('-struc1', help='reference ARG database structure, default is SARG database', default=workingdir+"/DB/single-component_structure.txt")
    optional_two.add_argument('-struc2', help='reference ARG database structure, default is SARG database', default=workingdir+"/DB/multi-component_structure.txt")
    optional_two.add_argument('-struc3', help='reference ARG database structure, default is SARG database', default=workingdir+"/DB/two-component_structure.txt")

    parser_stage_two.set_defaults(func=stage_two)

    if len(argv) < 2:
        parser.print_help()
        sys.exit(0)

    if len(argv) < 3:
        if argv[1] == 'stage_one':
            parser_stage_one.print_help()
        elif argv[1] == 'stage_two':
            parser_stage_two.print_help()
        else:
            parser.print_help()
        sys.exit(0)

    return(parser.parse_args(argv[1:]))

def main(argv=sys.argv):

    options = parse_options(argv)
    options.func(options)

if __name__ == "__main__":
    main(sys.argv)
