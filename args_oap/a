import sys
from argparse import ArgumentParser
import subprocess

from .stage_one import stage_one
from .stage_two import stage_two

def parse_options(argv):

    parser = ArgumentParser(prog='args_oap')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{stage_one, stage_two}')

    ## stage_one parameters
    parser_stage_one = subparsers.add_parser('stage_one', help='run stage_one')
    mandatory_one = parser_stage_one.add_argument_group('MANDATORY')
    mandatory_one.add_argument('-i','--i', help='input files directory', metavar='DIR', required=True)
    mandatory_one.add_argument('-o','--o', help='output files directory', metavar='DIR', required=True)
    mandatory_one.add_argument('-m','--m', help='meta data file', metavar='FILE', required=True)

    optional_one = parser_stage_one.add_argument_group('OPTIONAL')
    optional_one.add_argument('-n', '--n', help='number of threads, default 1', default=1)
    optional_one.add_argument('-f', '--f', help='the format of processed files, default fq', default='fq')
    optional_one.add_argument('-q', '--q', help='quality control of fastq sequences defualt not take effect, set to 1, then will do QC with fastp', default=0)
    optional_one.add_argument('-z', '--z', help='whether the fq files were .gz format, if -z, then firstly gzip -d, default(none) ', default=False, action='store_true')
    optional_one.add_argument('-x', '--x', help='evalue for searching 16S in usearch default 1e-10', default=1e-10)
    optional_one.add_argument('-y', '--y', help='evalue for searching universal single copy marker gene default 3', default=3)
    optional_one.add_argument('-v', '--v', help='the identity value for diamond to search the USCMGs default  0.45', default=0.45)

    parser_stage_one.set_defaults(func=stage_one)

    ## stage_one parameters
    parser_stage_two = subparsers.add_parser('stage_two', help='run stage_two')
    mandatory_two = parser_stage_two.add_argument_group('MANDATORY')
    mandatory_two.add_argument('-i','--i', help='the potential arg reads from stage one', required=True)
    mandatory_two.add_argument('-o','--o', help='Output prefix', required=True)
    mandatory_two.add_argument('-m','--m', help='meta data online from stage one', required=True)

    optional_two = parser_stage_two.add_argument_group('OPTIONAL')
    optional_two.add_argument('-n', '--n', help='lnumber of threads used for blastx, default 1', default=1)
    optional_two.add_argument('-b', '--b', help='if set then process the blastx results directly [default off], useful if user want to accelerate the stage two by running blastx paralell', default=False, action='store_true')
    optional_two.add_argument('-l', '--l', help='length cutoff, default 25', default=25)
    optional_two.add_argument('-e', '--e', help='evalue cutoff, default 1e-7', default=1e-7)
    optional_two.add_argument('-d', '--d', help='identity cutoff, default 80', default=80)

    parser_stage_two.set_defaults(func=stage_two)

    if len(argv) < 2:
        parser.print_help()
        sys.exit()

    if len(argv) < 3:
        if argv[1] == 'stage_one':
            parser_stage_one.print_help()
        elif argv[1] == 'stage_two':
            parser_stage_two.print_help()
        else:
            parser.print_help()
        sys.exit()

    return(parser.parse_args(argv[1:]))

def main(argv=sys.argv):

    options = parse_options(argv)
    options.func(options)

if __name__ == "__main__":
    main(sys.argv)
