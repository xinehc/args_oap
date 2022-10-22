import os
import subprocess
import sys

from .utils import *
from .settings import logger, settings

def make_db(file):
    '''
    Make database. The type of the database will be determined based on diamond's return code.
    '''
    if not os.path.isfile(file):
        logger.critical('File <{}> cannot be found. Please check input file (-i/--infile).'.format(file))
        sys.exit(2)

    logger.info('Building database of <{}> ...'.format(file))
    
    ## diamond
    subp = subprocess.run(['diamond', 'makedb', '--in', file, '--db', file, '--quiet'], check=False, stderr=subprocess.PIPE)
    
    if subp.returncode == 1: # auto choose dbtype by checking diamond's return code
        if subp.stderr == b'Error: The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.\n':
            subprocess.run(['bwa', 'index', file], check=True, stderr=subprocess.DEVNULL)
            dbtype = 'nucl'
        else:
            logger.critical('Cannot build diamond database of <{}> ({}). Please check input file (-i/--infile).'.format(file, subp.stderr))
            sys.exit(2)
    else:
        dbtype = 'prot'

    subprocess.run(['makeblastdb', '-in', file, '-dbtype', dbtype, '-out', file], check=True, stdout=subprocess.DEVNULL)
    logger.info('Finished.')

def run_make_db(options):
    make_db(vars(options).get('infile'))