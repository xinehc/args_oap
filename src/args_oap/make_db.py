import os
import subprocess
import sys

from .utils import *
from .settings import settings

def make_db(file,logger):
    '''
    Make database. The type of the database will be determined based on diamond's return code.
    '''
    if not os.path.isfile(file):
        logger.critical('File <{}> cannot be found. Please check input file (-i/--infile).'.format(file))
        sys.exit(2)

    logger.info('Building database of <{}> ...'.format(file))
    
    ## diamond
    subp = subprocess.run(['diamond', 'makedb', '--in', file, '--db', file, '--quiet'], check=False, stderr=subprocess.PIPE, creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    
    if subp.returncode == 1: # auto choose dbtype by checking diamond's return code
        if subp.stderr == b'Error: The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.\n' or subp.stderr == b'Error: The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.\r\n':
            subprocess.run(['bwa', 'index', file], check=True, stderr=subprocess.DEVNULL, creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            dbtype = 'nucl'
        else:
            logger.critical('Cannot build diamond database of <{}> ({}). Please check input file (-i/--infile).'.format(file, subp.stderr))
            raise RuntimeError('Cannot build diamond database of <{}> ({}). Please check input file (-i/--infile).'.format(file, subp.stderr))
    else:
        dbtype = 'prot'

    subprocess.run(['makeblastdb', '-in', os.path.basename(file), '-dbtype', dbtype, '-out', os.path.basename(file)], cwd=os.path.dirname(file), check=True, stdout=subprocess.DEVNULL, creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    logger.info('Finished.')

def run_make_db(options):
    make_db(vars(options).get('infile'))