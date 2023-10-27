import os
import subprocess
import sys

from .settings import logger


def make_db(file):
    '''
    Make database, also determine the type.
    '''
    if not os.path.isfile(file):
        logger.critical(f'File <{file}> cannot be found. Please check input file (-i/--infile).')
        sys.exit(2)

    logger.info(f'Building database of <{file}> ...')

    ## try making db using diamond, if fail use bwa
    try:
        subprocess.run([
            'diamond', 'makedb',
            '--in', file,
            '--db', file,
            '--quiet'], check=True, stderr=subprocess.DEVNULL)
        dbtype = 'prot'

    except subprocess.CalledProcessError:
        try:
            subprocess.run([
                'bwa', 'index',
                file], check=True, stderr=subprocess.DEVNULL)
            dbtype = 'nucl'

        except subprocess.CalledProcessError:
            pass

    ## if makeblastdb fail then stop
    try:
        subprocess.run([
            'makeblastdb',
            '-in', file,
            '-dbtype', dbtype,
            '-out', file], check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    except subprocess.CalledProcessError:
        logger.critical(f'Cannot build database of <{file}>. Please check the format of input file (-i/--infile).')
        sys.exit(2)

    logger.info('Finished.')


def run_make_db(options):
    make_db(vars(options).get('infile'))
