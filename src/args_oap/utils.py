import subprocess
import os
import sys
import re
import gzip 

from .settings import logger

def get_filename(file, format, drop=False):
    if drop:
        return re.sub(rf'(_1|_2)?\.{format}(.gz)?$','',os.path.basename(file))
    else:
        return re.sub(rf'\.{format}(.gz)?$','',os.path.basename(file))
        
## https://stackoverflow.com/a/850962
def buffer_count(file):
    nlines = 0
    if re.search('.gz$', file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file)

    char = f.read(1) # get first character, > or @
    f.seek(0) # roll back

    buf_size = 1024 * 1024
    read_f = f.read
    
    buf = read_f(buf_size)
    while buf:
        nlines += buf.count(char)
        buf = read_f(buf_size)

    f.close()
    return nlines

def simple_count(file):
    nbps = 0
    nlines = 0
    f = open(file)

    char = f.read(1) # get first character, > or @
    f.seek(0) # roll back

    record = False
    for line in f:
        if line[0] == char:
            nlines += 1
            record = True # record next row

        elif record:
            nbps += len(line) - 1
            record = False
            
    f.close()
    return nbps, nlines

def make_db(file):
    '''
    Make database. The type of the database will be determined based on diamond'd return code.
    '''
    if not os.path.isfile(file):
        logger.critical('File <{}> cannot be found. Please check input file (-i/--infile).'.format(file))
        sys.exit(2)

    logger.info('Building database of <{}> ...'.format(file))
    
    ## diamond
    subp = subprocess.run(['diamond', 'makedb', '--in', file, '--db', file, '--quiet'], check=False, stderr=subprocess.DEVNULL)
    
    if subp.returncode == 1: # auto choose dbtype by checking diamond's return code
        subprocess.run(['bwa', 'index', file], check=True, stderr=subprocess.DEVNULL)
        dbtype = 'nucl'
    else:
        if os.path.isfile(file + '.dmnd'):
            dbtype = 'prot'
        else:
            logger.critical('Cannot build diamond database of <{}>'.format(file))
            sys.exit(2)

    subprocess.run(['makeblastdb', '-in', file, '-dbtype', dbtype, '-out', file], check=True, stdout=subprocess.DEVNULL)
    logger.info('Finished.')

def run_make_db(options):
    make_db(vars(options).get('infile'))
