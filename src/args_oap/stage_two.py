import os
import subprocess
import sys
import pandas as pd

from glob import glob
from .utils import *
from .settings import logger, settings

class StageTwo:
    def __init__(self, options):
        self.__dict__.update(options)
        self._db = settings._sarg if self.database is None else self.database
        self._extracted = os.path.join(self.indir, 'extracted.fa')
        self._metadata = os.path.join(self.indir, 'metadata.txt')
        self._blastout = os.path.join(self.indir, 'blastout.txt')

        ## sanity check
        if not os.path.isdir(self.indir):
            logger.critical('Input folder <{}> does not exist.'.format(self.indir))
            sys.exit(2)
        else:
            if not os.path.isfile(os.path.join(self.indir, 'metadata.txt')):
                logger.critical('File "metadata.txt" does not exist in input folder <{}>.'.format(self.indir))
                sys.exit(2)
            if not os.path.isfile(os.path.join(self.indir, 'extracted.fa')):
                logger.critical('File "extracted.fa" does not exist in input folder <{}>.'.format(self.indir))
                sys.exit(2)

        ## database related
        if os.path.isfile(self._db + '.pdb'):
            self._dbtype = 'prot'
        elif os.path.isfile(self._db + '.ndb'):
            self._dbtype = 'nucl'
        else:
            logger.critical('Cannot load database <{}>.'.format(self._db))
            sys.exit(2)

        for structure in [self.structure1, self.structure2, self.structure3]:
            if structure is not None and not os.path.isfile(structure):
                logger.critical('Cannot load structure file <{}>.'.format(structure))
                sys.exit(2)

        structure_list = []
        if all([structure is None for structure in [self.structure1, self.structure2, self.structure3]]):
            if self.database is not None:
                logger.warning('No structures given for customized database.')

            for file, count in zip([settings._sarg_structure1, settings._sarg_structure2, settings._sarg_structure3], [1, 1/2, 1/3]):
                structure_list.append(pd.read_table(file, dtype = str).assign(count = count))
        else:
            if self.database is None:
                logger.warning('No database given for customized structures.')

            for file, count in zip([self.structure1, self.structure2, self.structure3], [1, 1/2, 1/3]):
                if file is not None:
                    structure_list.append(pd.read_table(file, dtype = str).assign(count = count))
        self.structures = pd.concat(structure_list, ignore_index=True)

    def extract_seqs(self):
        '''
        Extract target sequences using more stringent cutoffs & blast.
        '''
        logger.info('Processing <{}> ...'.format(self._extracted))
        nbps, nlines = simple_count(self._extracted)
        mt_mode = '1' if nbps / self.thread >= 2500000 else '0' 
        blast = 'blastx' if self._dbtype == 'prot' else 'blastn'

        logger.info('Extracting target sequences using BLAST ...'.format(self._extracted))
        logger.info('BLAST settings: {} bps, {} reads, {} threads, mt_mode {}.'.format(nbps, nlines, self.thread, mt_mode))
        
        cmd = [
            blast,
            '-db', self._db,
            '-query', self._extracted,
            '-out', self._blastout,
            '-outfmt', ' '.join(['6']+settings.cols),
            '-evalue', str(self.e),
            '-max_target_seqs', '5',
            '-num_threads', str(self.thread),
            '-mt_mode', mt_mode
            ]
        subprocess.run(cmd, check=True)

    def merge_files(self):
        '''
        Join extracted target sequences with metadata and structures. Aggregate according to levels (type/subtype/gene).
        '''
        logger.info('Merging files ...')
        try:
            df = pd.read_table(self._blastout, header=None, names=settings.cols)
        except:
            logger.critical('No target sequence detected, no further normalization will be made.')
            sys.exit(2)

        ## further filtering
        df = df[( df["pident"] >= self.id)  & ( df["evalue"] <= self.e)]
        df['qcov'] = df['length'] * 3 / df['qlen'] if self._dbtype == 'prot' else df['length'] / df['qlen'] # *3 for aa
        df = df[df['qcov'] >= self.qcov / 100]

        ## deduplicate
        df = df.sort_values(['qseqid', 'evalue', 'bitscore', 'length'], ascending=[True, True, False, False])
        df = df.drop_duplicates(subset='qseqid', keep='first')
        
        df['Sample'] = df['qseqid'].str.rsplit('@', 1).str.get(0).str.replace('_[12]$','', regex=True)
        df['scov'] = df['length'] / df['slen']

        ## merge blastout with metadata and structure
        levels = self.structures.columns[:-1]
        metadata = pd.read_table(self._metadata)
        if any(metadata['nCells']==0) or any(metadata['nReads']==0) or any(metadata['n16s']==0):
            logger.critical('Found zeros of reads/16s/cells in metadata file <{}>. Cannot normalize.'.format(self._metadata))
            sys.exit(2)

        if not all(x in set(self.structures[levels[0]]) for x in df['sseqid'].unique()):
            logger.warning('Not all extracted target sequences can be found in the structure file. Please check the database and/or the structure file.')
        
        result = pd.merge(df, self.structures, left_on = 'sseqid', right_on = levels[0], how = 'inner')
        if len(result)==0:
            logger.warning('No sequences remained after merging structure files.')

        result['scov'] = result['scov'] * result['count']
        for level in levels:
            for measure, normalizer, name in zip(['scov', 'scov', 'count'], 
                                                 ['nCells', 'n16s', 'nReads'],
                                                 ['normalized_cells', 'normalized_16s', 'ppm']):
                result['Sample'] = result['Sample'].astype(str)
                out = pd.merge(result.groupby([level, 'Sample'])[measure].sum(), metadata.set_index('Sample'), left_index = True, right_index=True, how ='outer').reset_index()
                out['value'] = out[measure]/out[normalizer]

                if measure == 'count':
                    out['value'] = out['value'] * 1e6

                out.set_index(level).pivot(columns = 'Sample', values = 'value').fillna(0).sort_index().to_csv(
                    os.path.join(self.indir, name + '.' + level + '.txt'), sep ='\t') 

                ## save unnormalized subject coverage or count
                if normalizer!= 'nCells':
                    out.set_index(level).pivot(columns = 'Sample', values = measure).fillna(0).sort_index().to_csv(
                        os.path.join(self.indir, 'unnormalized_' + measure + '.' + level + '.txt'), sep ='\t') 

    def run(self):
        self.extract_seqs()
        self.merge_files()
        logger.info('Finished.')

def run_stage_two(options):
    StageTwo(vars(options)).run()