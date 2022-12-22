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

        ## in indir does not exists
        if not os.path.isdir(self.indir):
            logger.critical('Input folder <{}> does not exist. Please check input folder (-i/--indir).'.format(self.indir))
            sys.exit(2)
        else:
            ## if metadata.txt missing
            if not os.path.isfile(os.path.join(self.indir, 'metadata.txt')):
                logger.critical('File "metadata.txt" does not exist in input folder <{}>. Please run <stage_one> first or check input folder (-i/--indir).'.format(self.indir))
                sys.exit(2)

            ## if extracted.fa missing
            if not os.path.isfile(os.path.join(self.indir, 'extracted.fa')):
                logger.critical('File "extracted.fa" does not exist in input folder <{}>. Please run <stage_one> first or check input folder (-i/--indir).'.format(self.indir))
                sys.exit(2)

        ## if outdir is not given, use indir instead
        self.outdir = self.indir if self.outdir is None else self.outdir
        self._blastout = os.path.join(self.outdir, 'blastout.txt')
        os.makedirs(self.outdir, exist_ok=True)

        ## bypass blast
        if self.blastout is not None:
            if not os.path.isfile(self.blastout):
                logger.critical('BLAST output file <{}> does not exist. Please check BLAST output file (--blastout)'.format(self.blastout))
                sys.exist(2)
            else:
                logger.info('BLAST output file <{}> given, skip BLAST'.format(self.blastout))
                self._blastout = self.blastout
        else:
            ## if blastout.txt exists in outdir, raise warning that it will be overwritten
            if os.path.isfile(self._blastout):
                logger.warning('Output folder <{}> contains <blastout.txt>, it will be overwritten.'.format(self.outdir))
                os.remove(self._blastout)

        ## metadata if zero entries then cannot normalise
        self.metadata = pd.read_table(self._metadata, index_col='Sample')
        if (self.metadata[['nRead', 'n16S', 'nCell']]==0).any(axis=None):
            logger.warning('Found zero reads/16s/cells in some samples in metadata file <{}>. These samples will be ignored.'.format(self._metadata))
            self.metadata = self.metadata[~(self.metadata==0).any(axis=1)]

        if len(self.metadata)==0:
            logger.critical('No valid sample in metadata file <{}>, no further normalization will be made.'.format(self._metadata))
            sys.exit(2)

        ## database check
        if os.path.isfile(self._db + '.pdb'):
            self._dbtype = 'prot'
        elif os.path.isfile(self._db + '.ndb'):
            self._dbtype = 'nucl'
        else:
            logger.critical('Cannot find database <{}>. Please run <make_db> first or check database (--database)'.format(self._db))
            sys.exit(2)

        ## structures check
        for structure in [self.structure1, self.structure2, self.structure3]:
            if structure is not None and not os.path.isfile(structure):
                logger.critical('Cannot load structure file <{}>. Please check structures (--strucutre1/--strucutre2/--strucutre3)'.format(structure))
                sys.exit(2)

        structure_list = []
        if all([structure is None for structure in [self.structure1, self.structure2, self.structure3]]):
            if self.database is not None:
                logger.warning('No structures given for customized database. Use default SARG structures')

            for file, count in zip([settings._sarg_structure1, settings._sarg_structure2, settings._sarg_structure3], [1, 1/2, 1/3]):
                structure_list.append(pd.read_table(file, dtype = str).assign(count = count))
        else:
            if self.database is None:
                logger.warning('No database given for customized structures. Use default SARG database')

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
        df = pd.read_table(self._blastout, header=None, names=settings.cols)
        if len(df)==0:
            logger.critical('No target sequence detected in <{}>, no further normalization will be made.'.format(self._extracted))
            sys.exit(2)

        if df.isnull().any(axis=None):
            logger.critical('BLAST output file <{}> cannot be parsed. Please check BLAST output file (--blastout).'.format(self.blastout))
            sys.exit(2)

        ## further filtering
        self.qcov = self.qcov/3/100 if self._dbtype == 'prot' else self.qcov/100 # for prot need to lower qcov as qlen is in bp
        self.length = self.length if self._dbtype == 'prot' else self.length * 3 # for nucl need to convert aa cut to bp

        df = df[(df["pident"] >= self.id)  & ( df["evalue"] <= self.e) & (df['length'] >= self.length)]
        df['qcov'] = df['length'] / df['qlen']
        df = df[(df['qcov'] >= self.qcov)]

        ## deduplicate
        df = df.sort_values(['qseqid', 'evalue', 'bitscore', 'length'], ascending=[True, True, False, False])
        df = df.drop_duplicates(subset='qseqid', keep='first')
        
        df['Sample'] = df['qseqid'].str.rsplit(pat='@', n=2).str.get(0).str.replace('_[12]$','', regex=True)
        df['scov'] = df['length'] / df['slen']

        ## merge blastout with metadata and structure
        levels = self.structures.columns[:-1]
        if not all(x in set(self.structures[levels[0]]) for x in df['sseqid'].unique()):
            logger.warning('Not all extracted target sequences can be found in the structure file. Please check the database and/or the structure files.')
        
        result = pd.merge(df, self.structures, left_on = 'sseqid', right_on = levels[0], how = 'inner')
        if len(result)==0:
            logger.critical('No target sequence remained after merging structure files, no further normalization will be made.')
            sys.exit(2)

        result['copy'] = result['scov'] * result['count']
        result['rpk'] = result['count']/ (result['slen'] / 1000)

        for level in levels:
            for measure, normalizer, name in zip(['copy', 'copy', 'count', 'rpk', 'rpk'], 
                                                 ['nCell', 'n16S', 'nRead', 'nRead', 'nRead'],
                                                 ['normalized_cell', 'normalized_16S', 'ppm', 'rpkm', 'tpm']):

                result['Sample'] = result['Sample'].astype(str)
                out = pd.merge(self.metadata, result.groupby([level, 'Sample'])[measure].sum(), left_index = True, right_index=True, how ='left').reset_index()
                out['value'] = out[measure]/out[normalizer]

                if name in {'ppm', 'rpkm'}:
                    out['value'] = out['value'] * 1e6

                if name in {'tpm'}:
                    out = pd.merge(out, out.groupby('Sample')['value'].sum().reset_index(), on='Sample')
                    out['value'] = out['value_x'] / out['value_y']

                out.set_index(level).pivot(columns = 'Sample', values = 'value').fillna(0).sort_index().to_csv(
                    os.path.join(self.outdir, name + '.' + level + '.txt'), sep ='\t') 

                ## save unnormalized subject coverage or count
                if name in {'normalized_cell', 'ppm'}:
                    out.set_index(level).pivot(columns = 'Sample', values = measure).fillna(0).sort_index().to_csv(
                        os.path.join(self.outdir, 'unnormalized_' + measure + '.' + level + '.txt'), sep ='\t')

    def run(self):
        if self.blastout is None:
            self.extract_seqs()
        self.merge_files()
        logger.info('Finished.')

def run_stage_two(options):
    StageTwo(vars(options)).run()