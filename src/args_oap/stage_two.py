import os
import subprocess
import sys
import pandas as pd

from .utils import simple_count
from .settings import Setting, logger


class StageTwo:
    def __init__(self, options):
        self.__dict__.update(options)
        if self.outdir is None:
            self.outdir = self.indir
        else:
            os.makedirs(self.outdir, exist_ok=True)
        self.setting = Setting(self.indir, self.outdir)

        ## if given customized database then switch
        self.db = self.setting.sarg if self.database is None else self.database

        ## in indir does not exists
        if not os.path.isdir(self.indir):
            logger.critical(f'Input folder <{self.indir}> does not exist. Please check input folder (-i/--indir).')
            sys.exit(2)
        else:
            ## if metadata.txt missing
            if not os.path.isfile(self.setting.metadata):
                logger.critical(f'File <metadata.txt> does not exist in <{self.indir}>. Please check input folder (-i/--indir).')
                sys.exit(2)

            ## if extracted.fa missing
            if not os.path.isfile(self.setting.extracted):
                logger.critical(f'File <extracted.fa> does not exist in <{self.indir}>. Please check input folder (-i/--indir).')
                sys.exit(2)

        ## bypass blast
        if self.blastout is not None:
            if not os.path.isfile(self.blastout):
                logger.critical(f'File <{self.blastout}> does not exist. Please check BLAST output file (--blastout)')
                sys.exist(2)
            else:
                logger.info(f'BLAST output file <{self.blastout}> given, skip BLAST')
                self.setting.blastout = self.blastout
        else:
            ## if blastout.txt exists in outdir, raise warning that it will be overwritten
            if (
                os.path.isfile(self.setting.blastout) or
                os.path.isfile(self.setting.extracted_filtered) or
                os.path.isfile(self.setting.blastout_filtered)
            ):
                logger.warning(f'Output folder <{self.outdir}> contains <blastout.txt> and/or <*.filtered.*>, they/it will be overwritten.')
            for file in [self.setting.blastout, self.setting.extracted_filtered, self.setting.blastout_filtered]:
                try:
                    os.remove(file)
                except OSError:
                    pass

        ## if zero entries in metadata then cannot normalize
        self.metadata = pd.read_table(self.setting.metadata, dtype={'sample': str})
        if (self.metadata[['nRead', 'n16S', 'nCell']] == 0).any(axis=None):
            logger.warning(f'Found zero reads/16s/cells in metadata file <{self.setting.metadata}>. These samples will be ignored.')
            self.metadata = self.metadata[~(self.metadata == 0).any(axis=1)]

        if len(self.metadata) == 0:
            logger.critical(f'No valid sample in metadata file <{self.setting.metadata}>. No further normalization will be made.')
            sys.exit(2)

        ## database check
        if os.path.isfile(f'{self.db}.pdb'):
            self.dbtype = 'prot'
        elif os.path.isfile(f'{self.db}.ndb'):
            self.dbtype = 'nucl'
        else:
            logger.critical(f'Cannot find database <{self.db}>. Please run <makedb> first or check database (--database)')
            sys.exit(2)

        ## structures check
        for structure in [self.structure1, self.structure2, self.structure3]:
            if structure is not None and not os.path.isfile(structure):
                logger.critical(f'Cannot load structure file <{structure}>. Please check structures (--strucutre1/--strucutre2/--strucutre3)')
                sys.exit(2)

        structure_list = []
        if all(structure is None for structure in [self.structure1, self.structure2, self.structure3]):
            if self.database is not None:
                logger.warning('No structures given for customized database. Use default SARG structures')

            for file, count in zip([self.setting.sarg_structure1, self.setting.sarg_structure2, self.setting.sarg_structure3], [1, 1 / 2, 1 / 3]):
                structure_list.append(pd.read_table(file, dtype=str).assign(count=count))

        else:
            if self.database is None:
                logger.warning('No database given for customized structures. Use default SARG database')

            for file, count in zip([self.structure1, self.structure2, self.structure3], [1, 1 / 2, 1 / 3]):
                if file is not None:
                    structure_list.append(pd.read_table(file, dtype=str).assign(count=count))

        self.structures = pd.concat(structure_list, ignore_index=True)

    def extract_seqs(self):
        '''
        Extract target sequences using more stringent cutoffs & blast.
        '''
        logger.info(f'Processing <{self.setting.extracted}> ...')
        nbps, nlines = simple_count(self.setting.extracted)
        mt_mode = '1' if nbps / self.thread >= 2500000 else '0'
        blast_mode = 'blastx' if self.dbtype == 'prot' else 'blastn'

        logger.info('Extracting target sequences using BLAST ...')
        logger.info(f'BLAST settings: {nbps} bps, {nlines} reads, {self.thread} threads, mt_mode {mt_mode}.')

        subprocess.run([
            blast_mode,
            '-db', self.db,
            '-query', self.setting.extracted,
            '-out', self.setting.blastout,
            '-outfmt', ' '.join(['6'] + self.setting.columns),
            '-evalue', str(self.e),
            '-max_target_seqs', '5',
            '-num_threads', str(self.thread),
            '-mt_mode', mt_mode], check=True)

    def merge_files(self):
        '''
        Join extracted target sequences with metadata and structures. Aggregate according to levels (type/subtype/gene).
        '''
        logger.info('Merging files ...')
        df = pd.read_table(self.setting.blastout, header=None, names=self.setting.columns, dtype={'sseqid': str})
        if len(df) == 0:
            logger.critical(f'No target sequence detected in file <{self.setting.extracted}>, no further normalization will be made.')
            sys.exit(2)

        if df.isnull().any(axis=None):
            logger.critical(f'File <{self.blastout}> cannot be parsed. Please check BLAST output file (--blastout).')
            sys.exit(2)

        ## further filtering
        self.qcov = self.qcov / 100 / 3 if self.dbtype == 'prot' else self.qcov / 100  # for prot need to lower qcov as qlen is in bp
        self.length = self.length if self.dbtype == 'prot' else self.length * 3  # for nucl need to convert aa cut to bp

        df['qcov'] = df['length'] / df['qlen']
        df = df[(df['pident'] >= self.id) & (df['evalue'] <= self.e) & (df['length'] >= self.length) & (df['qcov'] >= self.qcov)]

        ## deduplicate
        df = df.sort_values(['qseqid', 'evalue', 'bitscore', 'length'], ascending=[True, True, False, False])
        df = df.drop_duplicates(subset='qseqid', keep='first')

        df['sample'] = df['qseqid'].str.rsplit(pat='@', n=3).str.get(0)
        df['scov'] = df['length'] / df['slen']

        ## merge blastout with metadata and structure
        levels = self.structures.columns[:-1]
        if not all(x in set(self.structures[levels[0]]) for x in df['sseqid'].unique()):
            logger.warning('Not all extracted target sequences can be found in structure files. Please check the database and/or the structure files.')

        result = pd.merge(df, self.structures, left_on='sseqid', right_on=levels[0], how='inner')
        if len(result) == 0:
            logger.critical('No target sequence remained after merging structure files, no further normalization will be made.')
            sys.exit(2)

        ## compute normalization metrics
        result['copy'] = result['scov'] * result['count']
        result['rpk'] = result['count'] / (result['slen'] / 1000)
        result['rpksum'] = result.groupby('sample')['rpk'].transform('sum')

        ## save filtered results
        result.to_csv(self.setting.blastout_filtered, sep='\t', index=False)
        filtered_qseqid = set(result.qseqid)

        with open(self.setting.extracted_filtered, 'w') as w:
            with open(self.setting.extracted) as f:
                record = False
                for line in f:
                    if line[0] == '>':
                        if line[1:].rstrip().split(' ')[0] in filtered_qseqid:
                            w.write(line)
                            record = True
                    elif record:
                        w.write(line)
                        record = False

        for level in levels:
            for numerator, denominator, metric in zip(
                ['copy', 'copy', 'count', 'rpk', 'rpk'],
                ['nCell', 'n16S', 'nRead', 'nRead', 'rpksum'],
                ['normalized_cell', 'normalized_16S', 'ppm', 'rpkm', 'tpm']
            ):

                out = pd.merge(self.metadata, result, on='sample')
                out['value'] = out[numerator] / out[denominator]
                if metric in {'ppm', 'rpkm', 'tpm'}:
                    out['value'] = out['value'] * 1e6

                ## aggregate over levels
                out_agg = out.groupby([level, 'sample'])['value'].sum().unstack().fillna(0)
                out_agg.sort_index().to_csv(os.path.join(self.outdir, f'{metric}.{level}.txt'), sep='\t')

                ## save unnormalized subject coverage or count
                if metric in {'normalized_cell', 'ppm'}:
                    out_agg = out.groupby([level, 'sample'])[numerator].sum().unstack().fillna(0)
                    out_agg.sort_index().to_csv(os.path.join(self.outdir, f'unnormalized_{numerator}.{level}.txt'), sep='\t')

    def run(self):
        if self.blastout is None:
            self.extract_seqs()
        self.merge_files()
        logger.info('Finished.')


def run_stage_two(options):
    StageTwo(vars(options)).run()
