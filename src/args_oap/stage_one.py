import os
import subprocess
import sys
import glob
import pandas as pd

from .utils import buffer_count, nucleotide_count, simple_count
from .settings import File, Setting, logger
from .make_db import make_db


class StageOne:
    def __init__(self, options):
        self.__dict__.update(options)
        self.setting = Setting(None, self.outdir)

        ## if given customized database then switch
        self.db = self.setting.sarg if self.database is None else self.database

        ## setup metadata
        self.ko30 = pd.read_table(self.setting.ko30_structure, header=None, names=['sseqid', 'ko30'])

        ## setup input/output files
        self.files = glob.glob(
            os.path.join(self.indir, '*.' + self.format)) + glob.glob(
            os.path.join(self.indir, '*.' + self.format + '.gz'))

        ## if indir does not exist
        if not os.path.isdir(self.indir):
            logger.critical(f'Input folder <{self.indir}> does not exist. Please check input folder (-i/--indir).')
            sys.exit(2)

        ## if glob no files in indir
        if len(self.files) == 0:
            logger.critical(f'No <*.{self.format}> or <*.{self.format}.gz> file found in input folder <{self.indir}>. Please check input folder (-i/--indir) or format (-f/--format).')
            sys.exit(2)

        ## if outdir same as indir
        if self.indir == self.outdir:
            logger.critical('Folder for input/output cannot be identical. Please change input folder (-i/--indir) or output folder (-o/--outdir).')
            sys.exit(2)

        ## if extracted.fa or metadata.txt exists in outdir, raise warning that they will be overwritten
        os.makedirs(self.outdir, exist_ok=True)
        if os.path.isfile(self.setting.extracted) or os.path.isfile(self.setting.metadata):
            logger.warning(f'Output folder <{self.outdir}> contains <extracted.fa> and/or <metadata.txt>, they/it will be overwritten.')
            for file in [
                self.setting.extracted,
                self.setting.metadata,
                self.setting.absolute_metadata,
                self.setting.metadata_with_absolute_quantification,
                os.path.join(self.outdir, '.absolute_quantification.tmp'),
                os.path.join(self.outdir, 'absolute_metadata.txt')
            ]:
                try:
                    os.remove(file)
                except OSError:
                    pass

        ## first time running, build database manually
        if not (
            os.path.isfile(f'{self.setting.gg85}.ndb') or
            os.path.isfile(f'{self.setting.ko30}.pdb') or
            os.path.isfile(f'{self.setting.sarg}.pdb')
        ):
            logger.info('Building databases ...')
            logger.disabled = True  # skip logging
            for file in [self.setting.gg85, self.setting.ko30, self.setting.sarg]:
                make_db(file)
            logger.disabled = False

        ## check database type of customized database
        if os.path.isfile(f'{self.db}.pdb'):
            self.dbtype = 'prot'
        elif os.path.isfile(f'{self.db}.ndb'):
            self.dbtype = 'nucl'
        else:
            logger.critical(f'Cannot find database <{self.db}>. Please run <make_db> first or check database (--database).')
            sys.exit(2)

    def count_16s(self, file):
        '''
        Count 16S (GreenGenes 16S rRNA Database 85%) copy number using bwa (pre-filtering) and blastn (post-filtering).
        '''
        ## pre-filtering using bwa
        subprocess.run([
            'bwa', 'mem',
            '-t', str(self.thread),
            '-o', file.tmp_16s_sam,
            self.setting.gg85, file.file], check=True, stderr=subprocess.DEVNULL)

        ## convert sam to fasta for later usage, note that reads can be duplicated
        with open(file.tmp_16s_fa, 'w') as f:
            subprocess.run([
                'samtools',
                'fasta',
                '-F', '2308',
                file.tmp_16s_sam], check=True, stderr=subprocess.DEVNULL, stdout=f)

        ## post-filter using blastn
        ## switch mt_mode if too little queries or too many threads, blast raises error if <2,500,000 bases per thread
        mt_mode = '1' if simple_count(file.tmp_16s_fa)[0] / self.thread >= 2500000 else '0'
        subprocess.run([
            'blastn',
            '-db', self.setting.gg85,
            '-query', file.tmp_16s_fa,
            '-out', file.tmp_16s_txt,
            '-outfmt', ' '.join(['6'] + self.setting.columns),
            '-evalue', str(self.e1),
            '-max_hsps', '1',
            '-max_target_seqs', '1',
            '-mt_mode', mt_mode,
            '-num_threads', str(self.thread)], check=True, stderr=subprocess.DEVNULL)

        ## process blastn results, store subject cover
        df = pd.read_table(file.tmp_16s_txt, header=None, names=self.setting.columns)
        if len(df) == 0:
            logger.warning(f'No 16S-like sequences found in file <{file.file}>.')
        else:
            if df['qseqid'].duplicated().sum() > 0:
                logger.warning('Duplicated sequences in 16S copy number calculation.')
                df = df[~df['qseqid'].duplicated()]

            return (df['length'] / df['slen']).sum()

    def count_cells(self, file):
        '''
        Count Essential Single Copy Marker Genes (cell number) using diamond.
        '''
        ## filter using diamond
        subprocess.run([
            'diamond', 'blastx',
            '--db', f'{self.setting.ko30}.dmnd',
            '--query', file.file,
            '--out', file.tmp_cells_txt,
            '--outfmt', '6'] + self.setting.columns + [
            '--evalue', str(self.e2),
            '--id', str(self.id),
            '--query-cover', str(self.qcov),
            '--max-hsps', '1',
            '--max-target-seqs', '1',
            '--threads', str(self.thread),
            '--quiet'], check=True, stderr=subprocess.DEVNULL)

        ## process blastx results, store subject coverage
        df = pd.merge(pd.read_table(file.tmp_cells_txt, header=None, names=self.setting.columns), self.ko30, on='sseqid', how='left')
        if len(df) == 0:
            logger.warning(f'No marker-like sequences found in file <{file.file}>.')
        else:
            if df['qseqid'].duplicated().sum() > 0:
                logger.warning('Duplicated sequences in cell number calculation.')
                df = df[~df['qseqid'].duplicated()]

            return df.groupby('ko30').apply(lambda x: sum(x['length'] / x['slen'])).sum() / 30

    def extract_seqs(self, file):
        '''
        Prefilter ARGs, or other target sequences.
        '''
        ## diamond or bwa
        if self.dbtype == 'prot':
            subprocess.run([
                'diamond', 'blastx',
                '--db', f'{self.db}.dmnd',
                '--query', file.file,
                '--out', file.tmp_seqs_txt,
                '--outfmt', '6', 'qseqid', 'full_qseq',
                '--evalue', '10',
                '--id', '60',
                '--query-cover', '15',
                '--max-hsps', '1',
                '--max-target-seqs', '1',
                '--threads', str(self.thread),
                '--quiet'], check=True, stderr=subprocess.DEVNULL)
        else:
            subprocess.run([
                'bwa', 'mem',
                '-t', str(self.thread),
                '-o', file.tmp_seqs_sam,
                self.db, file.file], check=True, stderr=subprocess.DEVNULL)

            with open(file.tmp_seqs_fa, 'w') as f:
                subprocess.run([
                    'samtools',
                    'fasta',
                    '-F', '2308',
                    file.tmp_seqs_sam], check=True, stderr=subprocess.DEVNULL, stdout=f)

            ## convert fa to tab to make results consistent
            with open(file.tmp_seqs_txt, 'w') as w:
                with open(file.tmp_seqs_fa) as f:
                    for line in f:
                        if line[0] == '>':
                            qseqid = line[1:].rstrip().split(' ')[0]
                        else:
                            w.write(f'{qseqid}\t{line}')

        ## give a new header for each target sequences, merge all sequences to a single file
        with open(self.setting.extracted, 'a') as f:
            df = pd.read_table(file.tmp_seqs_txt, header=None, names=['qseqid', 'full_qseq'])
            if df['qseqid'].duplicated().sum() > 0:
                logger.warning('Duplicated sequences in sequence extraction.')
                df = df[~df['qseqid'].duplicated()]

            cnt = 0
            for qseqid, full_qseq in zip(df['qseqid'], df['full_qseq']):
                f.write(f'>{file.sample_name}@{file.file_name}@{cnt}@{qseqid}\n{full_qseq}\n')
                cnt += 1

    def run(self):
        metadata = []
        for i, file in enumerate(sorted(self.files)):
            logger.info(f'Processing <{file}> ({i + 1}/{len(self.files)}) ...')
            file = File(file, self.outdir, self.format)
            try:
                self.extract_seqs(file)

                ## skip 16S/cells calculation if necessary
                if not self.skip:
                    nread = buffer_count(file.file)
                    n16S = self.count_16s(file)
                    ncell = self.count_cells(file)
                    nucleotide = nucleotide_count(file.file)
                    metadata.append([nread, nucleotide['A'], nucleotide['T'], nucleotide['C'], nucleotide['G'], n16S, ncell, file.sample_name])
            except subprocess.CalledProcessError:
                logger.warning(f'Something is wrong with <{file.file}>, skip.')

            if not self.keep:
                for tmp in glob.glob(os.path.join(self.outdir, '*.tmp')):
                    os.remove(tmp)

        nucleotide_columns = ['An', 'Tn', 'Cn', 'Gn']
        metadata = pd.DataFrame(metadata, columns=['nRead'] + nucleotide_columns + ['n16S', 'nCell', 'sample']).groupby('sample').sum()

        ## optional metadata for absolute quantification
        meta_columns = ['SampleVolume_mL', 'ElutionVolume_uL', 'DNAConcentration_ng/uL']
        if self.absquantmeta is not None and os.path.isfile(self.absquantmeta):
            sample_metadata = pd.read_table(self.absquantmeta, dtype={'sample': str})
            required_columns = ['sample'] + meta_columns
            invalid = []
            if 'sample' not in sample_metadata.columns:
                logger.warning(f'File <{self.absquantmeta}> has no sample column. Skip absolute quantification for all samples.')
                sample_metadata = pd.DataFrame(columns=required_columns)
            else:
                for column in meta_columns:
                    if column not in sample_metadata.columns:
                        sample_metadata[column] = pd.NA
                duplicated = sample_metadata.loc[sample_metadata['sample'].duplicated(), 'sample'].tolist()
                if duplicated:
                    logger.warning(f'Duplicated absolute-quantification metadata found for sample(s): {", ".join(duplicated)}. Use the first row.')
                    sample_metadata = sample_metadata.drop_duplicates(subset='sample', keep='first')
                sample_metadata[meta_columns] = sample_metadata[meta_columns].apply(pd.to_numeric, errors='coerce')
                valid = sample_metadata[meta_columns].notnull().all(axis=1) & (sample_metadata[meta_columns] > 0).all(axis=1)
                invalid = sample_metadata.loc[~valid, 'sample'].dropna().tolist()
                if invalid:
                    logger.warning(f'Invalid or incomplete absolute-quantification metadata for sample(s): {", ".join(invalid)}. Skip absolute quantification for them.')
                sample_metadata = sample_metadata[valid]

            absolute_metadata = metadata.reset_index().merge(sample_metadata[required_columns], on='sample', how='inner').set_index('sample')
            missing = [sample for sample in metadata.index.difference(absolute_metadata.index) if sample not in invalid]
            if missing:
                logger.warning(f'No absolute-quantification metadata found for sample(s): {", ".join(missing)}. Absolute quantification will not be calculated for them.')
            if len(absolute_metadata) > 0:
                absolute_metadata[nucleotide_columns + meta_columns].to_csv(self.setting.absolute_metadata, sep='\t')
        else:
            if self.absquantmeta is not None:
                logger.warning(f'Absolute-quantification metadata file <{self.absquantmeta}> does not exist. Skip absolute quantification.')

        metadata = metadata.drop(columns=nucleotide_columns)

        metadata.to_csv(self.setting.metadata, sep='\t')
        logger.info('Finished.')


def run_stage_one(options):
    StageOne(vars(options)).run()
