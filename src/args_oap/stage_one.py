import os
import subprocess
import sys
import pandas as pd

from glob import glob
from .utils import *
from .settings import logger, settings

class StageOne:
    def __init__(self, options):
        self.__dict__.update(options)
        self._db = settings._sarg if self.database is None else self.database
        self._extracted = os.path.join(self.outdir, 'extracted.fa')
        self._metadata = os.path.join(self.outdir, 'metadata.txt')

        self.ko30 = pd.read_table(settings._ko30_structure, header=None, names=['sseqid', 'ko30'])
        self.files = glob(os.path.join(self.indir, '*.' + self.format)) + glob(os.path.join(self.indir, '*.' + self.format + '.gz'))

        ## if indir does not exist
        if not os.path.isdir(self.indir):
            logger.critical('Input folder <{}> does not exist. Please check input folder (-i/--indir).'.format(self.indir))
            sys.exit(2)
        else:
            ## if glob no files in indir
            if len(self.files)==0:
                logger.critical('No <*.{}> or <*.{}.gz> file found in input folder <{}>. Please check input folder (-i/--indir) or format (-f/--format).'.format(self.format, self.format, self.indir))
                sys.exit(2)

        ## if outdir same as indir
        if self.indir == self.outdir:
            logger.critical('Folder for input/output cannot be identical. Please check input folder (-i/--indir) or output folder (-o/--outdir).')
            sys.exit(2)

        ## if extracted.fa or metadata.txt exists in outdir, raise warning that they will be overwritten
        os.makedirs(self.outdir, exist_ok=True)
        if os.path.isfile(self._extracted) or os.path.isfile(self._metadata):
            logger.warning('Output folder <{}> contains <extracted.fa> and/or <metadata.txt>, they/it will be overwritten.'.format(self.outdir))
            for file in [self._extracted, self._metadata]:
                try:
                    os.remove(file)
                except OSError:
                    pass

        ## first time running, build database manually
        if not os.path.isfile(settings._gg85+'.ndb') or not os.path.isfile(settings._ko30+'.pdb') or not os.path.isfile(settings._sarg+'.pdb'):
            logger.info('Building databases ...')
            logger.disabled = True # skip logging
            for file in [settings._gg85, settings._ko30, settings._sarg]:
                make_db(file)
            logger.disabled = False

        ## database check
        if os.path.isfile(self._db + '.pdb'):
            self._dbtype = 'prot'
        elif os.path.isfile(self._db + '.ndb'):
            self._dbtype = 'nucl'
        else:
            logger.critical('Cannot find database <{}>. Please run <make_db> first or check database (--database)'.format(self._db))
            sys.exit(2)

    def count_16s(self, file):
        '''
        Extract 16S (GreenGenes 16S rRNA Database 85%) reads using bwa (pre-filtering) and blastn (post-filtering).
        '''
        filename = get_filename(file, self.format)
        _tmp_16s_fa = os.path.join(self.outdir, filename + '.16s.fa.tmp')
        _tmp_16s_txt = os.path.join(self.outdir, filename + '.16s.txt.tmp')

        ## pre-filtering using bwa
        subp = subprocess.run(['bwa', 'mem', '-t', str(self.thread), settings._gg85, file], check=True, capture_output=True)

        ## convert sam to fasta for later usage, note that reads can be duplicated
        with open(_tmp_16s_fa, 'w') as f:
            subprocess.run(['samtools', 'fasta', '-F4', '-F0x900', '-'], check=True, input=subp.stdout, stdout=f, stderr=subprocess.DEVNULL)

        ## post-filter using blastn
        ## switch mt_mode if too little queries or too many threads, blast raises error if <2,500,000 bases per thread
        mt_mode = '1' if simple_count(_tmp_16s_fa)[0] / self.thread >= 2500000 else '0' 
        subprocess.run([
            'blastn', 
            '-db', settings._gg85, 
            '-query', _tmp_16s_fa, 
            '-out', _tmp_16s_txt, 
            '-outfmt', ' '.join(['6']+settings.cols),
            '-evalue', str(self.e1), 
            '-max_hsps', '1', 
            '-max_target_seqs', '1',
            '-mt_mode', mt_mode, 
            '-num_threads', str(self.thread)], check=True, stderr=subprocess.DEVNULL)

        ## process blastn results, store subject cover
        df = pd.read_table(_tmp_16s_txt, header=None, names=settings.cols)
        if len(df)==0:
            logger.warning("No 16S-like sequences found in file <{}>.".format(file))
            return 0
        else:
            df['scov'] = df['length'] / df['slen']
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in 16S.')
                df = df[~df['qseqid'].duplicated()]
            return df['scov'].sum()

    def count_cells(self, file):
        '''
        Count Essential Single Copy Marker Genes (cells) using diamond.
        '''
        filename = get_filename(file, self.format)
        _tmp_cells_txt = os.path.join(self.outdir, filename + '.cells.txt.tmp')

        ## filter using diamond
        subprocess.run([
            'diamond', 'blastx', 
            '--db', settings._ko30 + '.dmnd', 
            '--query', file, 
            '--out', _tmp_cells_txt, 
            '--outfmt', '6'] + settings.cols + [
            '--evalue', str(self.e2), 
            '--id', str(self.id), 
            '--query-cover', str(self.qcov), 
            '--max-hsps', '1',
            '--max-target-seqs', '1', 
            '--threads', str(self.thread), 
            '--quiet'], check=True)

        ## process blastx results, store subject coverage
        df = pd.read_table(_tmp_cells_txt, header=None, names=settings.cols)
        if len(df)==0:
            logger.warning("No cells-like sequences found in file <{}>.".format(file))
            return 0
        else:
            df = pd.merge(df, self.ko30, on='sseqid', how='left')
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in cells {}.'.format(_tmp_cells_txt))
                df = df[~df['qseqid'].duplicated()]
            return df.groupby('ko30').apply(lambda x: sum(x['length'] / x['slen'])).sum()/30

    def extract_seqs(self, file):
        '''
        Prefilter ARGs, or other target sequences.
        '''
        filename = get_filename(file, self.format)
        _tmp_seqs_txt = os.path.join(self.outdir, filename + '.seqs.txt.tmp')

        ## diamond or bwa
        if self._dbtype == 'prot':
            subprocess.run([
                'diamond', 'blastx', 
                '--db', self._db + '.dmnd', 
                '--query', file, 
                '--out', _tmp_seqs_txt, 
                '--outfmt', '6', 'qseqid', 'full_qseq', 
                '--evalue', '10', 
                '--id', '60', 
                '--query-cover', '15', 
                '--max-hsps', '1',
                '--max-target-seqs', '1', 
                '--threads', str(self.thread), '--quiet'], check=True)
        else:
            subp = subprocess.run(['bwa', 'mem', '-t', str(self.thread), self._db, file], check=True, capture_output=True)
            subsubp = subprocess.run(['samtools', 'fasta', '-F4', '-F0x900', '-'], check=True, capture_output=True, input=subp.stdout)

            ## convert sam to tab for later usage
            with open(_tmp_seqs_txt, 'w') as f:
                subprocess.run(['awk', 'BEGIN{RS=">";OFS="\t"}NR>1{print "#"$1,$2}'], check=True, input=subsubp.stdout, stdout=f)

        ## give a new header for each target sequences, merge all sequences to a single file
        with open(self._extracted, 'a') as f:
            df = pd.read_table(_tmp_seqs_txt, header=None, names=['qseqid', 'full_qseq'])
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in target sequences.')
                df = df[~df['qseqid'].duplicated()]

            cnt = 0
            for qseqid, full_qseq in zip(df['qseqid'], df['full_qseq']):
                f.write('>' + filename + '@' + str(cnt) + '@' + qseqid + '\n' + full_qseq + '\n')
                cnt += 1

    def run(self):
        metadata = []
        for i, file in enumerate(sorted(self.files)):
            logger.info('Processing <{}> ({}/{}) ...'.format(file, i+1, len(self.files)))
            self.extract_seqs(file)

            ## skip 16S/cells calculation if necessary
            if not self.skip:
                nread=buffer_count(file)
                n16S=self.count_16s(file)
                ncell=self.count_cells(file)
                metadata.append([nread, n16S, ncell, get_filename(file, self.format, drop=True)])

            if not self.keep:
                for tmp in glob(os.path.join(self.outdir, '*.tmp')):
                    os.remove(tmp)
                    
        if metadata:
            pd.DataFrame(metadata, columns=['nRead','n16S','nCell','Sample']).groupby('Sample').sum().to_csv(
                self._metadata, sep='\t')

        logger.info('Finished.')

def run_stage_one(options):
    StageOne(vars(options)).run()