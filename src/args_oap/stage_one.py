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

        self.ko30 = pd.read_table(settings._ko30_structure, header=None, names=['sseqid', 'ko30_name', 'ko30_length'])
        self.files = glob(os.path.join(self.indir, '*.' + self.format))

        ## sanity check
        if not os.path.isdir(self.indir):
            logger.critical('Input folder <{}> does not exist.'.format(self.indir))
            sys.exit(2)

        if len(self.files)==0:
            logger.critical('No file <*.{}> found in input folder <{}>. Please check the formats and change -f.'.format(self.format, self.indir))
            sys.exit(2)

        if os.path.isdir(self.outdir):
            logger.warning('Output folder <{}> exists, overwrite all files.'.format(self.outdir))
            for file in [self._extracted, self._metadata]:
                if os.path.isfile(file):
                    os.remove(file)
        else:
            os.makedirs(self.outdir, exist_ok=True)

        ## first time running, build database manually
        if not os.path.isfile(settings._gg85+'.ndb') or not os.path.isfile(settings._ko30+'.pdb') or not os.path.isfile(settings._sarg+'.pdb'):
            logger.info('Building databases ...')
            logger.disabled = True
            for file in [settings._gg85, settings._ko30, settings._sarg]:
                make_db(file)
            logger.disabled = False

        ## database related
        if os.path.isfile(self._db + '.pdb'):
            self._dbtype = 'prot'
        elif os.path.isfile(self._db + '.ndb'):
            self._dbtype = 'nucl'
        else:
            logger.critical('Cannot load database <{}>.'.format(self._db))
            sys.exit(2)

    def count_reads(self, file):
        '''
        Count reads (headers) of a fa/fq file.
        '''
        return buffer_count(file)

    def count_16s(self, file):
        '''
        Extract 16s (Greengenes 16S rRNA Database 85%) reads using bwa (pre-filtering) and blastn (post-filtering).
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
        try:
            df = pd.read_table(_tmp_16s_txt, header=None, names=settings.cols)
            df['scov'] = df['length'] / df['slen']
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in 16s.')
                df = df[~df['qseqid'].duplicated()]
                
            return df['scov'].sum()
        except:
            logger.warning("No 16s-like sequences found in file <{}>.".format(file))
            return 0

    def count_cells(self, file):
        '''
        Count Universial Essencial Single Copy Marker Genes (cells) using diamond.
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
            '--max-target-seqs', '1', 
            '--threads', str(self.thread), 
            '--quiet'], check=True)

        ## process blastx results, store subject coverage
        try:
            df = pd.read_table(_tmp_cells_txt, header=None, names=settings.cols)
            df = pd.merge(df, self.ko30, on='sseqid', how='left')
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in cells.')
                df = df[~df['qseqid'].duplicated()]

            return df.groupby('ko30_name').apply(lambda x: sum(x['length'] / x['ko30_length'])).sum()/30
        except:
            logger.warning("No cells-like sequences found in file <{}>.".format(file))
            return 0

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
                '--max-target-seqs', '1', 
                '--threads', str(self.thread), '--quiet'], check=True)
        else:
            subp = subprocess.run(['bwa', 'mem', '-t', str(self.thread), self._db, file], check=True, capture_output=True)
            subsubp = subprocess.run(['samtools', 'fasta', '-F4', '-F0x900', '-'], check=True, capture_output=True, input=subp.stdout)

            ## convert sam to tab for later usage, note that reads can be duplicated
            with open(_tmp_seqs_txt, 'w') as f:
                subprocess.run(['awk', 'BEGIN{RS=">";OFS="\t"}NR>1{print "#"$1,$2}'], check=True, input=subsubp.stdout, stdout=f)

        ## give a new header for each target sequences, merge all sequences to a single file
        with open(self._extracted, 'a') as f:
            df = pd.read_table(_tmp_seqs_txt, header=None, names=['qseqid', 'full_qseq'])
            if df['qseqid'].duplicated().sum()>0:
                logger.warning('Duplicated qseqid in target sequences.')
                df = df[~df['qseqid'].duplicated()]

            for i, line in enumerate(df['full_qseq']):
                f.write('>' + filename + '@' + str(i) + '\n' + line + '\n')

    def run(self):
        metadata = []
        for i, file in enumerate(sorted(self.files)):
            logger.info('Processing <{}> ({}/{}) ...'.format(file, i+1, len(self.files)))
            self.extract_seqs(file)

            ## skip 16s/cells calculation if necessary
            if not self.skip:
                nreads=self.count_reads(file)
                n16s=self.count_16s(file)
                ncells=self.count_cells(file)
                metadata.append([nreads,n16s,ncells, get_filename(file, self.format, drop=True)])

        if metadata:
            pd.DataFrame(metadata, columns=['nReads','n16s','nCells','Sample']).groupby('Sample').sum().to_csv(
                self._metadata, sep='\t')

        if not self.keep:
            for tmp in glob(os.path.join(self.outdir, '*.tmp')):
                os.remove(tmp)

        logger.info('Finished.')

def run_stage_one(options):
    StageOne(vars(options)).run()