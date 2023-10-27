import os
import re
import logging

from dataclasses import dataclass

## setup logger
logging.basicConfig(
    level="INFO",
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")

## bold format
BOLD = "\033[1m"
RESET = "\033[0m"

## add color
logging.addLevelName(
    logging.WARNING,
    BOLD + "\x1b[33;20m%s\033[1;0m" % logging.getLevelName(logging.WARNING))

logging.addLevelName(
    logging.CRITICAL,
    BOLD + "\x1b[31;20m%s\033[1;0m" % logging.getLevelName(logging.CRITICAL))

logger = logging.getLogger(__name__)


## setup file format
@dataclass
class File:
    file: str
    outdir: str
    format: str

    @property
    def file_name(self) -> str:
        return re.sub(rf'\.{self.format}(.gz)?$', '', os.path.basename(self.file))

    @property
    def sample_name(self) -> str:
        return re.sub(r'(_R1|_R2|_1|_2|_fwd|_rev)$', '', self.file_name)

    @property
    def tmp_16s_fa(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.16s.fa.tmp')

    @property
    def tmp_16s_txt(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.16s.txt.tmp')

    @property
    def tmp_16s_sam(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.16s.sam.tmp')

    @property
    def tmp_cells_txt(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.cells.txt.tmp')

    @property
    def tmp_seqs_fa(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.seqs.fa.tmp')

    @property
    def tmp_seqs_txt(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.seqs.txt.tmp')

    @property
    def tmp_seqs_sam(self) -> str:
        return os.path.join(self.outdir, self.file_name + '.seqs.sam.tmp')


## setup database
@dataclass
class Setting:
    indir: str
    outdir: str
    db: str = os.path.join(os.path.dirname(__file__), 'db')

    @property
    def sarg(self) -> str:
        return os.path.join(self.db, 'sarg.fasta')

    @property
    def sarg_structure1(self) -> str:
        return os.path.join(self.db, 'single-component_structure.txt')

    @property
    def sarg_structure2(self) -> str:
        return os.path.join(self.db, 'two-component_structure.txt')

    @property
    def sarg_structure3(self) -> str:
        return os.path.join(self.db, 'multi-component_structure.txt')

    @property
    def gg85(self) -> str:
        return os.path.join(self.db, 'gg85.fasta')

    @property
    def ko30(self) -> str:
        return os.path.join(self.db, 'ko30.fasta')

    @property
    def ko30_structure(self) -> str:
        return os.path.join(self.db, 'ko30_structure.txt')

    @property
    def extracted(self) -> str:
        if self.indir is None:
            return os.path.join(self.outdir, 'extracted.fa')
        else:
            return os.path.join(self.indir, 'extracted.fa')

    @property
    def metadata(self) -> str:
        if self.indir is None:
            return os.path.join(self.outdir, 'metadata.txt')
        else:
            return os.path.join(self.indir, 'metadata.txt')

    @property
    def blastout(self) -> str:
        return os.path.join(self.outdir, 'blastout.txt')

    @property
    def blastout_filtered(self) -> str:
        return os.path.join(self.outdir, 'blastout.filtered.txt')

    @property
    def extracted_filtered(self) -> str:
        return os.path.join(self.outdir, 'extracted.filtered.fa')

    @property
    def columns(self) -> list:
        return ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']
