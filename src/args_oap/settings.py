import logging
import os

## setup logger
logging.basicConfig(
    level="INFO",
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")

## bold format
BOLD = "\033[1m"
RESET = "\033[0m"

## add color
logging.addLevelName( logging.WARNING, BOLD+"\x1b[33;20m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName( logging.CRITICAL, BOLD+"\x1b[31;20m%s\033[1;0m" % logging.getLevelName(logging.CRITICAL))

logger = logging.getLogger(__name__)

## setup default filepath
class Settings:
    def __init__(self):
        self._path = os.path.dirname(os.path.realpath(__file__))
        self._sarg = os.path.join(self._path, 'db', 'sarg.fasta')
        self._sarg_structure1 = os.path.join(self._path, 'db', 'single-component_structure.txt')
        self._sarg_structure2 = os.path.join(self._path, 'db', 'two-component_structure.txt')
        self._sarg_structure3 = os.path.join(self._path, 'db', 'multi-component_structure.txt')
        self._gg85 = os.path.join(self._path, 'db', 'gg85.fasta')
        self._ko30 = os.path.join(self._path, 'db', 'all_KO30.fasta')
        self._ko30_structure = os.path.join(self._path, 'db', 'all_KO30_name.list')
        self.cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore'] # blast6 output format

settings = Settings()