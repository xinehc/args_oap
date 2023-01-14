import logging
import os

## setup default filepath
class Settings:
    def __init__(self):
        self._path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../resources/args_oap'))
        self._sarg = os.path.join(self._path, 'db', 'sarg.fasta')
        self._sarg_structure1 = os.path.join(self._path, 'db', 'single-component_structure.txt')
        self._sarg_structure2 = os.path.join(self._path, 'db', 'two-component_structure.txt')
        self._sarg_structure3 = os.path.join(self._path, 'db', 'multi-component_structure.txt')
        self._gg85 = os.path.join(self._path, 'db', 'gg85.fasta')
        self._ko30 = os.path.join(self._path, 'db', 'ko30.fasta')
        self._ko30_structure = os.path.join(self._path, 'db', 'ko30_structure.txt')
        self.cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore'] # blast6 output format
        print("args_oap res dir:",self._path)

settings = Settings()