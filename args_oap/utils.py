import os
import logging


def count_bp(file, skipchar):
    """
    Count base pairs in a fa/fp file
    """
    nbp = 0
    nseq = 0
    with open(file, "r") as f:
        for line in f:
            if line[0]!=skipchar:
                nbp += len(line.rstrip())
                nseq += 1
    return nbp, nseq

## set up the logger
logging.basicConfig(
    level="INFO",
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")

## set up bold format
BOLD = "\033[1m"
RESET = "\033[0m"

## add color
logging.addLevelName( logging.WARNING, BOLD+"\x1b[33;20m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName( logging.CRITICAL, BOLD+"\x1b[31;20m%s\033[1;0m" % logging.getLevelName(logging.CRITICAL))

logger = logging.getLogger(__name__)
