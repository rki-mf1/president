"""Alignment Modules used in president."""
import logging
import subprocess
import tempfile
import time

import pandas as pd

logger = logging.getLogger(__name__)


def pblat(threads, reference, query, verbose=0):
    """Perform blat alignment.

    Parameters
    ----------
    threads : int
        Number of threads to use.
    reference : str
        file location of reference FASTA.
    query : str
        file location of query FASTA.
    verbose : bool
        If True, print pblat command.

    Returns
    -------
    alignments.

    """
    _, alignments = tempfile.mkstemp()
    cmd = f'pblat -threads={threads} {reference} {query} {alignments}'
    logger.info(f'Running pblat with: {cmd}')
    start = time.time()
    _ = subprocess.check_output(cmd, shell=True)
    end = time.time()
    logger.info(f'Finished pblat. (took: {(end - start) / 60.:.6f} minutes)')
    return alignments


def parse_alignment(alignment_file):
    """
    Parse a given alignment file into a pandas dataframe.

    Parameters
    ----------
    alignment_file : str
        location of the alignment file from pblat.

    Returns
    -------
    None.

    """
    # Read the alignment(s) with pandas
    # Pandas can be replaced with split to reduce dependencies
    logger.info(f"Parsing alignments from file: {alignment_file}")
    try:
        alignments = pd.read_csv(alignment_file, header=None, sep='\t', skiprows=5)
        labels = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
                  'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
                  'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
                  'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
                  'QStarts', 'TStarts']
        alignments.columns = labels
        logger.info("Successfully parsed alignment.")
    except Exception:
        logger.error("Error reading pblat output. Perhaps it did not align anything.")
        exit(1)
    return alignments
