"""Alignment Modules used in president."""
import subprocess
import tempfile
import logging
import pandas as pd

logger = logging.getLogger(__name__)

PSL_LABELS = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
              'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
              'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
              'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
              'QStarts', 'TStarts']


def pblat(threads, reference, query):
    """Perform blat alignment.

    Parameters
    ----------
    threads : int
        Number of threads to use.
    reference : str
        file location of reference FASTA.
    query : str
        file location of query FASTA.

    Returns
    -------
    alignments.

    """
    _, alignments = tempfile.mkstemp(suffix="_pblat.tsv")
    cmd = f'pblat -threads={threads} {reference} {query} {alignments}'
    logger.info(f'Running pblat with: {cmd}')
    _ = subprocess.check_output(cmd, shell=True)
    logger.info('Finished pblat.')
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
    try:
        alignments = pd.read_csv(alignment_file,
                                 header=None, sep='\t', skiprows=5)
        alignments.columns = PSL_LABELS
        alignments.QName = alignments.QName.astype(str)

    except pd.errors.EmptyDataError:
        logger.error("Error reading pblat output. Perhaps it did not align anything (pandas).")
        return pd.DataFrame()

    except Exception:
        logger.error("Error reading pblat output. Perhaps it did not align anything (exception).")
        return pd.DataFrame()
    return alignments
