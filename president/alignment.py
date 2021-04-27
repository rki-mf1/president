"""Alignment Modules used in president."""
import subprocess
import tempfile
import logging
import pandas as pd
import os

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


def diamond(threads, reference, query):
    """Perform diamond alignment.

    Parameters
    ----------
    threads : int
        Number of threads to use.
    reference : str
        file location of reference FASTP.
    query : str
        file location of query FASTN.

    Returns
    -------
    alignments.

    """
    _, db = tempfile.mkstemp(suffix="_diamond_db")
    _, alignments = tempfile.mkstemp(suffix="_diamond.xml")
    makedb = f'diamond makedb --in {reference} -d {db}'
    logger.info(f'Running diamond makedb with: {makedb}')
    _ = subprocess.check_output(makedb, shell=True)
    blastx = f'diamond blastx -d {db} --threads {threads} --outfmt 5 --frameshift 15 -q {query} -o {alignments}'
    logger.info(f'Running diamond blastx with: {blastx}')
    _ = subprocess.check_output(blastx, shell=True)
    os.remove(db + ".dmnd")
    os.remove(db)
    logger.info('Finished diamond.')
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
