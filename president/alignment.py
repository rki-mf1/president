"""Alignment Modules used in president."""
import subprocess
import tempfile

import pandas as pd

PSL_LABELS = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
              'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
              'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
              'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
              'QStarts', 'TStarts']


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
    print('Running pblat ...')
    _, alignments = tempfile.mkstemp(suffix="_pblat.tsv")
    cmd = f'pblat -threads={threads} {reference} {query} {alignments}'
    if verbose:
        print(cmd)
    _ = subprocess.check_output(cmd, shell=True)
    print('Finished pblat.')
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

    except pd.errors.EmptyDataError:
        print("Error reading pblat output. Perhaps it did not align anything (pandas).")
        return pd.DataFrame()

    except Exception:
        print("Error reading pblat output. Perhaps it did not align anything (exception).")
        return pd.DataFrame()
    return alignments
