"""Alignment Modules used in president."""
import tempfile
import subprocess

import pandas as pd


def remove_spaces(fasta):
    """Remove white spaces from IDs in a FASTA file.

    Parameters
    ----------
    fasta : str
        FASTA file to remove white spaces.

    Returns
    -------
    processed FASTA file.

    """
    _, output = tempfile.mkstemp()
    with open(output, "w") as fout:
        with open(fasta, "r") as fin:
            for line in fin:
                line = line.strip()
                if line[0] == '>':
                    fout.write(line.replace(" ", "%space%") + "\n")
                else:
                    fout.write(line + "\n")
    return output


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
    _, alignments = tempfile.mkstemp()
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
    alignments = pd.read_csv(alignment_file, header=None, sep='\t', skiprows=5)
    labels = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
              'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
              'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
              'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
              'QStarts', 'TStarts']
    alignments.columns = labels
    return alignments
