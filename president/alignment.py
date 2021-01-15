"""Alignment Modules used in president."""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tempfile
import subprocess

def pblat(threads, reference, query, verbose=0):
    """
    Wrapper to perform pblat alignment.

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
