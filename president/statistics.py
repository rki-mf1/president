#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import screed
import pandas as pd


def nucleotide_identity(query, alignments, max_invalid):
    """
    Calculate nucleotide ident from a 2-sequence MSA.

    Parameters
    ----------
    query : str
        query FASTA location.
    alignments : str
        alignment FASTA location.
    max_invalid : int
        maximal invalid entries in alignment.

    Raises
    ------
    ValueError
        Raised when max_invalid criteria is not fullfilled.

    Returns
    -------
    tuple
        ident, ident_non_canonical, non_canonical, len(qry)

    """

    with screed.open(query) as file:
        qry = next(file).sequence


    # Consider only sites where the query has non-ACTG characters
    # Metric issue #2 (B)
    non_canonical = sum([1 for i in qry if i not in 'ACTG'])
    if non_canonical > max_invalid:
        raise ValueError('Too many non-canonical nucleotides, abort!')


    # Read the alignment(s) with pandas
    # Pandas can be replaced with split to reduce dependencies
    alignments = pd.read_csv(
        alignments, header=None, sep='\t', skiprows=5)
    labels = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
              'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
              'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
              'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
              'QStarts', 'TStarts']
    alignments.columns = labels


    # Metric issue #2 (A)
    # Ns in the query count as mismatch
    ident = round(alignments.at[0, 'Matches'] / len(qry), 4)


    # Ns in the query don't count
    ident_non_canonical = round(alignments.at[0, 'Matches'] / 
                                (len(qry) - non_canonical), 4)


    return ident, ident_non_canonical, non_canonical, len(qry)