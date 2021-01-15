"""Module to compute alignment statistics."""
import pandas as pd
import numpy as np
import screed
from president import alignment


def nucleotide_identity(query, alignment_file, max_invalid):
    """Calculate nucleotide ident from a 2-sequence MSA.

    The query can consists of a multi-fasta file.

    Parameters
    ----------
    query : str
        query FASTA location.

    alignment_file : str
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
    alignments = alignment.parse_alignment(alignment_file)
    # get number of sequences and init results
    n_seqs = alignments.shape[0]

    query_ids = np.empty(n_seqs, dtype="object")
    ambiguous_bases = np.zeros(n_seqs)
    identities = np.zeros(n_seqs)
    ambiguous_identities = np.zeros(n_seqs)
    query_lengths = np.zeros(n_seqs, dtype=np.uint32)

    with screed.open(query) as seqfile:
        for idx, qry in enumerate(seqfile):
            # basic sequence info
            query_ids[idx] = qry.name

            # Consider only sites where the query has non-ACTG characters
            # Metric issue #2 (B)
            ambiguous_bases[idx] = sum([1 for i in qry.sequence if i not in 'ACTG'])

            # Metric valid_sequences #2 (A)
            # Ns in the query count as mismatch
            identities[idx] = alignments.at[idx, 'Matches'] / len(qry.sequence)

            # Ns in the query don't count
            ambiguous_identities[idx] = \
                alignments.at[idx, 'Matches'] / (len(qry.sequence) - ambiguous_bases[idx])

            query_lengths[idx] = len(qry.sequence)

    # format to single result dataframe
    metrics = pd.DataFrame({
        'ID': query_ids,
        'Valid': ambiguous_bases > max_invalid,
        'Identity': identities,
        'Ambiguous Identity': ambiguous_identities,
        'Ambiguous Bases': ambiguous_bases,
        'Query Length': query_lengths
    })
    # sort by invalid sequences first
    metrics = metrics.sort_values(by="Valid")
    metrics = metrics.round(4)
    return metrics
