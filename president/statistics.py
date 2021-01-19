"""Module to compute alignment statistics."""
import pandas as pd
import numpy as np
import screed

from president import alignment


def nucleotide_identity(query, alignment_file, id_threshold=0.93):
    """Calculate nucleotide ident from a 2-sequence MSA.

    The query can consists of a multi-fasta file.

    Parameters
    ----------
    query : str
        query FASTA location.

    alignment_file : str
        alignment FASTA location.

    id_threshold : float
       minimal id threshold to pass

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
        idx = 0
        for qry in seqfile:
            # increase idx if pblat print multiple alignments
            # since the first one is the best match anyway
            while idx > 0 and idx < alignments.shape[0] and \
                    alignments.at[idx, 'QName'] == alignments.at[idx - 1, 'QName']:
                query_ids[idx] = np.nan
                idx = idx + 1
            if qry.name.startswith(alignments.at[idx, 'QName']):
                # basic sequence info
                query_ids[idx] = qry.name.replace("%space%", " ")

                # Consider only sites where the query has non-ACTG characters
                # Metric issue #2 (B)
                ambiguous_bases[idx] = sum([1 for i in qry.sequence if i not in 'ACTG'])

                # Metric valid_sequences #2 (A)
                # Ns in the query count as mismatch
                identities[idx] = \
                    alignments.at[idx, 'Matches'] / \
                    max(alignments.at[idx, 'QSize'],
                        alignments.at[idx, 'TSize'])

                # Ns in the query don't count
                ambiguous_identities[idx] = \
                    alignments.at[idx, 'Matches'] / \
                    (max(alignments.at[idx, 'QSize'],
                         alignments.at[idx, 'TSize']) - ambiguous_bases[idx])

                query_lengths[idx] = len(qry.sequence)

                idx = idx + 1

            else:
                print(qry.name.replace("%space%", " "), 'could not be aligned with pblat')

    # format to single result dataframe
    metrics = pd.DataFrame({
        'ID': query_ids,
        'Valid': identities >= id_threshold,
        'Identity': identities,
        'Ambiguous Identity': ambiguous_identities,
        'Ambiguous Bases': ambiguous_bases,
        'Query Length': query_lengths
    })
    # sort by invalid sequences first
    metrics = metrics.sort_values(by="Valid")
    # drop 2nd, 3rd best alignment rows for the input sequences
    metrics = metrics.dropna(subset=['ID'])
    metrics = metrics.round(4)
    return metrics


def estimate_max_invalid(reference, id_threshold=0.93):
    """
    Estimate the lower bound (floor) of the sequence length and id_threshold in number of bases.

    Parameters
    ----------
    reference : str
        location of FASTA reference.
    id_threshold : float, optional
        Identity threshold for basepairs matching between query and reference. The default is 0.93.

    Returns
    -------
    int, maximal number of Ns in the sequence

    """
    # read reference file
    with screed.open(reference) as seqfile:
        length_ref = [len(i.sequence) for i in seqfile][0]
    return length_ref - np.ceil(id_threshold * length_ref)


def count_reference_sequences(reference):
    """
    Test the number of reference sequences and check if # entries is equal to 1.

    Parameters
    ----------
    reference : str
        location of FASTA reference.

    """
    # read reference file
    with screed.open(reference) as seqfile:
        nrefs = np.sum([1 for i in seqfile])

    if nrefs == 1:
        print("Number of references is equal to one (qc passed)")

    elif (nrefs == 0) | (nrefs) > 1:
        raise ValueError(f"Number of reference sequences ({nrefs}) is not equal to 1.")


def split_valid_sequences(query, reference, id_threshold=0.93):
    """
    Perform simple qc checks on the query sequence and discard if input does not met qc.

    Based on the length of the reference sequence, a query sequence with with (l1 - n) / l1,
    will never reach the id_threshold (l1 = reference length, n= #N in query).

    Parameters
    ----------
    query : str
        location of FASTA query.
    reference : str
        location of FASTA reference.
    id_threshold : float, optional
        Identity threshold for basepairs matching between query and reference. The default is 0.93.

    Returns
    -------
    triple,
        query location, diagnostic, identifier failed

    """
    # read reference file
    with screed.open(reference) as seqfile:
        length_ref = [len(i.sequence) for i in seqfile][0]

    # read query file
    with screed.open(query) as seqfile:
        queries = [seq for seq in seqfile]

    # get number of Ns per query / seqlength
    Ns_queries = np.array([i.sequence.count("N") for i in queries])
    length_queries = np.array([len(i.sequence) for i in queries])

    # compute if N/length threshold is reached
    # round, generously to avoid border line cases
    valid_Ns = ((length_ref - Ns_queries) / length_ref) >= id_threshold

    # take the lower bound here (floor)
    valid_lengths = length_queries >= np.floor(length_ref * id_threshold)

    # concat conditions
    all_valid = valid_Ns & valid_lengths

    if all_valid.sum() == len(Ns_queries):
        # all good, all sequences pass qc
        return query, "all_valid", []

    elif all_valid.sum() == 0:
        # all good, all sequences pass qc
        return query+"_invalid.fasta", "all_invalid", \
            [queries[idx].name for idx, cond in enumerate(all_valid) if not cond]

    else:
        valid_sequences = [queries[idx] for idx, cond in enumerate(all_valid) if cond]
        invalid_sequences = [queries[idx] for idx, cond in enumerate(all_valid) if not cond]
        invalid_identifier = [
            queries[idx].name.replace("%space%", " ")
            for idx, cond in enumerate(all_valid) if not cond]

        # write valid sequences only
        valid_loc = query+"_valid.fasta"
        with open(valid_loc, "w") as fout:
            for vseq in valid_sequences:
                fout.write(f">{vseq['name']} {vseq['description']}\n")
                fout.write(f"{vseq['sequence']}\n")

        # write valid sequences only
        with open(query+"_invalid.fasta", "w") as fout:
            for vseq in invalid_sequences:
                fout.write(f">{vseq['name']} {vseq['description']}\n")
                fout.write(f"{vseq['sequence']}\n")

        return valid_loc, "mixed", invalid_identifier
