"""Module to compute alignment statistics."""
from collections import Counter
from datetime import datetime
import tempfile

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
    acgt_bases = np.zeros(n_seqs)
    no_iupac_bases = np.zeros(n_seqs)

    identities = np.zeros(n_seqs)
    ambiguous_identities = np.zeros(n_seqs)
    iupac_ambiguous_identities = np.zeros(n_seqs)

    query_lengths = np.zeros(n_seqs, dtype=np.uint32)

    # collect failed ids here
    failed_ids = []
    with screed.open(query) as seqfile:
        idx = 0
        for qry in seqfile:
            # nucleotide counts
            acgts_ct, iupacs_ct, nonupac_ct = count_nucleotides(qry.sequence)

            # increase idx if pblat print multiple alignments
            # since the first one is the best match anyway
            while idx > 0 and idx < alignments.shape[0] and \
                    alignments.at[idx, 'QName'] == alignments.at[idx - 1, 'QName']:
                query_ids[idx] = np.nan
                idx = idx + 1

            if qry.name.startswith(alignments.at[idx, 'QName']):
                # basic sequence info
                query_ids[idx] = qry.name.replace("%space%", " ")

                query_lengths[idx] = len(qry.sequence)
                acgt_bases[idx] = acgts_ct
                no_iupac_bases[idx] = iupacs_ct

                # Consider only sites where the query has non-ACTG characters
                # Metric issue #2 (B)
                ambiguous_bases[idx] = sum([1 for i in qry.sequence if i not in 'ACTG'])

                # Metric valid_sequences #2 (A)
                # Ns in the query count as mismatch
                identities[idx] = \
                    alignments.at[idx, 'Matches'] / \
                    max(len(qry.sequence),
                        alignments.at[idx, 'TSize'])

                # Ns in the query don't count
                # q=query, t=target sequence length
                ambiguous_identities[idx] = \
                    alignments.at[idx, 'Matches'] / \
                    (max(len(qry.sequence),
                         alignments.at[idx, 'TSize']) - ambiguous_bases[idx])

                # Ns in the query don't count
                # q=query, t=target sequence length
                iupac_ambiguous_identities[idx] = \
                    alignments.at[idx, 'Matches'] / \
                    (max(len(qry.sequence),
                         alignments.at[idx, 'TSize']) - ambiguous_bases[idx] - no_iupac_bases[idx])
                idx = idx + 1

            else:
                failed_ids.append(qry.name.replace("%space%", " "))
                print(f"The sequence: '{failed_ids[-1]}' could not be aligned with pblat.")

    # format to single result dataframe
    metrics = pd.DataFrame({
        'ID': query_ids,
        'Valid': identities >= id_threshold,
        'ACGT Nucleotide identity': identities,
        'ACGT Nucleotide identity (ignoring Ns)': ambiguous_identities,
        'ACGT Nucleotide identity (ignoring non-ACGTNs)': iupac_ambiguous_identities,
        'Ambiguous Bases': ambiguous_bases,
        'Query Length': query_lengths,
        'Query #ACGT': acgt_bases,
        'Query #IUPAC-ACGT': no_iupac_bases,
        'Query #non-IUPAC': query_lengths - no_iupac_bases - acgt_bases,
        'aligned': True,
        'passed_initial_qc': True
    })

    # add not-aligned sequences to result
    if len(failed_ids) > 0:
        invalid_df = pd.DataFrame({"ID": failed_ids})
        invalid_df["passed_initial_qc"] = True
        invalid_df["aligned"] = False
        metrics = pd.concat([metrics, invalid_df]).reset_index(drop=True)

    # sort by invalid sequences first
    metrics = metrics.sort_values(by="Valid")
    # drop 2nd, 3rd best alignment rows for the input sequences
    metrics = metrics.dropna(subset=['ID'])
    metrics = metrics.round(4)
    # add processing date
    metrics['Date'] = [datetime.today().strftime('%Y-%m-%d')] * metrics.shape[0]
    return metrics


def metrics_all_invalid(query, n_seqs):
    """Calculate query metrics for case of invalid only sequences.

    Parameters
    ----------
        query : str
            query FASTA location.

        n_seqs : int
            number of sequences in query

    Returns
    -------
    pandas Dataframe
        ID, Valid, ACGT Nucleotide identity, ACGT Nucleotide identity (ignoring Ns),
        ACGT Nucleotide identity (ignoring non-ACGTNs), Ambiguous Bases, Query Length,
        Query #ACGT, Query #IUPAC-ACGT, Query #non-IUPAC, aligned, passed_initial_qc
    """

    query_ids = np.empty(n_seqs, dtype="object")

    ambiguous_bases = np.zeros(n_seqs)
    acgt_bases = np.zeros(n_seqs)
    no_iupac_bases = np.zeros(n_seqs)

    ambiguous_identities = np.zeros(n_seqs)
    iupac_ambiguous_identities = np.zeros(n_seqs)

    query_lengths = np.zeros(n_seqs, dtype=np.uint32)

    with screed.open(query) as seqfile:
        idx = 0
        for qry in seqfile:
            # nucleotide counts
            acgts_ct, iupacs_ct, nonupac_ct = count_nucleotides(qry.sequence)
            
            # basic sequence info
            query_ids[idx] = qry.name.replace("%space%", " ")

            query_lengths[idx] = len(qry.sequence)
            acgt_bases[idx] = acgts_ct
            no_iupac_bases[idx] = iupacs_ct

            # Consider only sites where the query has non-ACTG characters
            # Metric issue #2 (B)
            ambiguous_bases[idx] = Counter(qry.sequence)["N"]

            idx = idx + 1

    # format to single result dataframe
    metrics = pd.DataFrame({
        'ID': query_ids,
        'Valid': False,
        'ACGT Nucleotide identity': "",
        'ACGT Nucleotide identity (ignoring Ns)': "",
        'ACGT Nucleotide identity (ignoring non-ACGTNs)': "",
        'Ambiguous Bases': ambiguous_bases,
        'Query Length': query_lengths,
        'Query #ACGT': acgt_bases,
        'Query #IUPAC-ACGT': no_iupac_bases,
        'Query #non-IUPAC': query_lengths - no_iupac_bases - acgt_bases,
        'aligned': False,
        'passed_initial_qc': False
    })

    # add processing date
    metrics['Date'] = [datetime.today().strftime('%Y-%m-%d')] * metrics.shape[0]
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

    # get number of invalid sequences, and create indicator
    valid_nucleotides = np.array([count_nucleotides(seq.sequence)[2] for seq in queries])
    valid_nucleotides = (valid_nucleotides == 0)

    # compute if N/length threshold is reached
    # round, generously to avoid border line cases
    valid_Ns = ((length_ref - Ns_queries) / length_ref) >= id_threshold

    # take the lower bound here (floor)
    valid_lengths = length_queries >= np.floor(length_ref * id_threshold)

    # concat conditions
    # only if all conditions true, align
    all_valid = valid_Ns & valid_lengths & valid_nucleotides

    if all_valid.sum() == len(Ns_queries):
        # all good, all sequences pass qc
        return query, "all_valid", []

    elif all_valid.sum() == 0:
        # all bad, none of the sequences pass the qc
        invalid_sequences = [queries[idx] for idx, cond in enumerate(all_valid) if not cond]
        invalid_loc = tempfile.mkstemp(suffix="_invalid.fasta")[1]
        # write invalid sequences
        with open(invalid_loc, "w") as fout:
            for vseq in invalid_sequences:
                fout.write(f">{vseq['name']} {vseq['description']}\n")
                fout.write(f"{vseq['sequence']}\n")
        return invalid_loc, "all_invalid", \
            [qry.name.replace("%space%", " ") for qry in invalid_sequences]

    else:
        # some sequences pass qc, others not
        valid_sequences = [queries[idx] for idx, cond in enumerate(all_valid) if cond]
        invalid_sequences = [queries[idx] for idx, cond in enumerate(all_valid) if not cond]
        invalid_identifier = [
            queries[idx].name.replace("%space%", " ")
            for idx, cond in enumerate(all_valid) if not cond]

        # write valid sequences only
        valid_loc = tempfile.mkstemp(suffix="_valid.fasta")[1]
        with open(valid_loc, "w") as fout:
            for vseq in valid_sequences:
                fout.write(f">{vseq['name']} {vseq['description']}\n")
                fout.write(f"{vseq['sequence']}\n")

        invalid_loc = tempfile.mkstemp(suffix="_invalid.fasta")[1]
        # write invalid sequences only
        with open(invalid_loc, "w") as fout:
            for vseq in invalid_sequences:
                fout.write(f">{vseq['name']} {vseq['description']}\n")
                fout.write(f"{vseq['sequence']}\n")

        return valid_loc, "mixed", invalid_identifier


def init_UPAC_dictionary(acgt=False):
    """
    Generate a count template dictionary for IUPAC characters excluding ACGT.

    Definitions taken from here: https://www.bioinformatics.org/sms/iupac.html.
    However, "-", "."

    Parameters
    ----------
    acgt : bool, optional
        If True include ACGT in the dictionary. The default is False.

    Returns
    -------
    dict
        IUPAC (char:0) dictionary.

    """
    if acgt:
        return {i: 0 for i in "ACGTURYSWKMBDHVN."}
    else:
        return {i: 0 for i in "URYSWKMBDHVN."}


def count_nucleotides(sequence):
    """
    Compute IUPAC nucleotide statistics and return counts.

    Parameters
    ----------
    sequence : str
        sequence string.

    Returns
    -------
    triple,
        acgt counts, iupac counts, non-iupac counts
    """
    # init count dict
    iupac_dic = init_UPAC_dictionary(acgt=True)
    # count nucleotides
    iupac_dic.update(Counter(sequence))

    # ref characters
    valid_nucleotides = init_UPAC_dictionary(acgt=False)

    acgt_counts = sum([iupac_dic[nt] for nt in "ACGT"])
    iupac_counts = sum([iupac_dic[nt] for nt in iupac_dic.keys() if nt in valid_nucleotides])
    return acgt_counts, iupac_counts, len(sequence) - acgt_counts - iupac_counts
