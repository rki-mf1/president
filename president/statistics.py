"""Module to compute alignment statistics."""
from collections import Counter
from datetime import datetime
import tempfile
import subprocess

import pandas as pd
import numpy as np
import screed
import re

from president import alignment, writer


def get_largest_N_gap(sequence):
    """
    Compute the largest gap size with consecutive Ns.

    Parameters:
        sequence: str, sequence
    Returns:
        int,
            size of largest NGAP

    Code from https://github.com/connor-lab/ncov2019-artic-nf/blob/master/bin/qc.py#L75.
    # n_pos = [i for i, letter in enumerate(sequence.lower()) if letter == 'n']
    # n_pos = [0] + n_pos + [len(sequence)]
    # n_gaps = [j - i for i, j in zip(n_pos[:-1], n_pos[1:])]
    # sorted(n_gaps)[-1]
    """
    # easier and correct, but maybe slow
    ngap_sizes = [len(i) for i in re.findall("N+", sequence)]
    if len(ngap_sizes) == 0:
        return 0
    else:
        return max(ngap_sizes)


def summarize_query(query):
    """
    Compute statistics for the sequences in the query.

    Parameters:
    ----------
    query: str,
        location of the FASTA query.
    """
    # get number of sequences and init results
    cmd = f'grep ">" {query} | wc -l'
    n_seqs = int(subprocess.check_output(cmd, shell=True))

    # init data table
    statistics_ar = writer.init_metrics(n_seqs)

    with screed.open(query) as seqfile:
        for idx, qry in enumerate(seqfile):
            seq = qry.sequence

            # sequence info
            statistics_ar["query_name"][idx] = qry.name
            statistics_ar["query_description"][idx] = qry.description
            statistics_ar["query_index"][idx] = idx

            # get nucleotide counts
            acgts_ct, iupacs_ct, nonupac_ct, Ns = count_nucleotides(seq)

            statistics_ar["acgt_bases"][idx] = acgts_ct
            statistics_ar["iupac_bases"][idx] = iupacs_ct
            statistics_ar["non_upac_bases"][idx] = nonupac_ct
            statistics_ar["N_bases"][idx] = Ns
            statistics_ar["length_query"][idx] = len(seq.rstrip())
            statistics_ar["Ngap"][idx] = get_largest_N_gap(seq)

    return pd.DataFrame(statistics_ar)


def qc_check(reference, query_stats, id_threshold=0.9):
    """
    Perform logic checks on already computed qc metrics (length, Ns).

    Parameters
    ----------
    reference : str or int
        Either the location of the reference or the length of the sequence.

    query_stats : df
        dataframe with already computed qc metrics for the sequences, see summarize_query.
    id_threshold : float, optional
        Identity threshold for assigning valid filters. The default is 0.9.

    Returns
    -------
    None.

    """
    if type(reference) == str:
        # read reference file
        with screed.open(reference) as seqfile:
            length_ref = [len(i.sequence) for i in seqfile][0]
    else:
        length_ref = reference

    # just store
    query_stats["length_reference"] = length_ref

    # valid characters?
    query_stats["qc_valid_nucleotides"] = query_stats["non_upac_bases"] == 0

    # compute if N/length threshold is reached
    # round, generously to avoid border line cases
    query_stats["qc_valid_number_n"] = \
        ((length_ref - query_stats["N_bases"]) / length_ref) >= id_threshold

    # take the lower bound here (floor)
    query_stats["qc_valid_length"] = query_stats["length_query"] >= \
        np.floor(length_ref * id_threshold)

    query_stats["qc_all_valid"] = query_stats["qc_valid_nucleotides"] & \
        query_stats["qc_valid_number_n"] & \
        query_stats["qc_valid_length"]


def nucleotide_identity(alignment_file, summary_stats_query, id_threshold=0.9):
    """Calculate nucleotide ident from a 2-sequence MSA.

    The query can consists of a multi-fasta file.

    Parameters
    ----------
    alignment_file : str
        alignment FASTA location.

    summary_stats_query : dataframe
        Summary statistics.

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
    # deal with no alignments
    if alignments.shape == (0, 0):
        return alignments

    # remove redundant alignments, since first entry is best match
    alignments = alignments.drop_duplicates(subset=["QName"], keep="first").reset_index()
    alignments.columns = "pblat_" + alignments.columns

    idmap = {i: j for i, j in
             zip(summary_stats_query["query_name"].str.replace("%space%", " "),
                 summary_stats_query["query_index"])}

    # retrieve index for joining data
    alignments["query_index"] = alignments["pblat_QName"].str.replace("%space%", " ").map(idmap)

    president_df = summary_stats_query.merge(alignments, left_on="query_index",
                                             right_on="query_index", how="left",
                                             suffixes=("", "_pblat"))

    president_df["identities"] = np.nan
    president_df["ambiguous_identities"] = np.nan
    president_df["iupac_ambiguous_identities"] = np.nan
    president_df["qc_post_aligned"] = False

    # get ids of values to overwrite
    idx = alignments["pblat_QName"].str.replace("%space%", " ").map(idmap).values

    president_df.at[idx, "qc_post_aligned"] = True
    # Metric valid_sequences #2 (A)
    # Ns in the query count as mismatch
    president_df.at[idx, "identities"] = \
        president_df["pblat_Matches"] / president_df[["pblat_TSize", "length_query"]].max(axis=1)

    # Ns in the query don't count
    # q=query, t=target sequence length
    president_df.at[idx, "ambiguous_identities"] = \
        president_df["pblat_Matches"] / \
        (president_df[["pblat_TSize", "length_query"]].max(axis=1) - president_df["N_bases"])

    # Ns in the query don't count
    # q=query, t=target sequence length
    president_df.at[idx, "iupac_ambiguous_identities"] = \
        president_df["pblat_Matches"] / \
        (president_df[["pblat_TSize", "length_query"]].max(axis=1)
         - president_df[["iupac_bases", "non_upac_bases"]].sum(axis=1))

    # passed id threshold
    president_df["qc_post_align_pass_threshold"] = president_df["identities"] >= id_threshold

    president_df["qc_post_aligned_all_valid"] = president_df["qc_post_align_pass_threshold"] &\
        president_df["qc_all_valid"]
    # sort by invalid sequences first
    president_df = president_df.sort_values(by="qc_post_aligned_all_valid")
    president_df = president_df.round(4)
    # add processing date
    president_df['Date'] = datetime.today().strftime('%Y-%m-%d')

    # cleanup dataframe and improve reporting format
    president_df = president_df.rename(
        columns={"pblat_Matches": "Matches",
                 "pblat_Mismatches": "Mismatches",
                 "pblat_TName": "reference_name",
                 "identities": "ACGT Nucleotide identity",
                 "ambiguous_identities": "ACGT Nucleotide identity (ignoring Ns)",
                 "iupac_ambiguous_identities": "ACGT Nucleotide identity (ignoring non-ACGTNs)"})

    president_df["query_name"] = president_df["query_name"].str.replace("%space%", " ")
    president_df["reference_name"] = president_df["reference_name"].str.replace("%space%", " ")

    president_df = president_df.drop(president_df.filter(regex="pblat").columns, axis=1)
    president_df = president_df[sorted(president_df.columns)]
    return president_df


def estimate_max_invalid(reference, id_threshold=0.9):
    """
    Estimate the lower bound (floor) of the sequence length and id_threshold in number of bases.

    Parameters
    ----------
    reference : str
        location of FASTA reference.
    id_threshold : float, optional
        Identity threshold for basepairs matching between query and reference. The default is 0.9.

    Returns
    -------
    int, maximal number of Ns in the sequence

    """
    # read reference file
    with screed.open(reference) as seqfile:
        length_ref = [len(i.sequence) for i in seqfile][0]
    return length_ref - np.ceil(id_threshold * length_ref)


def count_sequences(fasta_file, kind="reference"):
    """
    Test the number of reference sequences and check if # entries is equal to 1.

    Parameters
    ----------
    fasta_file : str
        location of FASTA reference.
    kind: str,
          Parameter indicating if reference or query is used as input. (default: reference)

    Returns:
    -------
    bool, (1 if the input is valid, 0 if the input is invalid)
    """
    # read reference file
    with screed.open(fasta_file) as seqfile:
        nseqs = int(np.sum([1 for i in seqfile]))

    # valid, only 1 file
    if nseqs == 1:
        print("Number of references is equal to one (qc passed)")
        return 1

    elif (nseqs == 0) | (nseqs > 1):
        # if reference, we exactly need one file
        if kind == "reference":
            raise ValueError(f"Number of reference sequences ({nseqs}) is not equal to 1.")
        else:
            # if query, only 0 seqs is bad
            if nseqs == 0:
                return 0
            else:
                return 1


def split_valid_sequences(query, query_stats):
    """
    Perform simple qc checks on the query sequence and discard if input does not met qc.

    Based on the length of the reference sequence, a query sequence with with (l1 - n) / l1,
    will never reach the id_threshold (l1 = reference length, n= #N in query).

    Parameters
    ----------
    query : str
        location of FASTA query.
    query_stats : df
        Dataframe with qc checks of the query sequences.

    Returns
    -------
    triple,
        query location, diagnostic, identifier failed

    """
    all_valid = query_stats["qc_all_valid"]
    if all_valid.sum() == query_stats.shape[0]:
        # all good, all sequences pass qc
        return query, "all_valid", []

    elif all_valid.sum() == 0:
        # all bad, none of the sequences pass the qc
        invalid_loc = tempfile.mkstemp(suffix="_valid.fasta")[1]
        return invalid_loc, "all_invalid", \
            query_stats[~query_stats["qc_all_valid"]]["query_name"].values

    else:
        # some sequences pass qc, others not
        valid_sequences_idx = set(query_stats[query_stats["qc_all_valid"]]["query_index"].values)
        invalid_identifier = query_stats[~query_stats["qc_all_valid"]]["query_name"].values

        # write valid sequences
        valid_loc = tempfile.mkstemp(suffix="_valid.fasta")[1]
        valid_fout = open(valid_loc, "w")

        # write invalid sequences
        invalid_loc = tempfile.mkstemp(suffix="_invalid.fasta")[1]
        invalid_fout = open(invalid_loc, "w")

        # iterate over sequence again and split into valid / invalid
        with screed.open(query) as seqfile:
            for idx, qry in enumerate(seqfile):
                # if its a valid sequence, write to valid ...
                if idx in valid_sequences_idx:
                    seqwriter = valid_fout
                # else, we have an invalid sequence
                else:
                    seqwriter = invalid_fout
                writer.write_fasta(seqwriter, qry, format_sequence=True)
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
        acgt counts, iupac counts, non-iupac counts, Ns
    """
    # init count dict
    iupac_dic = init_UPAC_dictionary(acgt=True)
    # count nucleotides
    iupac_dic.update(Counter(sequence))

    # ref characters
    valid_nucleotides = init_UPAC_dictionary(acgt=False)

    acgt_counts = sum([iupac_dic[nt] for nt in "ACGT"])
    iupac_counts = sum([iupac_dic[nt] for nt in iupac_dic.keys() if nt in valid_nucleotides])
    return acgt_counts, iupac_counts, len(sequence) - acgt_counts - iupac_counts, iupac_dic["N"]
