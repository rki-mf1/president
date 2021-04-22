"""Module for writing and organizing results."""
import logging
import os
from datetime import datetime

import numpy as np
import screed

from president import alignment
from president import sequence

logger = logging.getLogger(__name__)

PSL_LABELS = alignment.PSL_LABELS

COL_ORDER = ["query_name", "reference_name", "file_in_query", "file_in_ref",
             "ACGT Nucleotide identity", "ACGT Nucleotide identity (ignoring Ns)",
             "ACGT Nucleotide identity (ignoring non-ACGTNs)", "qc_all_valid",
             "qc_is_empty_query", "qc_post_align_pass_threshold", "qc_post_aligned",
             "qc_post_aligned_all_valid", "qc_valid_length", "qc_valid_nucleotides",
             "qc_valid_pass_nthreshold", "acgt_bases", "iupac_bases",
             "non_iupac_bases", "N_bases", "length_query", "length_reference",
             "LongestNGap", "Matches", "Mismatches", "query_description",
             "query_index", "Date"]


def init_metrics(n_seqs, extend_cols=False, metrics_df=None, store_alignment=False):
    """
    Init output data structure.

    Parameters:
    ----------
    n_seqs: int,
            number of sequences to process / store meta values
    extend_cols: bool,
                if True, extends the columns to presidents output format.
                Needed for consistent writing in cases where no sequences were aligned.
    metrics_df:
        Potential Dataframe to be extended (default: None). Must be used with extend_cols.

    Returns:
    -------
    np structured array / pandas dataframe
    """
    statistic_values = np.dtype([("query_description", "object"),
                                 ("query_name", "object"),
                                 ("query_index", "object"),
                                 ("acgt_bases", np.uint32),
                                 ("iupac_bases", np.uint32),
                                 ("non_iupac_bases", np.uint32),
                                 ("N_bases", np.uint32),
                                 ("length_query", np.uint32),
                                 ("LongestNGap", np.uint32)])
    stats_ar = np.zeros(n_seqs, dtype=statistic_values)

    if extend_cols:
        # extend output data to be equal to the standard output
        needed_cols = \
            ["ACGT Nucleotide identity", "ACGT Nucleotide identity (ignoring Ns)",
             "ACGT Nucleotide identity (ignoring non-ACGTNs)", "Matches", "Mismatches",
             "qc_post_align_pass_threshold", "qc_post_aligned",
             "qc_post_aligned_all_valid", "reference_name", "length_reference"]
        if store_alignment:
            needed_cols += \
                ["PSL_"+c for c in PSL_LABELS]

        for coli in needed_cols:
            if "qc_post" in coli:
                if metrics_df.shape[0] == 0:
                    metrics_df[coli] = [False]
                else:
                    metrics_df[coli] = metrics_df.shape[0] * [False]
            else:
                metrics_df[coli] = np.nan
        metrics_df['Date'] = datetime.today().strftime('%Y-%m-%d')
        return metrics_df

    else:
        return stats_ar


def get_filename(file):
    """
    Extract filename of the file.

    Parameters
    ----------
    file : str
        file location.
    Returns
    -------
    str, filename without extension
    """
    return os.path.splitext(os.path.basename(file))[0]


def write_fasta(fileobj, seq, format_sequence=False, format_pblat=False):
    """
    Write fasta sequence to file.

    Parameters
    ----------
    fileobj: obj,
        file object (opened)
    seq: seq object,
        screed sequence object
    format_sequence: bool,
        If True, transform sequence to valid upper case sequence (default: False)
    Returns
    -------
        None
    """
    # pblat cuts away the description, fiddle around this issue by having the " " decoded as %space%
    if not format_pblat:
        sequence_header = f">{seq.name} {seq.description}".replace("%space%", " ")
    else:
        sequence_header = f">{seq.name}{seq.description}"
    fileobj.write(sequence_header.rstrip()+"\n")
    if format_sequence:
        # make sure to only store upper case, ACGT symbols. replace all others with "N"s
        fileobj.write(sequence.to_valid_upper(seq.sequence) + "\n")
    else:
        fileobj.write(f"{seq.sequence}\n")


def write_sequences(query, metrics, out_dir, evaluation, write_mode="w"):
    """
    Write sequences from the query into valid / invalid FASTA files based on metrics.

    Parameters
    ----------
    query: str,
        file location of the query FASTA
    metrics: dataframe,
            pandas dataframe filled with valid and aligned columnns
    out_dir: str,
            Output directory to store FASTAs

    Returns
    -------
        None
    """
    prefix = os.path.basename(out_dir)
    if evaluation != "all_invalid":
        # valid ids were aligned and pass the qc
        valid_ids = set(metrics[metrics["qc_post_aligned_all_valid"]]["query_name"].values)
    else:
        valid_ids = set()

    # split valid and invalid file writers
    valid_name = os.path.join(os.path.dirname(out_dir), f"{prefix}valid.fasta")
    valid_fout = open(valid_name, write_mode)
    invalid_name = os.path.join(os.path.dirname(out_dir), f"{prefix}invalid.fasta")
    invalid_fout = open(invalid_name, write_mode)

    # iterate over query and split into valid / invalid fasta
    nvalid = 0
    ninvalid = 0
    with screed.open(query) as seqfile:
        for seq in seqfile:
            if seq.name.replace("%space%", " ") in valid_ids:
                nvalid += 1
                write_fasta(valid_fout, seq, format_sequence=True)
            else:
                write_fasta(invalid_fout, seq)
                ninvalid += 1

    logger.info(f"Valid (n={nvalid}) sequences written to: {valid_name}")
    logger.info(f"Invalid (n={ninvalid}) sequences written to: {invalid_name}")
    valid_fout.close()
    invalid_fout.close()
