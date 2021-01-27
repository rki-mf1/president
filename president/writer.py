"""Module for writing and organizing results."""
import os

import screed

from president import sequence


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


def write_fasta(fileobj, seq, format_sequence=False):
    """
    Write fasta sequence to file.

    Parameters
    ----------
    fileobj: obj,
        file object (opened)
    seq: seq object,
        screed sequence object

    Returns
    -------
        None
    """
    fileobj.write(f">{seq.name} {seq.description}\n")
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
        valid_ids = set(metrics[metrics["Valid"] & metrics["aligned"]]["ID"].values)
    else:
        valid_ids = set()

    valid_name = os.path.join(os.path.dirname(out_dir), f"{prefix}valid.fasta")
    valid_fout = open(valid_name, write_mode)
    invalid_name = os.path.join(os.path.dirname(out_dir), f"{prefix}invalid.fasta")
    invalid_fout = open(invalid_name, write_mode)

    # iterate over query and split into valid / invalid fasta
    with screed.open(query) as seqfile:
        for seq in seqfile:
            if seq.name in valid_ids:
                write_fasta(valid_fout, seq, format_sequence=True)
            else:
                write_fasta(invalid_fout, seq)

    print(f"Valid sequences written to: {valid_name}")
    print(f"Invalid sequences written to: {invalid_name}")
    valid_fout.close()
    invalid_fout.close()
