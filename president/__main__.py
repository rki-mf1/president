"""President's main module.

Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage:
-----
president --query tiny_test.masked_consensus.fasta \
    --reference NC_045512.2.fasta -x 0.93 -p 8 -o report.tsv

Notes:
-----
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line
(below) to ignore Ns

Requirements:

- conda install -y -c bioconda python=3.8 pblat screed pandas

ANI definition:

- https://pubmed.ncbi.nlm.nih.gov/17220447/

Spec:
-----
> Sequence identity to a reference genome sequence (e.g. NC_045512.2) should
be calculated at the nucleotide level relative to the entire length of the
reference. Only informative nucleotides (A,T,G,C) are considered identical
to each other.
"""
import argparse
import os
# authors:
# RKI MF1;  Martin Hoelzer with great initial help of @phiweger (UKL Leipzig)
# HPI;      Fabio Malcher Miranda, Sven Giese, Alice Wittig
import sys
from shutil import which

import pandas as pd
import screed

from president import alignment, __version__, statistics, writer, sequence


def is_available(name="pblat"):
    """
    Check whether `name` is on PATH and marked as executable.

    Parameters
    ----------
    name : str
        tool name executable (default: pblat).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # tool installed
    if not which(name):
        raise ValueError(f'{name} not on PATH or marked as executable.')


def aligner(reference_in, query_in, prefix, id_threshold=0.93, threads=4):  # pragma: no cover
    """
    Align query to the reference and extract qc metrics.

    Parameters
    ----------
    reference_in : str
        reference FASTA location
    query_in : str
        query FASTA location.
    prefix : str
        Prefix where to store the results.
    id_threshold : float, optional
        Identity threshold after aligment that must be achieved. The default is 0.93.
    threads : int, optional
        Number of threads to use. The default is 4.

    Returns
    -------
    datframe,
        result metrics
    """
    # Files exist?
    assert os.path.isfile(reference_in)
    assert os.path.isfile(query_in)

    # handle path / prefix input
    out_dir = os.path.dirname(os.path.abspath(prefix))
    # dir is not created when prefix ends with /
    # this is a fix for the bug
    if prefix.endswith('/'):
        out_dir = prefix
    if prefix != "" and not prefix.endswith('/'):
        prefix += "_"
    file_prefix = os.path.basename(prefix)
    if not os.path.exists(out_dir):
        print("Creating output directory...")
        os.makedirs(out_dir, exist_ok=True)

    print(f"Writing files to: {out_dir}")
    print(f"Using the prefix: {file_prefix}_* to store results.")

    # remove white spaces from fasta files
    reference_tmp = sequence.preprocess(reference_in)
    query_tmp = sequence.preprocess(query_in)

    # check reference fasta
    statistics.count_reference_sequences(reference_tmp)

    # perform initial sequence check
    query_tmp, evaluation, invalid_ids = \
        statistics.split_valid_sequences(query_tmp, reference_tmp, id_threshold=id_threshold)
    # if none of the sequences pass the qc filter, exit.
    # else just perform the alignment witht he seqs passing qc
    if evaluation == "all_invalid":
        print("None of the sequences can pass the identity threshold. No alignment done.")
        print("Exiting president alignment.")
        sys.exit()
    else:
        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")

    # perform alignment with pblat
    alignment_file = alignment.pblat(threads, reference_tmp, query_tmp, verbose=1)

    # parse statistics from file
    metrics = statistics.nucleotide_identity(query_tmp, alignment_file, id_threshold)

    with screed.open(reference_in) as seqfile:
        refseq = [i for i in seqfile][0]

    # add invalid
    if len(invalid_ids) > 0:
        invalid_df = pd.DataFrame({"ID": invalid_ids})
        invalid_df["passed_initial_qc"] = False
        invalid_df["aligned"] = False
        metrics = pd.concat([metrics, invalid_df]).reset_index(drop=True)

    # store sequences
    writer.write_sequences(query_in, metrics, prefix)

    # store reference data
    metrics["reference_length"] = len(refseq.sequence)
    metrics["reference"] = os.path.basename(reference_in)
    metrics["query"] = os.path.basename(query_in)
    # remove temporary files
    os.remove(alignment_file)
    os.remove(reference_tmp)
    os.remove(query_tmp)
    metrics.to_csv(os.path.join(out_dir, f"{file_prefix}report.tsv"), index=False, sep='\t')
    print(metrics)
    return metrics


def main():  # pragma: no cover
    """
    Presidents main function for sequence alignment.

    Returns
    -------
    None.

    """
    parser = argparse.ArgumentParser(description='Calculate pairwise nucleotide identity.')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome.')
    parser.add_argument('-q', '--query', required=True, help='Query genome.')
    parser.add_argument('-x', '--id_threshold', type=float, default=0.93,
                        help='ACGT nucleotide identity threshold after alignment (percentage). '
                             'A query sequence is reported as valid based on this threshold '
                             '(def: 0.93)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use for pblat.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Prefix to be used to store results and FASTA files.')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    aligner(args.reference, args.query, args.prefix, args.id_threshold, args.threads)


if __name__ == "__main__":
    main()
