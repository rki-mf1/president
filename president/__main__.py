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
# authors:
# RKI MF1;  Martin Hoelzer with great initial help of @phiweger (UKL Leipzig)
# HPI;      Fabio Malcher Miranda, Sven Giese, Alice Wittig

import os
import argparse
from shutil import which

import pandas as pd
import screed

from president import alignment
from president import statistics


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
    parser.add_argument('-p', '--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file to write report.')
    args = parser.parse_args()

    # Files exist?
    assert os.path.isfile(args.reference)
    assert os.path.isfile(args.query)

    # remove white spaces from fasta files
    reference = alignment.remove_spaces(args.reference)
    query = alignment.remove_spaces(args.query)

    # check reference fasta
    statistics.count_reference_sequences(reference)

    # perform initial sequence check
    query, evaluation, invalid_ids = \
        statistics.split_valid_sequences(query, reference, id_threshold=args.id_threshold)

    # if none of the sequences pass the qc filter, exit.
    # else just perform the alignment witht he seqs passing qc
    if evaluation == "all_invalid":
        assert "None of the sequences can pass the identity threshold. No alignment done."
    else:
        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")

    # perform alignment with pblat
    print(reference)
    print(query)
    alignment_file = alignment.pblat(args.threads, reference, query)

    # parse statistics from file
    metrics = statistics.nucleotide_identity(query, alignment_file, args.id_threshold)

    with screed.open(reference) as seqfile:
        refseq = [i for i in seqfile][0]

    # add invalid
    if len(invalid_ids) > 0:
        invalid_df = pd.DataFrame({"ID": invalid_ids})
        invalid_df["passed_initial_qc"] = False
        invalid_df["aligned"] = False
        metrics = pd.concat([metrics, invalid_df]).reset_index(drop=True)

    # store reference data
    metrics["reference_length"] = len(refseq.sequence)
    # remove temporary files
    os.remove(alignment_file)
    os.remove(reference)
    os.remove(query)
    metrics.to_csv(args.output, index=False, sep='\t')
    print(metrics)


if __name__ == "__main__":
    main()
