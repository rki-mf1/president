"""President's main module.

Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage
-----
president --query tiny_test.masked_consensus.fasta \
    --reference NC_045512.2.fasta -x 3000 -p 8 -o report.tsv

Notes
-----
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line
(below) to ignore Ns

Requirements:

- conda install -y -c bioconda python=3.8 pblat screed pandas

ANI definition:

- https://pubmed.ncbi.nlm.nih.gov/17220447/

Spec
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
    None
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
    parser = argparse.ArgumentParser(description='Calculate pairwise nucleotide identity')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome')
    parser.add_argument('-q', '--query', required=True, help='Query genome')
    parser.add_argument('-x', '--max_invalid', type=int, default=3000,
                        help='Maximum number of non-ACTG positions allowed in query')
    parser.add_argument('-p', '--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file to write report')
    args = parser.parse_args()

    # Files exist?
    assert os.path.isfile(args.reference)
    assert os.path.isfile(args.query)

    # remove white spaces from fasta files
    reference = alignment.remove_spaces(args.reference)
    query = alignment.remove_spaces(args.query)

    # perform alignment with pblat
    alignment_file = alignment.pblat(args.threads, reference, query)

    # parse statistics from file
    metrics = statistics.nucleotide_identity(query, alignment_file, args.max_invalid)

    # remove temporary files
    os.remove(alignment_file)
    os.remove(reference)
    os.remove(query)
    metrics.to_csv(args.output, index=False, sep='\t')


if __name__ == "__main__":
    main()
