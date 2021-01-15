#! /usr/bin/env python

# author: @phiweger, @martinhoelzer

'''
Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage:

python pairwise_nucleotide_identity.py --query tiny_test.masked_consensus.fasta --reference NC_045512.2.fasta -x 3000 -p 8

Notes:

- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line
(below) to ignore Ns

Requirements:

- conda install -y -c bioconda python=3.8 pblat screed pandas

ANI definition:

- https://pubmed.ncbi.nlm.nih.gov/17220447/

Spec:

> Sequence identity to a reference genome sequence (e.g. NC_045512.2) should be calculated at the nucleotide level relative to the entire length of the reference. Only informative nucleotides (A,T,G,C) are considered identical to each other.
'''

import os
import argparse

from president import alignment
from president import statistics


def is_available(name="pblat"):
    """
    Check whether `name` is on PATH and marked as executable.'

    Parameters
    ----------
    name : str
        tool name executable (default: pblat).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    '''Check whether `name` is on PATH and marked as executable.'''
    from shutil import which
    
    # tool installed
    if not which(name):
        raise ValueError(f'{name} not on PATH or marked as executable.')


def main():
    parser = argparse.ArgumentParser(
        description='Calculate pairwise nucleotide identity')
    parser.add_argument('-r', '--reference', required=True,
        help='Reference genome')
    parser.add_argument('-q', '--query', required=True,
        help='Query genome')
    parser.add_argument('-x', '--max_invalid', type=int, default=3000,
        help='Maximum number of non-ACTG positions allowed in query')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads to use')
    args = parser.parse_args()


    # Files exist?
    assert os.path.isfile(args.reference) 
    assert os.path.isfile(args.query) 

    # perform alignment with pblat
    alignment_file = alignment.pblat(args.threads, args.reference, args.query)

    # parse statistics from file
    ident, identN, nc, len_  = \
        statistics.nucleotide_identity(args.query, alignment_file, 
                                       args.max_invalid)
    
    # remove temporary result file
    os.remove(alignment_file)
    
    print(f"Results comparing: {args.reference} and {args.query}")
    print(f'{ident} nucleotide identity')
    print(f'{identN} nucleotide identity (excluding non-ACTG)')
    print(f'{nc} non-ACTG symbols out of {len_}')


if __name__ == "__main__":
    main()