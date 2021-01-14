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

- conda install -y -c bioconda python=3.8 mafft screed

ANI definition:

- https://pubmed.ncbi.nlm.nih.gov/17220447/

Spec:

> Sequence identity to a reference genome sequence (e.g. NC_045512.2) should be calculated at the nucleotide level relative to the entire length of the reference. Only informative nucleotides (A,T,G,C) are considered identical to each other.
'''


import argparse
import os
import subprocess
import tempfile

import screed


def is_tool(name):
    '''Check whether `name` is on PATH and marked as executable.'''
    from shutil import which
    return which(name) is not None


def calculate_nucleotide_identity(path, max_invalid):
    '''
    Calculate nucleotide ident from a 2-sequence MSA
    '''
    with screed.open(path) as file:
        ref = next(file).sequence
        qry = next(file).sequence
    
    
    # Consider only sites where the query has non-ACTG characters
    # Metric issue #2 (B)
    non_canonical = sum([1 for i in qry if i not in 'ACTG'])
    if non_canonical > max_invalid:
        raise ValueError('Too many non-canonical nucleotides, abort!')
    
    same = 0
    for r, q in zip(ref, qry):
        if (q in 'ACTG' and r == q):
            same += 1
    
    # Metric issue #2 (A)
    # Ns in the query count as mismatch
    ident = round(same / len(qry), 4)
    
    # Ns in the query don't count
    ident_non_canonical = round(same / (len(qry) - non_canonical), 4)
    return ident, ident_non_canonical, non_canonical, len(qry)


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


    # Mafft installed?
    if not is_tool('mafft'):
        raise ValueError('Mafft aligner not on PATH or marked as executable.')

    # Files exist?
    assert os.path.isfile(args.reference) 
    assert os.path.isfile(args.query) 

    print('Running Mafft ...')
    _, path = tempfile.mkstemp()
    cmd = f'linsi --quiet --addfragments {args.query} --thread {args.threads} --adjustdirectionaccurately --keeplength --anysymbol {args.reference} > {path}'
    _ = subprocess.check_output(cmd, shell=True)


    ident, identN, nc, len_  = calculate_nucleotide_identity(
        path, args.max_invalid)
    os.remove(path)
    print(f'{ident} nucleotide identity')
    print(f'{identN} nucleotide identity (excluding non-ACTG)')
    print(f'{nc} non-ACTG symbols out of {len_}')


if __name__ == "__main__":
    main()
