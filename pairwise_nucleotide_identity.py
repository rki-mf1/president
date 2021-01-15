#! /usr/bin/env python

# authors: 
# RKI MF1;  Martin Hoelzer with great initial help of @phiweger (UKL Leipzig)
# HPI;      Fabio Malcher Miranda, Sven Giese, Alice Wittig

'''
Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage:

python pairwise_nucleotide_identity.py --query tiny_test.masked_consensus.fasta --reference NC_045512.2.fasta -x 3000 -p 8 -o report.tsv

Notes:

- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line
(below) to ignore Ns

Requirements:

- conda install -y -c bioconda python=3.8 pblat=2.5 screed pandas

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
import pandas as pd


def is_tool(name):
    '''Check whether `name` is on PATH and marked as executable.'''
    from shutil import which
    return which(name) is not None


def calculate_nucleotide_identity(query, alignments, max_invalid):
    # Read the alignment(s) with pandas
    # Pandas can be replaced with split to reduce dependencies
    alignments = pd.read_csv(
        alignments, header=None, sep='\t', skiprows=5)
    labels = ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
              'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
              'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
              'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
              'QStarts', 'TStarts']
    alignments.columns = labels


    '''
    Calculate nucleotide ident from a 2-sequence MSA
    '''
    query_ids = []
    valid_sequences = []
    ambiguous_bases = []
    identities = []
    ambiguous_identities = []
    query_lengths = []
    index = 0
    with screed.open(query) as seqfile:
        for qry in seqfile:
            query_ids.append(qry.name)
            # Consider only sites where the query has non-ACTG characters
            # Metric issue #2 (B)
            non_canonical = sum([1 for i in qry.sequence if i not in 'ACTG'])
            ambiguous_bases.append(non_canonical)
            if non_canonical > max_invalid:
                valid_sequences.append(False)
            else:
                valid_sequences.append(True)
            # Metric issue #2 (A)
            # Ns in the query count as mismatch
            identities.append(round(alignments.at[index, 'Matches'] / len(qry.sequence), 4))
            # Ns in the query don't count
            ambiguous_identities.append(round(alignments.at[index, 'Matches'] / (len(qry.sequence) - non_canonical), 4))
            query_lengths.append(len(qry.sequence))
            index = index + 1

    metrics = pd.DataFrame({
        'ID': query_ids,
        'Valid': valid_sequences,
        'Identity': identities,
        'Ambiguous Identity': ambiguous_identities,
        'Ambiguous Bases': ambiguous_bases,
        'Query Length': query_lengths
    })
    return metrics


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
    parser.add_argument('-o', '--output', required=True,
        help='Output TSV file to write report')
    args = parser.parse_args()


    # pblat installed?
    if not is_tool('pblat'):
        raise ValueError('pblat aligner not on PATH or marked as executable.')

    # Files exist?
    assert os.path.isfile(args.reference) 
    assert os.path.isfile(args.query) 

    print('Running pblat ...')
    _, alignments = tempfile.mkstemp()
    cmd = f'pblat -threads={args.threads} {args.reference} {args.query} {alignments}'
    _ = subprocess.check_output(cmd, shell=True)


    metrics  = calculate_nucleotide_identity(
        args.query, alignments, args.max_invalid)
    os.remove(alignments)
    metrics.to_csv(args.output, index=False, sep='\t')


if __name__ == "__main__":
    main()
