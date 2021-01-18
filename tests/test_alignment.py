#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import filecmp
import pytest
import tempfile
import subprocess
import time

import numpy as np
import pandas as pd
import screed

from president import alignment

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_pblat_simple_10MM():
    # read sample data
    reference = os.path.join(fixtures_loc, "100bp_0N_sample_reference.fasta")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    exp_matches = 95
    exp_mismatches = 5
    exp_gaps = 0
    exp_Ns = 0

    # run pblat
    alignment_f = alignment.pblat(4, reference, query, verbose=1)

    alignments = pd.read_csv(alignment_f, header=None, sep='\t', skiprows=5)
    alignments.columns = \
        ['Matches', 'Mismatches', 'RepMatch', 'Ns', 'QGapCount',
         'QGapBases', 'TGapCount', 'TGapBases', 'Strand',
         'QName', 'QSize', 'QStart', 'QEnd', 'TName', 'TSize',
         'TStart', 'TEnd', 'BlockCount', 'BlockSizes',
         'QStarts', 'TStarts']

    assert exp_matches == alignments["Matches"].iloc[0]
    assert exp_mismatches == alignments["Mismatches"].iloc[0]
    assert exp_gaps == alignments["QGapCount"].iloc[0]
    assert exp_Ns == alignments["Ns"].iloc[0]


def test_parse_alignment():
    alignment_f = os.path.join(fixtures_loc, "100bp_pblat_results.txt")
    metrics = alignment.parse_alignment(alignment_f)
    assert metrics.shape == (1, 21)


def test_parse_alignment_multi():
    alignment_f = os.path.join(fixtures_loc, "100pb_pblat_results_multi.txt")
    metrics = alignment.parse_alignment(alignment_f)
    assert metrics.shape == (2, 21)


def test_remove_spaces():
    fin = os.path.join(fixtures_loc, "test_remove_spaces.fasta")
    expected = os.path.join(fixtures_loc, "test_remove_spaces_expected.fasta")
    fout = alignment.remove_spaces(fin)
    assert filecmp.cmp(fout, expected, shallow=True)


def test_pblat_simple_gaps():
    reference = os.path.join(fixtures_loc, "50kbp_1000N_500MM_500G_ref.fasta")
    query = os.path.join(fixtures_loc, "50kbp_1000N_500MM_500G_query.fasta")
    alignment_f = alignment.pblat(4, reference, query)
    metrics = alignment.parse_alignment(alignment_f)
    exp_gaps = 1
    exp_gapbases = 500

    assert exp_gaps == metrics["QGapCount"].iloc[0]
    assert exp_gapbases == metrics["QGapBases"].iloc[0]


@pytest.mark.skip
def test_large_real_data():
    # download data from EBI and run large test
    url = "ftp://ftp.ebi.ac.uk/pub/databases/covid19dataportal/viral_sequences/sequences/"
    filename = "sequences_fasta_20210105-1030.fa.gz"
    url += filename

    # init temp file
    _, tmpfile = tempfile.mkstemp(suffix=".fa.gz")
    #  tmpfile = os.path.join(tmpfile, filename)

    # download sequences
    wget_cmd = f"wget {url} -O {tmpfile}"
    print("downloading files...")
    _ = subprocess.check_output(wget_cmd, shell=True)

    # extract sequences
    cmd = f"gzip -d {tmpfile}"
    print("unzipping files...")
    _ = subprocess.check_output(cmd, shell=True)

    print("sampling files ...")
    # down sample input sequence
    seqs = []
    n = 10000
    # open the file and add more sequences than needed
    with screed.open(tmpfile.replace(".gz", "")) as seqfile:
        for idx, qry in enumerate(seqfile):
            if idx > 3*n:
                break
            seqs.append(qry)

    # init temp file for sampled fasta
    _, tmpfile_sampled = tempfile.mkstemp(suffix=".fa")
    # write the sampled fasta entries to file
    np.random.seed(42)
    idx = np.random.choice(np.arange(0, n), n, replace=False)
    with open(tmpfile_sampled, "w") as sample_out:
        for ii in idx:
            sample_out.write(f">{seqs[ii].name}\n")
            sample_out.write(f"{seqs[ii].sequence}\n")

    # remove gz extension for calling alignments on the extracted fasta
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    query = tmpfile_sampled

    # perform alignment
    a = time.time()
    alignment_f = alignment.pblat(4, reference, query)
    metrics = alignment.parse_alignment(alignment_f)
    b = time.time()
    duration = np.round((b-a) / 60, 2)
    print(f"alignment of {n} real sequences took: {duration} minutes")
    exp_shape = (n, 21)
    try:
        os.remove(tmpfile)
        os.remove(tmpfile_sampled)
    except FileNotFoundError:
        print("Nothing to delete.")

    assert metrics.shape[0] > exp_shape[0] * 0.9
    assert duration <= 15


# @pytest.mark.benchmark(min_rounds=3)
# benchmark
def test_scale_alignments_10():
    reference = os.path.join(fixtures_loc, "30kbp_500N_1000MM_250G_ref_nseqs.fasta")
    query = os.path.join(fixtures_loc, "30kbp_500N_1000MM_250G_query_nseqs10.fasta")
    # alignment_f = benchmark(alignment.pblat, 4, reference, query)
    alignment_f = alignment.pblat(4, reference, query)
    metrics = alignment.parse_alignment(alignment_f)
    assert 10 == metrics.shape[0]


# @pytest.mark.benchmark(min_rounds=3)
# benchmark
def test_scale_alignments_100():
    reference = os.path.join(fixtures_loc, "30kbp_500N_1000MM_250G_ref_nseqs.fasta")
    query = os.path.join(fixtures_loc, "30kbp_500N_1000MM_250G_query_nseqs100.fasta")
    # alignment_f = benchmark(alignment.pblat, 4, reference, query)
    alignment_f = alignment.pblat(4, reference, query)
    metrics = alignment.parse_alignment(alignment_f)
    assert 100 == metrics.shape[0]
