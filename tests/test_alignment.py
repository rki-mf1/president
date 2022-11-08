#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

import pandas as pd
from president import alignment

fixtures_loc = os.path.join(os.path.dirname(__file__), "fixtures")


def test_pblat_simple_10MM():
    # read sample data
    reference = os.path.join(fixtures_loc, "100bp_0N_sample_reference.fasta")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    exp_matches = 95
    exp_mismatches = 5
    exp_gaps = 0
    exp_Ns = 0

    # run pblat
    alignment_f = alignment.pblat(4, reference, query)

    alignments = pd.read_csv(alignment_f, header=None, sep="\t", skiprows=5)
    alignments.columns = [
        "Matches",
        "Mismatches",
        "RepMatch",
        "Ns",
        "QGapCount",
        "QGapBases",
        "TGapCount",
        "TGapBases",
        "Strand",
        "QName",
        "QSize",
        "QStart",
        "QEnd",
        "TName",
        "TSize",
        "TStart",
        "TEnd",
        "BlockCount",
        "BlockSizes",
        "QStarts",
        "TStarts",
    ]

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
