#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 10:15:48 2021

@author: hanjo
"""
import os

from president import statistics

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_statistics():
    alignment_f = os.path.join(fixtures_loc, "100bp_pblat_results.txt")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")

    ident, ident_non_canonical, non_canonical, length_query = \
        statistics.nucleotide_identity(query, alignment_f, max_invalid=3000)

    exp_ident = 0.95
    exp_ident_non_canonical = 0.95
    exp_non_canonical = 0
    exp_length_query = 100

    assert exp_ident == ident
    assert exp_ident_non_canonical == ident_non_canonical
    assert exp_non_canonical == non_canonical
    assert exp_length_query == length_query
