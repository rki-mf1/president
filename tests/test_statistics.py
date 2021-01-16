#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 10:15:48 2021

@author: hanjo
"""
import os
import numpy as np

from president import statistics

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_statistics():
    alignment_f = os.path.join(fixtures_loc, "100bp_pblat_results.txt")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, max_invalid=3000)

    exp_ident = 0.95
    exp_ambig_identity = 0.95
    exp_ambig_bases = 0
    exp_length_query = 100
    exp_valid = True

    assert exp_ident == metrics["Identity"].iloc[0]
    assert exp_valid == metrics["Valid"].iloc[0]
    assert exp_ambig_identity == metrics["Ambiguous Identity"].iloc[0]
    assert exp_ambig_bases == metrics["Ambiguous Bases"].iloc[0]
    assert exp_length_query == metrics["Query Length"].iloc[0]


def test_multi_statistics():
    alignment_f = os.path.join(fixtures_loc, "100pb_pblat_results_multi.txt")
    query = os.path.join(fixtures_loc, "100bp_multi.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, max_invalid=4)
    metrics = metrics.sort_values(by="ID")

    # sorting gives first the one without Ns
    exp_invalid = [True, False]
    exp_ident = [0.95, 0.9]
    # 100 - 10 mismatches = 90; 100 - 5 Ns = 95 --> non-canocical
    exp_ambig_identity = np.round([0.95, 90 / 95], 4)
    exp_ambig_bases = [0, 5]
    exp_length_query = [100, 100]
    exp_IDS = ["100bp_0N_5MM", "100bp_5N_5MM"]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_ident == metrics["Identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["Ambiguous Identity"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)
