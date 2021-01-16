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


def test_unmapped():
    alignment_f = os.path.join(fixtures_loc, "test_unmapped.tsv")
    query = os.path.join(fixtures_loc, "test_unmapped.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, max_invalid=3000)

    exp_invalid = [False, True, True, True, True]
    exp_ident = [0.4063, 0.9953, 0.9952, 0.9947, 0.9796]
    exp_ambig_identity = [0.9991, 0.9994, 0.9993, 0.9989, 0.9993]
    exp_ambig_bases = [17742, 123, 121, 123, 590]
    exp_length_query = [29903, 29902, 29903, 29896, 29902]
    exp_IDS = ["FAO96286_barcode67/ARTIC/medaka_MN908947.3",
               "FAO96286_barcode07/ARTIC/medaka_MN908947.3",
               "FAO96286_barcode21/ARTIC/medaka_MN908947.3",
               "FAO96286_barcode33/ARTIC/medaka_MN908947.3",
               "FAO96286_barcode44/ARTIC/medaka_MN908947.3"]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_ident == metrics["Identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["Ambiguous Identity"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)
