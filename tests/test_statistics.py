#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 10:15:48 2021

@author: hanjo
"""
import os
import numpy as np

import screed
import pytest
from president import statistics, alignment

fixtures_loc = os.path.join(os.path.dirname(__file__), "fixtures")


def test_statistics_nthreshold():
    query = os.path.join(fixtures_loc, "100bp_5N_05MM_sample_query.fasta")

    qc_stats_pass = statistics.summarize_query(query)
    statistics.qc_check(100, qc_stats_pass, id_threshold=0.93, n_threshold=0.05)

    qc_stats_fail = statistics.summarize_query(query)
    statistics.qc_check(100, qc_stats_fail, id_threshold=0.93, n_threshold=0.04)
    assert qc_stats_pass["qc_valid_pass_nthreshold"].iloc[0]
    assert not qc_stats_fail["qc_valid_pass_nthreshold"].iloc[0]


def test_statistics():
    alignment_f = os.path.join(fixtures_loc, "100bp_pblat_results.txt")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(100, qc_stats, id_threshold=0.93)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats, id_threshold=0.9)

    exp_ident = 0.95
    exp_ambig_identity = 0.95
    exp_ambig_bases = 0
    exp_length_query = 100
    exp_valid = True

    assert exp_ident == metrics["ACGT Nucleotide identity"].iloc[0]
    assert exp_valid == metrics["qc_post_aligned_all_valid"].iloc[0]
    assert (
        exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].iloc[0]
    )
    assert exp_ambig_bases == metrics["iupac_bases"].iloc[0]
    assert exp_length_query == metrics["length_query"].iloc[0]


def test_multi_statistics():
    alignment_f = os.path.join(fixtures_loc, "100pb_pblat_results_multi.txt")
    query = os.path.join(fixtures_loc, "100bp_multi.fasta")

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(100, qc_stats, id_threshold=0.91)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats, id_threshold=0.91)
    metrics = metrics.sort_values(by="query_name")

    # sorting gives first the one without Ns
    exp_invalid = [True, False]
    exp_ident = [0.95, 0.9]
    # 100 - 10 mismatches = 90; 100 - 5 Ns = 95 --> non-canocical
    exp_ambig_identity = np.round([0.95, 90 / 95], 4)
    exp_ambig_bases = [0, 5]
    exp_length_query = [100, 100]
    exp_IDS = ["100bp_0N_5MM", "100bp_5N_5MM"]

    assert np.all(exp_IDS == metrics["query_name"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["qc_post_aligned_all_valid"].values)
    assert np.all(
        exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values
    )
    assert np.all(exp_ambig_bases == metrics["iupac_bases"].values)
    assert np.all(exp_length_query == metrics["length_query"].values)


def test_unmapped():
    alignment_f = os.path.join(fixtures_loc, "test_unmapped.tsv")
    query = os.path.join(fixtures_loc, "test_unmapped.fasta")

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(29903, qc_stats, id_threshold=0.93)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats, id_threshold=0.9)
    metrics = metrics[metrics["qc_post_aligned"]]

    exp_IDS = np.array(
        [
            "FAO96286_barcode67/ARTIC/medaka MN908947.3",
            "FAO96286_barcode07/ARTIC/medaka MN908947.3",
            "FAO96286_barcode21/ARTIC/medaka MN908947.3",
            "FAO96286_barcode33/ARTIC/medaka MN908947.3",
            "FAO96286_barcode44/ARTIC/medaka MN908947.3",
        ]
    )
    exp_idx = np.argsort(exp_IDS)
    exp_IDS = exp_IDS[exp_idx]
    exp_invalid = np.array([False, True, True, True, True])[exp_idx]
    exp_ident = np.array([0.4063, 0.9953, 0.9952, 0.9947, 0.9796])[exp_idx]
    exp_ambig_identity = np.array([0.9991, 0.9994, 0.9993, 0.9989, 0.9993])[exp_idx]
    exp_ambig_bases = np.array([17742, 123, 121, 123, 590])[exp_idx]
    exp_length_query = np.array([29903, 29902, 29903, 29896, 29902])[exp_idx]

    metrics = metrics.sort_values(by="query_name", ascending=True)
    assert np.all(exp_IDS == metrics["query_name"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["qc_post_aligned_all_valid"].values)
    assert np.all(
        exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values
    )
    assert np.all(exp_ambig_bases == metrics["iupac_bases"].values)
    assert np.all(exp_length_query == metrics["length_query"].values)


def test_global_identity():
    alignment_f = os.path.join(fixtures_loc, "global_identity.tsv")
    query = os.path.join(fixtures_loc, "global_identity.fasta")

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(query, qc_stats, id_threshold=0.5)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats, id_threshold=0.5)

    exp_invalid = [True]
    exp_ident = [0.5]
    exp_ambig_identity = [0.5]
    exp_ambig_bases = [0]
    exp_length_query = [50]
    exp_IDS = ["query"]

    assert np.all(exp_IDS == metrics["query_name"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["qc_post_aligned_all_valid"].values)
    assert np.all(
        exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values
    )
    assert np.all(exp_ambig_bases == metrics["iupac_bases"].values)
    assert np.all(exp_length_query == metrics["length_query"].values)


def test_get_largest_N_gap():
    seq1 = "NNATCGACTAGTCTGAC"
    seq2 = "NNATCGANNNCTAGTCTGAC"
    seq3 = "NNATCGACTANNNGTCTGACNNNN"
    exp_n1 = 2
    exp_n2 = 3
    exp_n3 = 4
    ngap1 = statistics.get_largest_N_gap(seq1)
    ngap2 = statistics.get_largest_N_gap(seq2)
    ngap3 = statistics.get_largest_N_gap(seq3)

    assert exp_n1 == ngap1
    assert exp_n2 == ngap2
    assert exp_n3 == ngap3


def test_repeated_pblat():
    alignment_f = os.path.join(fixtures_loc, "test_repeated_pblat.tsv")
    query = os.path.join(fixtures_loc, "test_repeated_pblat.fasta")

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(29903, qc_stats, id_threshold=0.93)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats, id_threshold=0.05)

    exp_invalid = [False, True, True, True, True]
    exp_ident = [0.4063, 0.9953, 0.9952, 0.9947, 0.9796]
    exp_ambig_identity = [0.9991, 0.9994, 0.9993, 0.9989, 0.9993]
    exp_ambig_bases = [17742, 123, 121, 123, 590]
    exp_length_query = [29903, 29902, 29903, 29896, 29902]
    exp_IDS = [
        "FAO96286_barcode67/ARTIC/medaka MN908947.3",
        "FAO96286_barcode07/ARTIC/medaka MN908947.3",
        "FAO96286_barcode21/ARTIC/medaka MN908947.3",
        "FAO96286_barcode33/ARTIC/medaka MN908947.3",
        "FAO96286_barcode44/ARTIC/medaka MN908947.3",
    ]

    assert np.all(exp_IDS == metrics["query_name"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["qc_post_aligned_all_valid"].values)
    assert np.all(
        exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values
    )
    assert np.all(exp_ambig_bases == metrics["iupac_bases"].values)
    assert np.all(exp_length_query == metrics["length_query"].values)


def test_split_valid_sequences():
    reference = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    # this query has 5N / 10N with 5MM to the reference
    query = os.path.join(fixtures_loc, "100bp_5N_05MM_sample_query_multi.fasta")

    # all good
    # no id threshold
    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats, id_threshold=0.0, n_threshold=1.0)
    outfile1, case1, ids1 = statistics.split_valid_sequences(query, qc_stats)

    # all invalid
    # > 5 MM
    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats, id_threshold=0.95, n_threshold=0.0)
    outfile2, case2, ids2 = statistics.split_valid_sequences(query, qc_stats)

    # mixed
    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats, id_threshold=0.95, n_threshold=0.05)
    outfile3, case3, ids3 = statistics.split_valid_sequences(query, qc_stats)

    assert case1 == "all_valid"
    assert ids1 == []

    assert case2 == "all_invalid"
    assert np.all(ids2 == ["100bp_5N_5MM_reference", "100bp_10N_5MM_reference"])

    assert case3 == "mixed"
    assert ids3 == ["100bp_10N_5MM_reference"]


def test_split_valid_sequences_uneven():
    reference = os.path.join(fixtures_loc, "101bp_5N_5MM_5G_ref.fasta")
    query = os.path.join(fixtures_loc, "101bp_5N_5MM_5G_query_3seqs.fasta")

    # mixed
    # 0.9 * 101 = 90.9
    # 0.94 * 101 = 94.94
    # 0.95 * 101 = 95.949
    # we have 101 - 5 = 96 = 96 / 101 = 0.95
    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats, id_threshold=0.95, n_threshold=1.0)
    outfile1, case1, ids1 = statistics.split_valid_sequences(query, qc_stats)

    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats, id_threshold=0.96, n_threshold=0.04)
    outfile2, case2, ids2 = statistics.split_valid_sequences(query, qc_stats)

    assert case1 == "all_valid"
    assert case2 == "mixed"


def test_number_reference_check_fail():
    reference = os.path.join(fixtures_loc, "100bp_multi.fasta")

    with pytest.raises(ValueError):
        statistics.count_sequences(reference)


def test_number_reference_check_pass():
    reference = os.path.join(fixtures_loc, "100bp_0N_sample_reference.fasta")
    statistics.count_sequences(reference)
    assert True


def test_estimate_max_invalid_even():
    reference = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    exp_maxinvalid = 5
    maxinvalid = statistics.estimate_max_invalid(reference, id_threshold=0.95)
    assert exp_maxinvalid == maxinvalid


def test_estimate_max_invalid_uneven():
    reference = os.path.join(fixtures_loc, "101bp_5N_5MM_5G_ref.fasta")
    # total - valid -> invalid
    exp_maxinvalid = (101 - 96) * 1.0
    maxinvalid = statistics.estimate_max_invalid(reference, id_threshold=0.95)
    assert exp_maxinvalid == maxinvalid


def test_init_UPAC_dictionary():
    iupac_dict = statistics.init_UPAC_dictionary()
    iupac_dict_acgt = statistics.init_UPAC_dictionary(acgt=True)

    assert "A" not in iupac_dict
    assert "C" not in iupac_dict
    assert "T" not in iupac_dict
    assert "G" not in iupac_dict

    assert "A" in iupac_dict_acgt
    assert "C" in iupac_dict_acgt
    assert "T" in iupac_dict_acgt
    assert "G" in iupac_dict_acgt


def test_nucleotide_counts():
    query = os.path.join(fixtures_loc, "test_UPAC_statistics.fasta")
    with screed.open(query) as seqfile:
        results = [statistics.count_nucleotides(seq.sequence) for seq in seqfile]
    exp_tuple = [
        (0, 0, 60, 0),  # only X
        (4, 1, 0, 1),  # ACGT + N
        (96, 1, 0, 0),  # ACGT.... + "."]
        (4, 11, 3, 1),  # ACGT + "RYSWKMBDHVN" + ?!-]
        (96, 0, 1, 0),
    ]  # ACGT.... + "-"]
    assert exp_tuple == results


def test_not_aligned():
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    query = os.path.join(fixtures_loc, "test_unmapped.fasta")

    # entries in fasta
    # one is invalid
    exp_rows = 6

    # run pblat
    qc_stats = statistics.summarize_query(query)
    statistics.qc_check(reference, qc_stats)

    alignment_f = alignment.pblat(4, reference, query)
    metrics = statistics.nucleotide_identity(alignment_f, qc_stats)

    assert exp_rows == metrics.shape[0]
    assert metrics["qc_post_aligned"].sum() == exp_rows - 1
