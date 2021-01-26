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

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_statistics():
    alignment_f = os.path.join(fixtures_loc, "100bp_pblat_results.txt")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, id_threshold=0.93)

    exp_ident = 0.95
    exp_ambig_identity = 0.95
    exp_ambig_bases = 0
    exp_length_query = 100
    exp_valid = True

    assert exp_ident == metrics["ACGT Nucleotide identity"].iloc[0]
    assert exp_valid == metrics["Valid"].iloc[0]
    assert exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].iloc[0]
    assert exp_ambig_bases == metrics["Ambiguous Bases"].iloc[0]
    assert exp_length_query == metrics["Query Length"].iloc[0]


def test_multi_statistics():
    alignment_f = os.path.join(fixtures_loc, "100pb_pblat_results_multi.txt")
    query = os.path.join(fixtures_loc, "100bp_multi.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, id_threshold=0.91)
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
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)


def test_unmapped():
    alignment_f = os.path.join(fixtures_loc, "test_unmapped.tsv")
    query = os.path.join(fixtures_loc, "test_unmapped.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, id_threshold=0.93)
    metrics = metrics[metrics["aligned"]]

    exp_invalid = [False, True, True, True, True]
    exp_ident = [0.4063, 0.9953, 0.9952, 0.9947, 0.9796]
    exp_ambig_identity = [0.9991, 0.9994, 0.9993, 0.9989, 0.9993]
    exp_ambig_bases = [17742, 123, 121, 123, 590]
    exp_length_query = [29903, 29902, 29903, 29896, 29902]
    exp_IDS = ["FAO96286_barcode67/ARTIC/medaka MN908947.3",
               "FAO96286_barcode07/ARTIC/medaka MN908947.3",
               "FAO96286_barcode21/ARTIC/medaka MN908947.3",
               "FAO96286_barcode33/ARTIC/medaka MN908947.3",
               "FAO96286_barcode44/ARTIC/medaka MN908947.3"]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)


def test_global_identity():
    alignment_f = os.path.join(fixtures_loc, "global_identity.tsv")
    query = os.path.join(fixtures_loc, "global_identity.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, id_threshold=0.5)

    exp_invalid = [True]
    exp_ident = [0.5]
    exp_ambig_identity = [0.5]
    exp_ambig_bases = [0]
    exp_length_query = [50]
    exp_IDS = ["query"]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)


def test_repeated_pblat():
    alignment_f = os.path.join(fixtures_loc, "test_repeated_pblat.tsv")
    query = os.path.join(fixtures_loc, "test_repeated_pblat.fasta")

    metrics = statistics.nucleotide_identity(query, alignment_f, id_threshold=0.93)

    exp_invalid = [False, True, True, True, True]
    exp_ident = [0.4063, 0.9953, 0.9952, 0.9947, 0.9796]
    exp_ambig_identity = [0.9991, 0.9994, 0.9993, 0.9989, 0.9993]
    exp_ambig_bases = [17742, 123, 121, 123, 590]
    exp_length_query = [29903, 29902, 29903, 29896, 29902]
    exp_IDS = ["FAO96286_barcode67/ARTIC/medaka MN908947.3",
               "FAO96286_barcode07/ARTIC/medaka MN908947.3",
               "FAO96286_barcode21/ARTIC/medaka MN908947.3",
               "FAO96286_barcode33/ARTIC/medaka MN908947.3",
               "FAO96286_barcode44/ARTIC/medaka MN908947.3"]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)


def test_split_valid_sequences():
    reference = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    query = os.path.join(fixtures_loc, "100bp_5N_05MM_sample_query_multi.fasta")

    # all good
    outfile1, case1, ids1 = statistics.split_valid_sequences(query, reference, id_threshold=0.0)

    # all invalid
    outfile2, case2, ids2 = statistics.split_valid_sequences(query, reference, id_threshold=1)

    # mixed
    outfile3, case3, ids3 = statistics.split_valid_sequences(query, reference, id_threshold=0.95)

    assert case1 == "all_valid"
    assert ids1 == []

    assert case2 == "all_invalid"
    assert ids2 == ['100bp_5N_5MM_reference', '100bp_10N_5MM_reference']

    assert case3 == "mixed"
    assert ids3 == ['100bp_10N_5MM_reference']


def test_split_valid_sequences_uneven():
    reference = os.path.join(fixtures_loc, "101bp_5N_5MM_5G_ref.fasta")
    query = os.path.join(fixtures_loc, "101bp_5N_5MM_5G_query_3seqs.fasta")

    # mixed
    # 0.93 * 101 = 93.93
    # 0.94 * 101 = 94.94
    # 0.95 * 101 = 95.949
    # we have 101 - 5 = 96 = 96 / 101 = 0.95
    outfile1, case1, ids1 = statistics.split_valid_sequences(query, reference, id_threshold=0.95)
    outfile2, case2, ids2 = statistics.split_valid_sequences(query, reference, id_threshold=0.96)

    assert case1 == "all_valid"
    assert case2 == "mixed"


def test_number_reference_check_fail():
    reference = os.path.join(fixtures_loc, "100bp_multi.fasta")

    with pytest.raises(ValueError):
        statistics.count_reference_sequences(reference)


def test_number_reference_check_pass():
    reference = os.path.join(fixtures_loc, "100bp_0N_sample_reference.fasta")
    statistics.count_reference_sequences(reference)
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
    exp_tuple = [(0, 0, 60),  # only X
                 (4, 1, 0),  # ACGT + N
                 (96, 1, 0),  # ACGT.... + "."]
                 (4, 11, 3),  # ACGT + "RYSWKMBDHVN" + ?!-]
                 (96, 0, 1)]  # ACGT.... + "-"]
    assert exp_tuple == results


def test_not_aligned():
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    query = os.path.join(fixtures_loc, "test_unmapped.fasta")

    # entries in fasta
    # one is invalid
    exp_rows = 6

    # run pblat
    alignment_f = alignment.pblat(4, reference, query, verbose=1)
    metrics = statistics.nucleotide_identity(query, alignment_f)

    assert exp_rows == metrics.shape[0]
    assert metrics["aligned"].sum() == exp_rows - 1


def test_queries_all_invalid():
    query = os.path.join(fixtures_loc, "nnnnnnn_2files.fasta")
    n_seqs = 2
    metrics = statistics.metrics_all_invalid(query, n_seqs)

    exp_IDS = [
        ("NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, "
            "complete genome"),
        "seq 2"]
    exp_invalid = [False, False]
    exp_ident = ["", ""]
    exp_ambig_identity = ["", ""]
    exp_non_ACTGN_identity = ["", ""]
    exp_ambig_bases = [4220, 4220]
    exp_length_query = [29693, 29693]
    exp_query_ACTG = [25473, 25473]
    exp_query_IUPAC = [4220, 4220]
    exp_query_non_IUPAC = [0, 0]
    exp_aligned = [False, False]
    exp_passed_initial_qc = [False, False]

    assert np.all(exp_IDS == metrics["ID"])
    assert np.all(exp_invalid == metrics["Valid"].values)
    assert np.all(exp_ident == metrics["ACGT Nucleotide identity"].values)
    assert np.all(exp_ambig_identity == metrics["ACGT Nucleotide identity (ignoring Ns)"].values)
    assert np.all(exp_non_ACTGN_identity == metrics["ACGT Nucleotide identity (ignoring\
           non-ACGTNs)"].values)
    assert np.all(exp_ambig_bases == metrics["Ambiguous Bases"].values)
    assert np.all(exp_length_query == metrics["Query Length"].values)
    assert np.all(exp_query_ACTG == metrics["Query #ACGT"].values)
    assert np.all(exp_query_IUPAC == metrics["Query #IUPAC-ACGT"].values)
    assert np.all(exp_query_non_IUPAC == metrics["Query #non-IUPAC"].values)
    assert np.all(exp_aligned == metrics["aligned"].values)
    assert np.all(exp_passed_initial_qc == metrics["passed_initial_qc"].values)
