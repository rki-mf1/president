#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from president import sequence
import filecmp
import os

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_remove_spaces():
    fin = os.path.join(fixtures_loc, "test_remove_spaces.fasta")
    expected = os.path.join(fixtures_loc, "test_remove_spaces_expected.fasta")
    fout = sequence.preprocess(fin)
    assert filecmp.cmp(fout, expected, shallow=True)


def test_to_valid_upper():
    dna_seq = "ACGTNXZY123!?acgt"
    exp_sequence = "ACGT" + "N" * len("NXZY123!?") + "ACGT"
    cm_sequence = sequence.to_valid_upper(dna_seq)
    assert cm_sequence == exp_sequence
