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
