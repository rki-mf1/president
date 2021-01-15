#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from president import alignment


fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')
fixtures_loc = '/home/hanjo/workspace/president/tests/fixtures/'

def test_pblat_simple_10MM():
    # read sample data
    reference = os.path.join(fixtures_loc, "100bp_0N_sample_reference.fasta")
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    
    # run pblat

    alignment_f = alignment.pblat(4, reference, query, verbose=1)


    exp_output = \


