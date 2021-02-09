import os
import tempfile

import numpy as np
import pytest

from president import __main__ as pm

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_isavailable():
    with pytest.raises(ValueError):
        pm.is_available("asdhasdasdasd_tool")


def test_aligner_combined():
    query = os.path.join(fixtures_loc, "test_combined.fasta")
    # has 19 sequences
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    tmppath = os.path.splitext(tmpfile)[0]
    president_df = pm.aligner(reference_in=reference, query_in_raw=query, id_threshold=0.0,
                              threads=4, path_out=tmppath)
    os.remove(tmpfile)
    assert president_df.shape == (19, 27)


def test_empty_alignment():
    query = os.path.join(fixtures_loc, "nnnnnnn.fasta")
    # has only failing sequence
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    tmppath = os.path.splitext(tmpfile)[0]
    president_df = pm.aligner(reference_in=reference, query_in_raw=query, id_threshold=0.9,
                              threads=4, path_out=tmppath)
    os.remove(tmpfile)
    assert not president_df["qc_all_valid"].iloc[0]
    assert president_df.shape == (1, 27)


def test_empty_fasta():
    query = os.path.join(fixtures_loc, "empty.fasta")
    # has only failing sequence
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    president_df = pm.aligner(reference_in=reference, query_in_raw=query, id_threshold=0.9,
                              threads=4, prefix_in=tmpfile)
    os.remove(tmpfile)
    assert president_df.shape == (1, 27)


def test_same_format():
    # test if aligned and non aligned outputs ahve the same format.
    query = os.path.join(fixtures_loc, "test_combined.fasta")
    # has 19 sequences
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    tmppath = os.path.splitext(tmpfile)[0]
    president_df = pm.aligner(reference_in=reference, query_in_raw=query, id_threshold=0.0,
                              threads=4, path_out=tmppath)
    os.remove(tmpfile)

    query = os.path.join(fixtures_loc, "nnnnnnn.fasta")
    # has only failing sequence
    reference = os.path.join(fixtures_loc, "NC_045512.2.fasta")
    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    tmppath = os.path.splitext(tmpfile)[0]
    president_df2 = pm.aligner(reference_in=reference, query_in_raw=query, id_threshold=0.9,
                               threads=4, path_out=tmppath)
    os.remove(tmpfile)

    assert np.all(president_df.columns == president_df2.columns)
