import os
import tempfile

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
    president_df = pm.aligner(reference_in=reference, query_in=query, id_threshold=0.0,
                              threads=4,
                              prefix=tmpfile)
    os.remove(tmpfile)
    assert president_df.shape == (19, 27)
