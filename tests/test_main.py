import os
import tempfile

import pytest

from president import __main__ as pm

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_isavailable():
    with pytest.raises(ValueError):
        pm.is_available("asdhasdasdasd_tool")


def test_aligner():
    query = os.path.join(fixtures_loc, "test_UPAC_statistics.fasta")
    reference = os.path.join(fixtures_loc, "test_UPAC_ref.fasta")

    tmpfile = tempfile.mkstemp(suffix=".csv")[1]
    metrics = pm.aligner(reference_in=reference, query_in=query, id_threshold=0.0, threads=4,
                         prefix=tmpfile)

    exp_rows = 5
    os.remove(tmpfile)

    # TODO:
    # very rudimentary test
    # check for values in dataframe, output files, sequences in the output files ...
    assert exp_rows == metrics.shape[0]
