import tempfile
import os
import shutil

import screed
import pandas as pd
import hashlib

from president import writer

fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_get_filename():
    filename = "/home/user/dir/myfasta.fasta"
    basename = writer.get_filename(filename)
    exp_basename = "myfasta"

    assert exp_basename == basename


def test_write_fasta():
    query = os.path.join(fixtures_loc, "100bp_0N_05MM_sample_query.fasta")
    with screed.open(query) as seqfile:
        sequence = [i for i in seqfile][0]

    # write output
    tmpfile = tempfile.mkstemp(suffix=".fasta")[1]
    fout = open(tmpfile, "w")
    writer.write_fasta(fout, sequence)
    fout.close()

    # read it again
    with screed.open(tmpfile) as seqfile:
        sequence_parsed = [i for i in seqfile][0]

    os.remove(tmpfile)
    assert sequence.name == sequence_parsed.name
    assert sequence.description == sequence_parsed.description
    assert sequence.sequence == sequence_parsed.sequence


def test_write_sequences():
    query = os.path.join(fixtures_loc, "test_UPAC_statistics.fasta")

    metrics = pd.DataFrame()
    metrics["qc_post_aligned_all_valid"] = [False, False, True, False, False]
    metrics["aligned"] = [False, False, True, False, True]
    metrics["query_name"] = ["all_x", "acgtn_short", "aligns", "has_iuepac_2", "has_iuepac_3"]

    tmpdir = os.path.join(tempfile.mkdtemp(), "")
    writer.write_sequences(query, metrics, tmpdir, "mixed")
    prefix = ""
    valid_name = os.path.join(tmpdir, f"{prefix}valid.fasta")
    invalid_name = os.path.join(tmpdir, f"{prefix}invalid.fasta")

    with screed.open(valid_name) as seqfile:
        valid_seqs = [seq.name for seq in seqfile]

    with screed.open(invalid_name) as seqfile:
        invalid_seqs = [seq.name for seq in seqfile]

    shutil.rmtree(tmpdir)

    assert valid_seqs == ["aligns"]
    assert invalid_seqs == ['all_x', 'acgtn_short', 'has_iuepac_2', 'has_iuepac_3']


def test_trailing_whitespace():
    query = os.path.join(fixtures_loc, "test_UPAC_statistics.fasta")

    metrics = pd.DataFrame()
    metrics["Valid"] = [False, False, True, False, False]
    metrics["qc_post_aligned_all_valid"] = [False, False, True, False, True]
    metrics["query_name"] = ["all_x", "acgtn_short", "aligns", "has_iuepac_2", "has_iuepac_3"]

    tmpdir = os.path.join(tempfile.mkdtemp(), "")
    writer.write_sequences(query, metrics, tmpdir, "mixed")

    with open(os.path.join(tmpdir, "invalid.fasta"), "r") as fobj:
        header = fobj.readline()

    # test if trailing whitespace
    assert not header[:-1].endswith(" ")


def test_download():
    expected_md5sum = "f89b59a16cac6fea465d41fb98e16fac"
    url = (
        "https://gitlab.com/RKIBioinformaticsPipelines/president/"
        "-/raw/master/tests/fixtures/NC_045512.2_cds.fasta"
    )
    suffix = "_cds.fasta"
    path = writer.download(url, suffix)
    sequences = open(path, "r", encoding='utf-8').read().encode()
    md5sum = hashlib.md5(sequences).hexdigest()
    assert md5sum == expected_md5sum
