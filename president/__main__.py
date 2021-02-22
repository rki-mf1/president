"""President's main module.

Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage:
-----
president --query tiny_test.masked_consensus.fasta \
    --reference NC_045512.2.fasta -x 0.9 -p output/report -t 4

Notes:
-----
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line
(below) to ignore Ns

Requirements:

- conda install -y -c bioconda python=3.8 pblat screed pandas

ANI definition:

- https://pubmed.ncbi.nlm.nih.gov/17220447/

Spec:
-----
> Sequence identity to a reference genome sequence (e.g. NC_045512.2) should
be calculated at the nucleotide level relative to the entire length of the
reference. Only informative nucleotides (A,T,G,C) are considered identical
to each other.
"""
import argparse
import os
# authors:
# RKI MF1;  Martin Hoelzer with great initial help of @phiweger (UKL Leipzig)
# HPI;      Fabio Malcher Miranda, Sven Giese, Alice Wittig
from shutil import which
import pandas as pd

from president import alignment, __version__, statistics, writer, sequence


def is_available(name="pblat"):
    """
    Check whether `name` is on PATH and marked as executable.

    Parameters
    ----------
    name : str
        tool name executable (default: pblat).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # tool installed
    if not which(name):
        raise ValueError(f'{name} not on PATH or marked as executable.')


def aligner(reference_in, query_in_raw, path_out, prefix_out="",
            id_threshold=0.9, n_threshold=0.05,
            threads=4, store_alignment=False):  # pragma: no cover
    """
    Align query to the reference and extract qc metrics.

    Parameters
    ----------
    reference_in : str
        reference FASTA location
    query_in_raw : str / list
        query FASTA location(s).
    path_out : str
        Path to be used to store the results.
    prefix_out : str
        File prefix to be used for each output file.
    id_threshold : float, optional
        Identity threshold after aligment that must be achieved. The default is 0.9.
    n_threshold: float, optional
        Percentage of allowed Ns, sequences with hier N% will be rejected
    threads : int, optional
        Number of threads to use. The default is 4.
    store_alignment : bool, optional
        Should we store the results of the alignements in the qc metrics table? default is False
    Returns
    -------
    datframe,
        result metrics
    """
    # Files exist?
    assert os.path.isfile(reference_in)
    print(query_in_raw)
    # make sure input is iterable
    if not isinstance(query_in_raw, list):
        query_in_raw = [query_in_raw]

    collect_dfs = []
    print("##################### Running President ##########################")
    # preprocess fasta files
    query_tmp, query_source = sequence.preprocess(query_in_raw, "_query.fasta")
    reference_tmp, _ = sequence.preprocess(reference_in, "_reference.fasta")

    prefix = prefix_out
    print(f"Running file: {query_tmp}")
    assert os.path.isfile(query_tmp)

    # handle path / prefix input
    out_dir = os.path.abspath(path_out)
    if not os.path.exists(out_dir):
        print("Creating output directory...")
        os.makedirs(out_dir, exist_ok=True)

    print(f"Writing files to: {out_dir}")
    print(f"Using the prefix: {prefix}* to store results.")

    # check reference fasta
    _ = statistics.count_sequences(reference_tmp)
    query_valid = statistics.count_sequences(query_tmp, "query")

    # check query data
    if query_valid:
        summary_stats_query = statistics.summarize_query(query_tmp)
        statistics.qc_check(reference_tmp, summary_stats_query, id_threshold=id_threshold,
                            n_threshold=n_threshold)

        # perform initial sequence check
        query_tmp, evaluation, invalid_ids = \
            statistics.split_valid_sequences(query_tmp, summary_stats_query)

        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")
    else:
        # make sure mock dataframe looks like regular one
        evaluation = "all_invalid"
        summary_stats_query = pd.DataFrame(writer.init_metrics(0))
        summary_stats_query = \
            summary_stats_query.assign(**{'qc_all_valid': [], 'qc_valid_length': [],
                                          'qc_valid_nucleotides': [], 'qc_valid_number_n': []})

    # align sequences if more than 1 sequence passes the initial qc
    if evaluation != "all_invalid":
        alignment_file = alignment.pblat(threads, reference_tmp, query_tmp, verbose=1)
        # parse statistics from file
        metrics = statistics.nucleotide_identity(alignment_file, summary_stats_query,
                                                 id_threshold, store_alignment)
        metrics["qc_is_empty_query"] = False
    else:
        # if no sequences are there to be aligned, create a pseudooutput that looks
        # exactly as the aligned output
        metrics = writer.init_metrics(1, extend_cols=True, metrics_df=summary_stats_query)
        metrics["qc_is_empty_query"] = True
    # store sequences
    writer.write_sequences(query_tmp, metrics, os.path.join(out_dir, f"{prefix}"), evaluation)

    # store reference data
    if len(query_source) > 0:
        metrics["file_in_query"] = query_source
    else:
        metrics["file_in_query"] = 'NaN'
    metrics["file_in_ref"] = os.path.basename(reference_in)
    metrics = metrics[metrics.columns.sort_values()]
    # putting the columns with PSL prefix at the end
    PSL_columns = metrics.filter(regex="^PSL_").columns.to_list()
    other_columns = [c for c in metrics.columns if c not in PSL_columns]
    metrics = metrics[other_columns + PSL_columns]

    metrics.to_csv(os.path.join(out_dir, f"{prefix}report.tsv"), index=False, sep='\t')

    # remove temporary files
    if evaluation != "all_invalid":
        os.remove(alignment_file)

    os.remove(query_tmp)
    os.remove(reference_tmp)
    print(metrics)
    print(metrics.shape)
    collect_dfs.append(metrics)

    # if there are more input files to iterate from, concat results
    metrics_all = pd.concat(collect_dfs)
    return metrics_all


def main():  # pragma: no cover
    """
    Presidents main function for sequence alignment.

    Returns
    -------
    None.

    """
    parser = argparse.ArgumentParser(description='Calculate pairwise nucleotide identity.')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome.')
    parser.add_argument('-q', '--query', required=True, help='Query genome(s).', nargs="+")
    parser.add_argument('-x', '--id_threshold', type=float, default=0.9,
                        help='ACGT nucleotide identity threshold after alignment (percentage). '
                             'A query sequence is reported as valid based on this threshold '
                             '(def: 0.9)')
    parser.add_argument('-n', '--n_threshold', type=float, default=0.05,
                        help='A query sequence is reported as valid, if the percentage of Ns '
                             'is smaller or equal the threshold (def: 0.05)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use for pblat.')
    parser.add_argument('-p', '--path', required=True,
                        help='Path to be used to store results and FASTA files.')
    parser.add_argument('-f', '--prefix', required=False, default="",
                        help='Prefix to be used t store results in the path')
    parser.add_argument('-a', '--store_alignment',
                        required=False, action="store_true",
                        help='add query alignment columns (PSL format)')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    aligner(args.reference, args.query, args.path, args.prefix, args.id_threshold,
            args.threads, args.store_alignment)


if __name__ == "__main__":
    main()
