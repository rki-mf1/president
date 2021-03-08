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
# authors:
# RKI MF1;  Martin Hoelzer with great initial help of @phiweger (UKL Leipzig)
# HPI;      Fabio Malcher Miranda, Sven Giese, Alice Wittig
import argparse
import os
from shutil import which
import sys
from datetime import datetime
import pandas as pd
import logging
from president import alignment, __version__, statistics, writer, sequence

logger = logging.getLogger(__name__)


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
            threads=4, store_alignment=False, verbose=False):  # pragma: no cover
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
    verbose: bool,
            if True, print logging information to screen
    Returns
    -------
    datframe,
        result metrics
    """
    # handle path / prefix input
    out_dir = os.path.abspath(path_out)
    preloginfo = ""
    if not os.path.exists(out_dir):
        preloginfo = "Creating output directory..."
        os.makedirs(out_dir, exist_ok=True)

    # create logger
    logger = logging.getLogger('president')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.FileHandler(os.path.join(path_out, f"{prefix_out}president_logger.log"), "w")
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(ch)

    if not verbose:
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(sh)

    logger.info("Starting Time: {}".format(datetime.now().strftime("%H:%M:%S")))
    logger.info("Using president version: {}".format(__version__))
    logger.info(f"Writing files to: {path_out}")
    logger.info(f"Using the prefix: {prefix_out}_* to store results.")

    logger.info(preloginfo)
    logger.info(f"Query: {query_in_raw}")
    logger.info(f"Reference: {reference_in}")
    # Files exist?
    if not os.path.isfile(reference_in):
        logger.error(f"Reference {reference_in} does not exist.")
        raise FileNotFoundError(reference_in)

    # make sure input is iterable
    if not isinstance(query_in_raw, list):
        query_in_raw = [query_in_raw]

    collect_dfs = []
    logger.info("Starting president processing")
    # preprocess fasta files
    query_tmp, query_source = sequence.preprocess(query_in_raw, "_query.fasta")
    reference_tmp, _ = sequence.preprocess(reference_in, "_reference.fasta")

    prefix = prefix_out
    assert os.path.isfile(query_tmp)

    # check reference fasta
    _ = statistics.count_sequences(reference_tmp)
    query_valid = statistics.count_sequences(query_tmp, "query")

    # check query data
    if query_valid:
        summary_stats_query = statistics.summarize_query(query_tmp)
        statistics.qc_check(reference_tmp, summary_stats_query, id_threshold=id_threshold,
                            n_threshold=n_threshold)

        # perform initial sequence check
        query_valid_tmp, evaluation, invalid_ids = \
            statistics.split_valid_sequences(query_tmp, summary_stats_query)

        logger.info(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")
    else:
        # make sure mock dataframe looks like regular one
        evaluation = "all_invalid"
        summary_stats_query = pd.DataFrame(writer.init_metrics(0))
        summary_stats_query = \
            summary_stats_query.assign(**{'qc_all_valid': [], 'qc_valid_length': [],
                                          'qc_valid_nucleotides': [], 'qc_valid_number_n': [],
                                          'qc_valid_pass_nthreshold': []})

    # align sequences if more than 1 sequence passes the initial qc
    if evaluation != "all_invalid":
        alignment_file = alignment.pblat(threads, reference_tmp, query_valid_tmp)
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
        metrics["file_in_query"] = query_source[metrics["query_index"].astype(int)]
    else:
        metrics["file_in_query"] = 'NaN'

    metrics["file_in_ref"] = os.path.basename(reference_in)

    # pretify output
    metrics = metrics[metrics.columns.sort_values()]
    # putting the columns with PSL prefix at the end
    PSL_columns = metrics.filter(regex="^PSL_").columns.to_list()
    metrics = metrics[writer.COL_ORDER + PSL_columns]
    metrics.to_csv(os.path.join(out_dir, f"{prefix}report.tsv"), index=False, sep='\t')

    # remove temporary files
    if evaluation != "all_invalid":
        os.remove(alignment_file)

    os.remove(query_tmp)
    os.remove(reference_tmp)
    logger.info(metrics)
    logger.info(metrics.shape)
    collect_dfs.append(metrics)

    # if there are more input files to iterate from, concat results
    metrics_all = pd.concat(collect_dfs)

    logger.info("president finished from command:")
    logger.info(f"Call:'president -r {reference_in} -q {query_in_raw} -p {path_out} "
                f"-f {prefix_out} -x {id_threshold} -n {n_threshold} -a {store_alignment} "
                f"-t {threads}'")
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
    parser.add_argument('-e', '--quite', dest="verbose", action="store_true", default=False,
                        help="Print log messages also to the screen (False)", required=False)
    args = parser.parse_args()

    aligner(args.reference, args.query, args.path, args.prefix, args.id_threshold,
            args.n_threshold, args.threads, args.store_alignment, args.verbose)


if __name__ == "__main__":
    main()
