"""President's main module.

Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence.

Usage:
-----
president --query tiny_test.masked_consensus.fasta \
    --reference NC_045512.2.fasta -x 0.93 -p 8 -o report.tsv

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


def aligner(reference_in, query_in_raw, prefix_in, id_threshold=0.93,
            threads=4):  # pragma: no cover
    """
    Align query to the reference and extract qc metrics.

    Parameters
    ----------
    reference_in : str
        reference FASTA location
    query_in_raw : str
        query FASTA location.
    prefix_in : str
        Prefix where to store the results.
    id_threshold : float, optional
        Identity threshold after aligment that must be achieved. The default is 0.93.
    threads : int, optional
        Number of threads to use. The default is 4.

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

    for qi, query_in in enumerate(query_in_raw):
        print("##################### Running President ##########################")
        prefix = prefix_in
        print(f"Running file: {query_in}")
        assert os.path.isfile(query_in)

        # handle path / prefix input
        out_dir = os.path.dirname(os.path.abspath(prefix))
        # dir is not created when prefix ends with /
        if prefix.endswith('/'):
            out_dir = prefix

        if prefix != "" and not prefix.endswith('/'):
            prefix += "_"

        file_prefix = os.path.basename(prefix)
        if not os.path.exists(out_dir):
            print("Creating output directory...")
            os.makedirs(out_dir, exist_ok=True)

        print(f"Writing files to: {out_dir}")
        print(f"Using the prefix: {file_prefix}_* to store results.")

    # check reference fasta
    statistics.count_reference_sequences(reference_tmp)
    # check query data
    summary_stats_query = statistics.summarize_query(query_in)
    statistics.qc_check(reference_tmp, summary_stats_query, id_threshold=id_threshold)

    # perform initial sequence check
    query_tmp, evaluation, invalid_ids = \
        statistics.split_valid_sequences(query_tmp, summary_stats_query)

    # if none of the sequences pass the qc filter, exit.
    # else just perform the alignment witht he seqs passing qc
#    if evaluation == "all_invalid":
#        print("None of the sequences can pass the identity threshold. No alignment done.")
#        print("Exiting president alignment.")
#        sys.exit()
#    else:
#        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")

        # perform initial sequence check
        query_tmp, evaluation, invalid_ids = \
            statistics.split_valid_sequences(query_tmp, reference_tmp, id_threshold=id_threshold)
        # if none of the sequences pass the qc filter, exit.
        # else just perform the alignment witht he seqs passing qc
    #    if evaluation == "all_invalid":
    #        print("None of the sequences can pass the identity threshold. No alignment done.")
    #        print("Exiting president alignment.")
    #        sys.exit()
    #    else:
    #        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")

        print(f"Performing alignment with valid sequences (excluding {len(invalid_ids)}).")

        # parse statistics from file
        metrics = statistics.nucleotide_identity(alignment_file, summary_stats_query, id_threshold)

    # store sequences
    writer.write_sequences(query_in, metrics, prefix, evaluation)

        print(metrics)
                           mode="a", header=False)
            metrics.to_csv(os.path.join(out_dir, f"{file_prefix}report.tsv"), index=False, sep='\t',
        else:
            metrics.to_csv(os.path.join(out_dir, f"{file_prefix}report.tsv"), index=False, sep='\t')
        if qi == 0:
        os.remove(query_tmp)
        os.remove(reference_tmp)
            os.remove(alignment_file)
        if evaluation != "all_invalid":
        # remove temporary files
    # store reference data
    metrics["file_in_query"] = os.path.basename(query_in)
    metrics["file_in_ref"] = os.path.basename(reference_in)
    return metrics


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
    parser.add_argument('-x', '--id_threshold', type=float, default=0.93,
                        help='ACGT nucleotide identity threshold after alignment (percentage). '
                             'A query sequence is reported as valid based on this threshold '
                             '(def: 0.93)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use for pblat.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Prefix to be used to store results and FASTA files.')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    aligner(args.reference, args.query, args.prefix, args.id_threshold, args.threads)


if __name__ == "__main__":
    main()
