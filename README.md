![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)
![](https://img.shields.io/badge/python-3.8-orange)

[![](https://img.shields.io/badge/ANI-definition-violet.svg)](https://pubmed.ncbi.nlm.nih.gov/17220447/)

[![pipeline status](https://gitlab.com/RKIBioinformaticsPipelines/president/badges/master/pipeline.svg)](https://gitlab.com/RKIBioinformaticsPipelines/president/-/commits/master)
[![coverage report](https://gitlab.com/RKIBioinformaticsPipelines/president/badges/master/coverage.svg)](https://gitlab.com/RKIBioinformaticsPipelines/president/-/commits/master)

#### PRESIDENT: PaiRwisE Sequence IDENtiTy
Calculate pairwise nucleotide identity with respect to a reference sequence.

Given a reference and a query sequence, calculate pairwise nucleotide identity with respect to the reference sequence relative to the entire length of the reference. In the main metric, only informative nucleotides (A, T, G, C) are considered identical to each other. The tool also provides some further metrics (e.g. regarding ambiguous 'N's) and splits the input FASTA into _valid_ and _failed_ FASTA files for further processing. 

#### Installation:
To install president with conda, run the commands below:

```
conda create -y -n president -c bioconda -c conda-forge president
conda activate president
```

Note that `pblat` is a dependency and only runs on Linux. Alternatively, president can be installed with pip in an environment where pblat is in the PATH:

```
pip install president
```

#### Usage:
pypresident installs the package and the pairwise alignment can be run with the following console call:

```
president --query query.fasta --reference reference.fasta -x identity_threshold -n n_threshold -t threads -p /path/to/output/ -f prefix
```

To run an example, download a [query](https://gitlab.com/RKIBioinformaticsPipelines/president/-/raw/master/examples/NC_045512.2.20mis.fasta) FASTA and
a [reference](https://gitlab.com/RKIBioinformaticsPipelines/president/-/raw/master/examples/NC_045512.2.fasta) FASTA from GitLab.

Run the alignment with the following command and identity of ACGT bases of 93% and maximal 5% Ns in the query. Note that multiple fasta sequences are allowed to be present in the query but **not** in the reference FASTA.

```
president -q NC_045512.2.20mis.fasta -r NC_045512.2.fasta -x 0.9 -n 0.05 -t 4 -p output/ -f test_
```


#### Output:
The script provides:

* a tab-separated file with the below listed columns
* a FASTA file with _valid_ sequences
* a FASTA file with _invalid_ sequences

The separation between the _valid_ and _invalid_ bin is mainly based on the defined identity/n thresholds (`-x`, default: 0.9; `-n`, default: 0.05) and further sanity checks (non-IUPAC characters).
Note that for multiple query inputs the valid and invalid sequence headers are NOT associate with their filename. This meta information must be retrieved from the report.
An example from the call described above can be found [online](https://gitlab.com/RKIBioinformaticsPipelines/president/-/raw/master/examples/report.csv).

The transposed version of the output is shown below.

| variable                                          | value                                                                                           |
|:--------------------------------------------------|:------------------------------------------------------------------------------------------------|
| ACGT Nucleotide identity                          | 0.9987                                                                                          |
| ACGT Nucleotide identity (ignoring Ns)            | 0.9994                                                                                          |
| ACGT Nucleotide identity (ignoring non-ACGTNs)    | 0.9994                                                                                          | 
| Date                                              | 2021-01-20                                                                                      |
| Matches                                           | 29864                                                                                           |
| Mismatches                                        | 19                                                                                              |
| N_bases                                           | 20                                                                                              |  
| Ngap                                              | 20                                                                                              |
| acgt_bases                                        | 29883                                                                                           |
| file_in_query	                                    | NC_045512.2.20mis.fasta                                                                         |
| file_in_ref                                       | NC_045512.2.fasta                                                                               |
| iupac_bases                                       | 20                                                                                              |
| length_query                                      | 29903                                                                                           |          
| length_reference                                  | 29903                                                                                           |
| non_iupac_bases                                   | 0                                                                                               |
| qc_all_valid                                      | True                                                                                            |
| qc_is_empty_query                                 | False                                                                                           |
| qc_post_align_pass_threshold                      | True                                                                                            |
| qc_post_aligned                                   | True                                                                                            |
| qc_post_aligned_all_valid                         | True                                                                                            |
| qc_valid_length                                   | True                                                                                            |
| qc_valid_nucleotides                              | True                                                                                            |
| qc_valid_number_n                                 | True                                                                                            |
| query_description                                 |                                                                                                 |
| query_index                                       | 0                                                                                               |
| query_name                                        | NC_045512.2 Severe acute [...]                                                                  |
| reference_name                                    | NC_045512.2 Severe acute [...]                                                                  |


##### Definitions:

1) ACGT Nucleotide identity - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths))
2) ACGT Nucleotide identity (ignoring Ns) - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths) - #Ns in query)
3) ACGT Nucleotide identity (ignoring non-ACGTNs) - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths) - #non-ACGTNs in query)
4) Date - the yyyy-mm-dd of the execution of the script
5) Matches - the number of matches in the alignment
6) Mismatches - the number of mismatches in the alignment
7) N_bases - the number of N bases
8) LongestNGap - Lenght of the longest N gap
9) acgt_bases - number of ACGT bases
10) file_in_query - the input query file name
11) file_in_ref - the input reference file name (note that for multiple input queries a tempororay file name is reported)
12) non_iupac_bases - Number of non-IUPAC nucleotides in the query
13) qc_all_valid - True if all checks below are True
14) qc_is_empty_query - True if input query file is not empty
15) qc_post_align_pass_threshold - True if 'ACGT Nucleotide identity' is greater than the `-x` parameter
16) qc_post_aligned - Set to True if the sequence was aligned (can be False either because of initial checks or failed alignments)
17) qc_post_aligned_all_valid - True if all checks before alignment are True, alignment is only performed if this value is True
18) qc_valid_length - True if the query is actually able to achieve the `-x` score even if the query is shorter/ longer
19) qc_valid_nucleotides - True if only valid IUPAC characters are in the query
20) qc_valid_pass_nthreshold - True if `-n` percentage of Ns in the query is not exceeded
21) query_description - query description, if available
22) query_index - the position of the sequence in the query fasta input file (for multiple files the counter is not resetted)
23) query_name - FASTA ID of the query sequence + file origin, format <header>:<filename>
24) reference_name - FASTA ID of the reference sequence

__Note__: max(sequence_lengths) is equal to max(length_query, length_reference).


##### Notes:
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides


#### ANI definition:
- https://pubmed.ncbi.nlm.nih.gov/17220447/
