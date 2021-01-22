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
conda create -y -n president -c bioconda president
conda activate president
```

Note that `pblat` is a dependency and only runs on Linux. Alternatively, president can be installed with pip in an environment where pblat is in the PATH:

```
pip install president
```

#### Usage:
pypresident installs the package and the pairwise alignment can be run with the following console call:

```
president --query query.fasta --reference reference.fasta -x identity_threshold -t threads -p /path/to/output/prefix
```

To run an example, download a [query](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/NC_045512.2.20mis.fasta) FASTA and
a [reference](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/NC_045512.2.fasta) FASTA from GitLab.

Run the alignment with the following command and identity of ACGT bases of 93%. Note that multiple fasta sequences are allowed to be present in the query but **not** in the reference FASTA.

```
president -q NC_045512.2.20mis.fasta -r NC_045512.2.fasta -x 0.93 -t 4 -p output/test
```


#### Output:
The script provides:

* a tab-separated file with the below listed columns
* a FASTA file with _valid_ sequences
* a FASTA file with _invalid_ sequences

The separation between the _valid_ and _invalid_ bin is mainly based on the defined identity threshold (`-x`, default: 0.93) and further sanity checks (non-IUPAC characters, amount of 'N's and query length that cause sequence identity to drop below `-x`).

An example from the call described above can be found [online](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/report.csv).

The transposed version of the output is shown below.

| variable                                          | value                                                                                           |
|:--------------------------------------------------|:------------------------------------------------------------------------------------------------|
| ID                                                | NC_045512.2                                                                                     |
| Valid                                             | True                                                                                            |
| ACGT Nucleotide identity                          | 0.9987                                                                                          |
| ACGT Nucleotide identity (ignoring Ns)            | 0.9994                                                                                          |
| ACGT Nucleotide identity (ignoring non-ACGTNs)    | 1.0                                                                                             | 
| Ambiguous Bases                                   | 20.0                                                                                            |
| Query Length                                      | 29903                                                                                           |
| Query #ACGT                                       | 29883                                                                                           |
| Query #IUPAC-ACGT                                 | 20.0                                                                                            |
| Query #non-IUPAC                                  | 0.0                                                                                             |
| aligned                                           | True                                                                                            |
| passed_initial_qc                                 | True                                                                                            |
| Date                                              | 2021-01-20                                                                                      |
| reference_length                                  | 29903                                                                                           |
| reference                                         | NC_045512.2.fasta                                                                               |
| query                                             | NC_045512.2.20mis.fasta                                                                         |


##### Definitions:

1) ID - the fasta sequence name in the query
2) Valid - True if 'ACGT Nucleotide identity' is greater than the -x parameter
3) ACGT Nucleotide identity - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths))
4) ACGT Nucleotide identity (ignoring Ns) - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths) - #Ns in query)
5) ACGT Nucleotide identity (ignoring non-ACGTNs) - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths) - #non-ACGTNs in query)
6) Ambiguous Bases - Number of N nucleotides in the query.
7) Query Length - Sequence length of the query.
8) Query #ACGT - Number of ACGT in the query.
9) Query #IUPAC-ACGT - Number of IUPAC characters that are not ACGT in the query.
10) Query #non-IUPAC - Number of non-IUPAC characters in the query.
11) aligned - True if the sequence got aligned with `pblat`
12) passed_initial_qc - True if the sequence is long enough / has enough ACGT nucleotides (instead of Ns) to reache the identiy threshold
13) Date - the yyyy-mm-dd of the execution of the script
14) reference_length -  Sequence length of the reference
15) reference - the reference file name
16) query - the query file name

__Note__: max(sequence_lengths) is equal to max(length_query, length_reference).


##### Notes:
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides


#### ANI definition:
- https: // pubmed.ncbi.nlm.nih.gov/17220447/
