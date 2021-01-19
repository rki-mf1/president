![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)
![](https://img.shields.io/badge/python-3.8-orange)

[![](https://img.shields.io/badge/ANI-definition-violet.svg)](https://pubmed.ncbi.nlm.nih.gov/17220447/)

[![pipeline status](https://gitlab.com/RKIBioinformaticsPipelines/president/badges/master/pipeline.svg)](https://gitlab.com/RKIBioinformaticsPipelines/president/-/commits/master)
[![coverage report](https://gitlab.com/RKIBioinformaticsPipelines/president/badges/master/coverage.svg)](https://gitlab.com/RKIBioinformaticsPipelines/president/-/commits/master)

#### PRESIDENT: PaiRwisE Sequence IDENtiTy
Calculate pairwise nucleotide identity with respect to a reference sequence.

Given a reference and a query sequence(which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence relative
to the entire length of the reference. Only informative nucleotides (A, T, G, C)
are considered identical to each other.

#### Requirements:
To get president running follow the steps below. Note that pbat only runs on linux.

```
conda create -y -n president -c bioconda python=3.8 pblat
conda activate president
pip install pypresident
```

#### Usage:
pypresident installs the package and the pairwise alignment can be run with
the following console call:

```
president --query query.fasta --reference reference.fasta -x identity_threshold -p threads -o output.tsv
```

To run an example, download a [query](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/NC_045512.2.20mis.fasta) FASTA and
a [reference](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/NC_045512.2.fasta) FASTA from GitLab.

Run the alignment with the following command and require and identity of ACGT bases of 93%. Note
that multiple fasta sequences are allowed to be present in the query but **not** in the reference.


```
president --query NC_045512.2.20mis.fasta --reference NC_045512.2.fasta -x 0.93 -p 8 -o report.tsv
```


#### Output:
The script provides a tab-separated file with the following columns. An example from the
call described above can be found [online](https://gitlab.com/RKIBioinformaticsPipelines/president/-/blob/master/examples/report.csv).


The transposed version of the output is shown below.

| variable           | value                                                                                           |
|:-------------------|:------------------------------------------------------------------------------------------------|
| ID                 | NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome |
| Valid              | True                                                                                            |
| Identity           | 0.9987                                                                                          |
| Ambiguous Identity | 0.9994                                                                                          |
| Ambiguous Bases    | 20.0                                                                                            |
| Query Length       | 29903                                                                                           |
| passed_initial_qc  | True                                                                                            |
| aligned            | True                                                                                            |
| reference_length   | 29903                                                                                           |


##### Definitions:

1) ID - the fasta sequence name in the query
2) Valid - True if 'Identity' is greater than the -x parameter
3) Identity - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths))
4) Ambiguous Identity - Percentage of ACGT matches to the reference (# matches / max(sequence_lengths) - #Ns in query)
5) Ambiguous Bases - Number of N nucleotides in the query.
6) Query Length - Sequence length of the query.
7) passed_initial_qc - True if the sequence is long enough / has enough ACGT nucleotides (instead of Ns) to reache the identiy threshold
8) aligned - True if the sequence got aligned with pblat
9) reference_length -  Sequence length of the reference.

Note: max(sequence_lengths) is equal to max(length_query, length_reference).


##### Notes:
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides


#### ANI definition:
- https: // pubmed.ncbi.nlm.nih.gov/17220447/
