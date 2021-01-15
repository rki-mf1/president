![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)
![](https://img.shields.io/badge/python-3.8-orange)

[![](https://img.shields.io/badge/ANI-definition-violet.svg)](https://pubmed.ncbi.nlm.nih.gov/17220447/)

# PRESIDENT: PaiRwisE Sequence IDENtiTy
Calculate pairwise nucleotide identity with respect to a reference sequence.

Given a reference and a query sequence (which can be fragmented), calculate
pairwise nucleotide identity with respect to the reference sequence relative 
to the entire length of the reference. Only informative nucleotides (A,T,G,C) 
are considered identical to each other.

#### Requirements:
```
conda create -y -n president -c bioconda python=3.8 pblat=2.5 screed pandas && conda activate president
```

#### Usage:
```
python pairwise_nucleotide_identity.py --query NC_045512.2.20mis.fasta --reference NC_045512.2.fasta -x 3000 -p 8
```

#### Output:
The script provides __three different output metrics__ to assess the quality of the query sequence in comparison to the reference sequence. The rationale behind this is, that ambiguous bases (Ns) can impact sequence identity differently depending on how they are counted. 

```
Running pblat ...
a) 0.9987 nucleotide identity
b) 0.9994 nucleotide identity (excluding non-ACTG)
c) 20 non-ACTG symbols out of 29903
```

a) Percentage identity based on hamming distance and including gaps/ Ns as mismatches 

b) Percentage identity only based on A, T, C, G characters `(ACTG matches / (len_query - nonACTGchars_query))`

c) Count all characters in the query sequence that are not A, T, C, G

#### Notes:
- nextstrain uses a quality threshold of < 3000 non-canonical nucleotides
- Ns in the query are treated as mismatches, uncomment the corresponding line directly in the code to ignore Ns

#### ANI definition:
- https://pubmed.ncbi.nlm.nih.gov/17220447/

