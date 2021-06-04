# I3 benchmark specification

## Synopsis

Benchmark to compare output of tool's RNAseq based site identification to an established atlas of polyA-sites (PolyASite, PolyA_DB). The goal is to calculate precision with respect to data from existing databases for several different distance thresholds.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites identified in PolyASite database
3. Poly(A) sites identified in PolyA_DB database

Based on the input data the following metrics are computed:

1. Precision = (TP/(TP+FP))

TP - true positives - PAS identified by the tool and present within X nucleotides of PAS identified in the database  
FP - false positives  - PAS identified by the tool and not present within X nucleotides of PAS identified in the database

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified in the database i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified in the database for the prediction to be considered true.  
The following distance thresholds should be used:

- 0 nt
- 20 nt
- 40 nt
- 60 nt
- 80 nt
- 100 nt
- 120 nt
- 140 nt
- 160 nt
- 180 nt
- 200 nt


## General info

* **Challenge:** Identification
* **Identifier:** I3

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | BED | [Specification][spec-bed-polyAsite] | [Link][in2] |
  | 3 | custom TSV | [Specification][spec-custom-polyAdb] | [Link][in3] |

### Additional info inputs
  
#### Format 1

This BED file contains PAS identified by the benchmarked tool from RNAseq data. Poly(A) sites should be single nucleotide.

#### Format 2

This BED file contains PAS from PolyASite database

#### Format 3

This TSV file contains PAS from PolyA_DB database


## Outputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | JSON | [Specification][spec-json] | [Link][out1] |
  | 2 | JSON | [Specification][spec-json] | [Link][out2] |
  
### Additional info outputs

#### Format 1
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with PolyASite database; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float` |


#### Format 2
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float` |

A vector of length=11 containing precision values of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float`
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Precision of PAS identification compared with PolyASite database | N/A | Output 1 | Read file, parse JSON and extract attribute `Precision` that has type `vector` of `float`. The values should be plotted against a vector of distance cutoff values from 0 to 200 nt with 20 nt intervals. | `matrix` | N/A |
  | 2 | Precision of PAS identification compared with PolyA_DB database | N/A | Output 2 | Read file, parse JSON and extract attribute `Precision` that has type `vector` of `float`. The values should be plotted against a vector of distance cutoff values from 0 to 200 nt with 20 nt intervals. | `matrix` | N/A |

### Additional info metrics
  
  N/A

[//]: # (References)
  
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[in3]: ./example_files/input3.tsv

[out1]: ./example_files/output1.json
[out2]: ./example_files/output2.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[spec-bed-polyAsite]: <https://polyasite.unibas.ch/atlas>
[spec-custom-polyAdb]: <https://exon.apps.wistar.org/polya_db/v3/download/3.2/readme.txt>
