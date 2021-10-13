# I5 benchmark specification

## Synopsis

Benchmark to compare RNAseq-based site identification with ground truth, i.e. the poly(A) sites of transcripts used to **simulate** RNAseq data, for several different distance thresholds.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites of transcripts used to simulate RNAseq reads

Based on the input data the following metrics are computed:

1. Precision = (TP/(TP+FP))

TP - true positives - PAS identified by the tool and present within X nucleotides of PAS of transcripts from which the reads were simulated
FP - false positives  - PAS identified by the tool and not present within X nucleotides of PAS of transcripts from which the reads were simulated

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified in the database i.e. the PAS identified by the tool should be within X nucleotides from the PAS of transcripts from which the reads were simulated for the prediction to be considered true.

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
* **Identifier:** I5

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | BED | [Specification][spec-bed] | [Link][in2] |

### Additional info inputs
  
#### Format 1

This BED file contains PAS identified by the benchmarked tool from RNAseq data. Poly(A) sites should be single nucleotide.

#### Format 2

This BED file contains PAS of transcripts that were used to **simulate** RNA-seq reads.

## Outputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | JSON | [Specification][spec-json] | [Link][out1] |
  
### Additional info outputs

#### Format 1
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with the poly(A) sites of transcripts from which the RNAseq reads were simulated; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float` |

  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Precision of PAS identification compared with ground truth, i.e. poly(A) sites of transcripts used to simulate the RNAseq data | N/A | Output 1 | Read file, parse JSON and extract attribute `Precision` that has type `vector` of `float`. The values should be plotted against a vector of distance cutoff values from 0 to 200 nt with 20 nt intervals. | `array` | N/A |

### Additional info metrics
  
  N/A

[//]: # (References)
  
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed

[out1]: ./example_files/output1.json

[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>

