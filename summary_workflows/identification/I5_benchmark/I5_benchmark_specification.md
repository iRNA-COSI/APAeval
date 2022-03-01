# I5 benchmark specification

## Synopsis

Benchmark to compare poly(A) site identification from **simulated** RNAseq data with ground truth, i.e. the poly(A) sites of transcripts used to simulate RNAseq data, for several different distance thresholds.

Input data:

1. Poly(A) sites identified based on simulated RNAseq data using the benchmarked tool
2. Poly(A) sites of transcripts used to simulate RNAseq reads

Based on the input data the following metrics are computed:

1. Precision = (TP/(TP+FP))

TP - true positives - PAS identified by the tool and present within X nucleotides of PAS of transcripts from which the reads were simulated  
FP - false positives - PAS identified by the tool and not present within X nucleotides of PAS of transcripts from which the reads were simulated

The metrics should be computed for different distance thresholds between PAS identified by the tool from simulated RNAseq dataset and PAS of transcripts used to simulate RNAseq reads, i.e. the PAS identified by the tool should be within X nucleotides from the PAS of transcripts from which the reads were simulated for the prediction to be considered true.  
Distance thresholds should be between 0 nt and 200 nt with 20 nt increments.

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

This BED file contains single-nucleotide position of poly(A) sites identified by the benchmarked method.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 2

This BED file contains PAS of transcripts that were used to **simulate** RNA-seq reads.  
Fields are the same as in Format 1.

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

