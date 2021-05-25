# I2 benchmark specification

## Synopsis

Benchmark to test sensitivity and False Discovery Rate (FDR) of RNAseq-based poly(A) site identification compared with ground truth (poly(A)-sites identified from orthogonal 3'end seq data).

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites identified based on orthogonal 3'end seq dataset

Based on the input data the following metrics are computed:

1. Sensitivity = (TP/(TP+FN))
2. FDR = (FP/(TP+FP))

TP - true positives - PAS identified by the tool and present in the orthogonal dataset  
FP - false positives - PAS identified by the tool and not present in the orthogonal dataset  
FN - false negatives - PAS not identified by the tool but present in the orthogonal dataset

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the orthogonal 3'end seq dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.  
The following distance thresholds should be used:

- 0 nt
- 10 nt
- 20 nt
- 30 nt
- 40 nt
- 50 nt
- 60 nt
- 70 nt
- 80 nt
- 90 nt
- 100 nt

## General info

* **Challenge:** Identification
* **Identifier:** I2

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
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 2

This BED file contains single-nucleotide position of poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

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
  | `sensitivity` | `vector` | N/A | A vector of length=11 containing sensitivity values of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `FDR` | `vector` | N/A | A vector of length=11 containing False Discovery Rate values of PAS identification compared with orthogonal dataset; FDR = (FP/(TP+FP)); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Sensitivity | N/A | Output 1 | Read file, parse JSON and extract attribute `sensitivity` that has type `vector` of `float`. The values should be plotted against a vector of distance cutoff values from 0 to 100 nt with 10 nt intervals. | `array` | N/A |
  | 2 | FDR | N/A | Output 1 | Read file, parse JSON and extract attribute `FDR` that has type `vector` of `float`. The values should be plotted against a vector of distance cutoff values from 0 to 100 nt with 10 nt intervals. | `array` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
  [in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
  [spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
  