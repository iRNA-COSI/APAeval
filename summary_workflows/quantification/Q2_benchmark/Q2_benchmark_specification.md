# Q2 benchmark specification

## Synopsis

Benchmark to test correlation between RNAseq-based poly(A) site quantification and poly(A)-sites quantification based on orthogonal 3'end seq data.

Input data:

1. Poly(A) sites quantification based on RNAseq data using the benchmarked tool
2. Poly(A) sites quantification based on orthogonal 3'end-seq dataset

Based on the input data the following metrics are computed:
 
1. Correlation between RNAseq-based quantification and 3'end-seq quantification.

To compute the metrics, the poly(A) sites identified from the RNAseq data using the benchmarked tool and from the orthogonal dataset should be first mapped to each other:

- find multiple ground truth sites overlapping one predicted site
  - weights are added for the predicted site based on distance to the ground truth site
  - expression can be calculated by `weight_i*expression_i` for all prediction sites having the same ground truth site assigned
- find multiple predicted sites overlapping one ground truth site
  - weigths are added for the predicted sites based on distance to the ground truth site
  - expression can be calculated by `sum(weight_i*expression_i)` to get a single value matching the one ground truth site OR just summed without considering weight


## General info

* **Challenge:** Quantification
* **Identifier:** Q2

## Inputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | BED | [Specification][spec-bed] | [Link][in1] |
| 2 | BED | [Specification][spec-bed] | [Link][in2] |

### Additional info inputs

#### Format 1

This BED file contains single-nucleotide positions of unique poly(A) sites as well as TPM values for each identified site quantified from RNAseq data by the benchmarked tool.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - the TPM value - quantified PAS
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 2

This BED file contains single-nucleotide positions of unique poly(A) sites as well as TPM values for each identified site quantified from the orthogonal 3'end-seq dataset.  
Fields (the same as format 1):

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - the TPM value - quantified PAS
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
 | `correlation_coefficient` | `float` | N/A | Correlation between RNAseq-based quantification and 3'end-seq quantification |

## Metrics
 
 | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
 | --- | --- | --- | --- | --- | --- | --- |
 | 1 | Correlation | N/A | Output 1 | Read file, parse JSON and extract attribute `correlation_coefficient` | `vector` | N/A |
 
### Additional info metrics
 
 N/A

[//]: # (References)
 
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
