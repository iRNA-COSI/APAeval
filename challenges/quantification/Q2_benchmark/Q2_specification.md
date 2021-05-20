# Q2 benchmark specification

## Synopsis

Benchmark to test correlation between RNAseq-based poly(A) site quantification and poly(A)-sites quantification based on orthogonal 3'end seq data.

Input data:

1. Poly(A) sites quantification based on RNAseq data using the benchmarked tool
2. Poly(A) sites quantification based on orthogonal 3'end seq dataset

Based on the input data the following metrics are computed:
 
1. Correlation between RNAseq-based quantification and 3'end seq quantification.

The metrics should be computed for the following distance threshold between PAS identified by the tool in RNAseq dataset and PAS identified in orthogonal 3'end seq dataset:

- 50 nt

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

This BED file contains positions of unique poly(A) sites as well as an additional column with TPM values for each identified site quantified from RNAseq data by the benchmarked tool.

#### Format 2

This BED file contains positions of unique poly(A) sites as well as an additional column with TPM values for each identified site quantified from the orthogonal 3'end-seq dataset

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
 | `correlation_coefficient` | `float` | N/A | Correlation between RNAseq-based quantification and 3'end seq quantification; calculated for distance threshold of 50 nt |

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
