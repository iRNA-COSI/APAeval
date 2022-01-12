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
  
The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the orthogonal 3'end seq dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.  
Distance thresholds should be between 0 nt and 100 nt with 10 nt increments.

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

This BED file contains genomic positions of unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from RNAseq data by the benchmarked tool.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - the TPM value - quantified PAS
- **strand** - defines the strand; either "." (=no strand) or "+" or "-"

#### Format 2

This BED file contains genomic positions of unique unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from the orthogonal 3'end-seq dataset.
Fields are the same as in format 1.

## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

1. **bar plot** visualizing **Correlation** of poly(A) site quantification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - benchmarked tool  
**Y axis** - correlation

Ranking: The best performing tool is the one with the highest correlation value.

Optional plots:

2.  **2D line plot** visualising **correlation** as a function of distance threshold.

**X axis** - correlation(d)  
**Y axis** - distance threshold _d_ 

Ranking: For each tool, the area under the curve (AUC) is calculated. The best performing tool is the one with the highest AUC.

Note: 2D line plot is not supported in OpenEBench yet. If it's not implemented, the data should be visualised outside of OpenEBench.

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
| `correlation_coefficient` | `vector` | N/A | A vector of length=11 containing values of correlation between RNAseq-based quantification and 3'end-seq quantification; calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |

## Metrics
 
| # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Correlation | N/A | Output 1 | Read file, parse JSON and extract attribute `correlation_coefficient` that has type `vector` of `float` | `array` | N/A |
 
### Additional info metrics
 
 N/A

[//]: # (References)
 
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
