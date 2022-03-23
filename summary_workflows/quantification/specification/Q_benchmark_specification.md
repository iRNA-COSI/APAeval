# Quantification benchmarks specification

## Synopsis

Benchmarks to assess the performance of RNA-Seq-based poly(A) site quantification compared to poly(A)-sites quantification based on orthogonal 3'end seq data.

### Input data:

Comparison of PAS predicted from RNA-Seq data with 3'end sequencing data

1. Poly(A) sites quantification based on RNAseq data using the benchmarked tool  
2. Ground truth: Poly(A) sites quantification based on orthogonal 3'end-seq dataset

Comparison of PAS predicted from simulated RNA-Seq data with dataset used for simulation
  
1. Poly(A) sites quantification based on simulated RNAseq data using the benchmarked tool  
2. Ground truth: Poly(A) sites quantification based on abundances of transcripts used to simulate the RNAseq data

### Metrics

Based on the input data the following metrics are computed:
 
1. Correlation between RNAseq-based quantification and 3'end-seq quantification.

To compute the metrics, the poly(A) sites identified from the RNAseq data by the benchmarked tool and from the ground-truth dataset should be first mapped to each other:

- find multiple ground truth sites overlapping one predicted site
  - weights are added for the predicted site based on distance to the ground truth site
  - expression can be calculated by `weight_i*expression_i` for all prediction sites having the same ground truth site assigned
- find multiple predicted sites overlapping one ground truth site
  - weigths are added for the predicted sites based on distance to the ground truth site
  - expression can be calculated by `sum(weight_i*expression_i)` to get a single value matching the one ground truth site OR just summed without considering weight
  
The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the ground truth dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.

### OpenEBench challenges

The metrics are visualised using 2D scatter plots and barplots, as described in _Plots_ section, which are then used for ranking the participating tools.
A plot for any given input dataset constitutes a benchmarking challenge as understood in the OpenEBench schema.
Not all plots have to be prepared for each dataset, as described in _Plots_ section.
All quantification challenges belong to the same benchmarking event.

## Inputs

### Examples

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :----------|
| 1 | BED | [Specification][spec-bed] | [Link][in1] | BED file with PAS quantification by the benchmarked tool |
| 2 | BED | [Specification][spec-bed] | [Link][in2] | BED file with PAS quantification from ground truth dataset |

### Additional info

#### Input 1

This BED file contains genomic positions of unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from RNAseq data by the benchmarked tool.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - the TPM value - quantified PAS
- **strand** - defines the strand; either "." (=no strand) or "+" or "-"

#### Input 2

This BED file contains genomic positions of unique unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from the orthogonal 3'end-seq dataset.
Fields are the same as in Input 1.

## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

1. **bar plot** visualizing **Correlation** of poly(A) site quantification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - benchmarked tool  
**Y axis** - correlation

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

Ranking: The best performing tool is the one with the highest correlation value.

Optional plots:

2.  **2D line plot** visualising **correlation** as a function of distance threshold.

**X axis** - correlation(d)  
**Y axis** - distance threshold _d_ 

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

Ranking: For each tool, the area under the curve (AUC) is calculated. The best performing tool is the one with the highest AUC.

Note: 2D line plot is not supported in OpenEBench yet. If it's not implemented, the data should be visualised outside of OpenEBench.

## Outputs

Calculated metrics are saved in JSON file adhering to OpenEBench schema. 
Assessment output is generated for each tool separately and contains values of calculated metrics for a given input dataset.
Consolidation output contains summarized data from all benchmarked tools within one challenge in a format suitable for plotting, e.g. single values for barplot or X,Y value pairs for 2D scatter plot.


| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :-------- |
| 1 | JSON | [Specification][spec-json] | [Link][assessment_out] | Assessment output JSON |
| 2 | JSON | [Specification][spec-json] | [Link][consolidation_out] | Consolidation output JSON |

### Additional info
 
#### Output 1

The OpenEBench assessment file contains the following attributes:

- **\_id** - follows the format: community:challenge\_metric\_tool
- **challenge_id** - describing the combination of input dataset with indication whether it is based on real or simulated data and any parameters used for metric calculation such as distance threshold (window size); e.g. datasetA\_simulated\_10nt
- **participant_id** - benchmarked tool
- **stderr** - standard error, set to zero
- **value** - metric value
- **metric_id** - metric name; metric names used in this benchmark are specified in the table below
 
The following tables list the metric names, value types and units, and a description of each attribute-value pair:

| Metric_id | Type | Unit | Description |
| :--- | :--- | :--- | :----------------- |
| `Correlation_coefficient_10nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 10 nt distance threshold (window) |
| `Correlation_coefficient_50nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 50 nt distance threshold (window) |
| `Correlation_coefficient_100nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 100 nt distance threshold (window) |

#### Output 2

The OpenEBench consolidation file specifies visualization types (2D scatter plot, barplot), descriptions of metrics used for X and Y axis and (X,Y) value pairs for each challenge.


[//]: # (References)
 
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[assessment_out]: ./example_files/assessment_out.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
