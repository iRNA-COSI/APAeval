# Quantification benchmarks specification

- [Synopsis](#synopsis)
  - [Input data:](#input-data)
  - [Metrics](#metrics)
    - [Poly(A) site matching](#polya-site-matching)
    - [Signature metrics](#signature-metrics)
    - [All metrics](#all-metrics)
  - [OpenEBench challenges](#openebench-challenges)
- [Inputs](#inputs)
  - [Input 1](#input-1)
  - [Input 2](#input-2)
  - [Input 3](#input-3)
- [Plots](#plots)
  - [Correlation vs. TPM FP](#correlation-vs-tpm-fp)
  - [Correlation of relative usage](#correlation-of-relative-usage)
  - [Optional plots](#optional-plots)
    - [Correlation by window size](#correlation-by-window-size)
- [Outputs](#outputs)
  - [Output 1](#output-1)
  - [Output 2](#output-2)
  - [Output 3](#output-3)
## Synopsis

Benchmarks to assess the performance of RNA-Seq-based **poly(A) site quantification** compared to poly(A)-sites quantification based on orthogonal 3'end seq data.

The metrics are computed for two types of datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### Input data:

Comparison of PAS predicted from **RNA-Seq data with 3'end sequencing data**

1. Poly(A) sites quantification (TPM) based on RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites quantification based on orthogonal 3'end-seq dataset
3. Genome annotation file in GTF format

Comparison of PAS predicted from **simulated RNA-Seq data with dataset used for simulation**
  
1. Poly(A) sites quantification (TPM) based on simulated RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites quantification based on abundances of transcripts used to simulate the RNAseq data
3. Genome annotation file in GTF format

### Metrics
#### Poly(A) site matching
To compute the metrics, the poly(A) sites identified from the RNAseq data by the benchmarked tool (PD) and from the ground-truth dataset (GT) should be first mapped to each other as [implemented in the `matchPAS` module][matchpas]. Briefly
- Matching of PD and GT sites is performed for different window sizes (specified below); i.e. a PD site is considered matched if its coordinates are within [GT - window, GT + window]
  - If a PD site matches multiple GT sites its score is split between the GT sites according to distance
  - If multiple PD sites match one GT they are merged and their score is summed
- False negatives (FN) are GT sites without a PD match
- False positives (FP) are PD sites without a GT match

#### Signature metrics
The following metrics represent the core set for the absolute quantification event:
 
1. Pearson correlation between RNAseq-based quantification and 3'end-seq quantification.
  
The correlation is calculated based on GT and PD sites that overlap for a given window size (referred to as "matched sites", or "intersection")

2. Percent expression of non-matched sites

The expression values of all PD sites that were not mapped to any GT site (= FP) are summed together. As the PAS expression levels are provided as TPM, these values are relative to the expression of all identified PAS. For readability however, "non-matched" TPM are reported as percentage of "total" TPM.
  
3. Correlation of relative PAS usage calculated from RNAseq-based PAS quantification and orthogonal 3'end seq data

The PD sites are mapped to GT sites as in point 1.
The identified PAS are then assigned to genes based on genome annotation, and relative PAS usage is calculated separately for each gene.
PAS that cannot be assigned to any genes are discarded.
Correlation is then calculated globally for all poly(A) sites assigned to genes.

#### All metrics

| Metric_id | Type | Unit | Description |
| :--- | :--- | :--- | :----------------- |
| `Sum_FP_TPM:10nt` | `float` | N/A | Total expression of non-matched sites quantified by benchmarked tool for 10 nt distance threshold (window) |
| `Sum_FP_TPM:50nt` | `float` | N/A | Total expression of non-matched sites quantified by benchmarked tool for 50 nt distance threshold (window) |
| `Sum_FP_TPM:100nt` | `float` | N/A | Total expression of non-matched sites quantified by benchmarked tool for 100 nt distance threshold (window) |
| `Percent_FP_TPM:10nt` | `float` | N/A | Percent total expression of non-matched sites for total expression of non-matched and matched sites quantified by benchmarked tool for 10 nt distance threshold (window) |
| `Percent_FP_TPM:50nt` | `float` | N/A | Percent total expression of non-matched sites for total expression of non-matched and matched sites quantified by benchmarked tool for 50 nt distance threshold |
| `Percent_FP_TPM:100nt` | `float` | N/A | Percent total expression of non-matched sites for total expression of non-matched and matched sites quantified by benchmarked tool for 100 nt distance threshold |
| `Sensitivity:10nt` | `float` | N/A | Sensitivity (recall, hit rate, true positive rate) = TP / (TP + FN) for 10 nt distance threshold |
| `Sensitivity:50nt` | `float` | N/A | Sensitivity (recall, hit rate, true positive rate) = TP / (TP + FN) for 50 nt distance threshold |
| `Sensitivity:100nt` | `float` | N/A | Sensitivity (recall, hit rate, true positive rate) = TP / (TP + FN) for 100 nt distance threshold |
| `Precision:10nt` | `float` | N/A | Precision (positive predictive value) = TP / (TP + FP) for 10 nt distance threshold |
| `Precision:50nt` | `float` | N/A | Precision (positive predictive value) = TP / (TP + FP) for 50 nt distance threshold |
| `Precision:100nt` | `float` | N/A | Precision (positive predictive value) = TP / (TP + FP) for 100 nt distance threshold |
| `F1_score:10nt` | `float` | N/A | F1 score, harmonic mean of precision and sensitivity for 10 nt distance threshold |
| `F1_score:50nt` | `float` | N/A | F1 score, harmonic mean of precision and sensitivity for 50 nt distance threshold |
| `F1_score:100nt` | `float` | N/A | F1 score, harmonic mean of precision and sensitivity for 100 nt distance threshold |
| `Jaccard_index:10nt` | `float` | N/A | Jaccard index (Jaccard similarity coefficient) = TP / (TP + FP + FN) for 10 nt distance threshold |
| `Jaccard_index:50nt` | `float` | N/A | Jaccard index (Jaccard similarity coefficient) = TP / (TP + FP + FN) for 50 nt distance threshold |
| `Jaccard_index:100nt` | `float` | N/A | Jaccard index (Jaccard similarity coefficient) = TP / (TP + FP + FN) for 100 nt distance threshold |
| `Pearson_r:all_GT:10nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 10 nt distance threshold (window) |
| `Pearson_r:all_GT:50nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 50 nt distance threshold (window) |
| `Pearson_r:all_GT:100nt` | `float` | N/A | Correlation between poly(A) site quantification by benchmarked tool and ground truth dataset calculated for 100 nt distance threshold (window) |
| `Pearson_r_relative:union:100nt` | `float` | N/A | Correlation between relative PAS usage calculated from RNAseq-based PAS quantification and orthogonal 3'end seq data for 100 nt distance threshold (window) |

> The above mentioned correlation metrics are actually calculated for three different matching types - `all_GT` (this is equivalent to TP + FN), `union` (TP + FP + FN) and `intersection` (TP), and for *Pearson*'s as well as *Spearman*'s r. 

### OpenEBench challenges

This document describes the specifications for the APAeval "absolute quantification **benchmarking event**", which is partially compatible with the [OpenEBench (OEB) benchmarking plattform][oeb]. The event consists of several **challenges**. In APAeval, a challenge as understood in the OEB terminology is determined by an input *dataset* (one or more *samples*). A participant submits one file ([Input 1](#input-1) as defined below) for each challenge it participates in. The participant's result is compared to the challenge's ground truth ([Input 2](#input-2)) and the metrics described here are stored in OEB specific JSON files.


## Inputs

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :----------|
| 1 | BED | [Specification][spec-bed] | [Link][in1] | BED file with PAS quantification by the benchmarked tool |
| 2 | BED | [Specification][spec-bed] | [Link][in2] | BED file with PAS quantification from ground truth dataset |
| 3 | GTF | [Specification][spec-gtf] | [Link][in3] | GTF file with genome annotation |

### Input 1

This BED file contains genomic positions of unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from RNAseq data by the benchmarked tool.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - the TPM value - quantified PAS
- **strand** - defines the strand; either "." (=no strand) or "+" or "-"

### Input 2

This BED file contains genomic positions of unique cleavage/polyadenylation sites as well as TPM values for each identified site quantified from the orthogonal 3'end-seq dataset.
Fields are the same as in Input 1.

### Input 3

The GTF file is a standard genome annotation file.

## Plots

The results of this benchmark will be visualised in OEB using the following plots:

### Correlation vs. TPM FP
2D scatter plot visualizing correlation of poly(A) site quantification and percentage of total expression of non-matched sites. 

**Plot type**: 2D scatter plot

**Metric X**: Percent expression of non-matched sites  
**Metric Y**: Correlation

**Ranking**: The best performing tool is the one with the highest correlation value and lowest total expression of non-matched sites (top-left part of the plot) and the worst performing tool is the one with the lowest correlation combined with the highest total expression of non-matched sites (bottom-right part of the plot). The plot should be divided into diagonal quartiles based on the distance from optimal performance.

Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### Correlation of relative usage
Bar plot visualizing correlation of relative PAS usage on gene level. 

**Plot type**: Bar plot

**Metric**: Correlation of relative PAS usage

**Ranking**: The best performing tool is the one with the highest correlation value value.

Separate plots should be prepared for different values of distance threshold:

- 100 nt

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### Optional plots

This section includes suggestions for plots that can be visualised outside of OEB. Additional suggestions, especially using the metrics that are computed already (see below), are always welcome.

#### Correlation by window size
2D line plot visualising correlation as a function of distance threshold.

**Plot type**: 2D line plot

**Metric X**: Correlation(d)  
**Metric Y**: Distance threshold _d_ 

**Ranking**: None

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

Note: 2D line plot is not supported in OEB yet. The metrics should be visualised only outside of OEB.

## Outputs

Calculated metrics are saved in JSON files adhering to OEB schema. 
Assessment output is generated for each tool separately and contains values of calculated metrics for a given input dataset.
Consolidation output contains summarized data from all benchmarked tools within one challenge in a format suitable for plotting, e.g. single values for barplot or X,Y value pairs for 2D scatter plot.


| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :-------- |
| 1 | JSON | [Specification][spec-json] | [Link][assessment_out] | Assessment output JSON |
| 2 | JSON | [Specification][spec-json] | [Link][aggregation_out] | Aggregation output JSON |
| 3 | JSON | [Specification][spec-json] | N/A | Consolidation output JSON |

 
### Output 1

The OEB assessment file contains the following attributes:

- **\_id** - follows the format: community:challenge\_metric\_tool
- **challenge_id** - describing the combination of input dataset with indication whether it is based on real or simulated data and any parameters used for metric calculation such as distance threshold (window size); e.g. datasetA\_simulated\_10nt
- **participant_id** - benchmarked tool
- **metrics**:
	- **value** - metric value
	- **metric_id** - metric name; metric names used in this benchmark are specified [below](#list-of-calculated-metrics).
  
 The assessment file will contain **all** metric variants desribed [above](#all-metrics).



### Output 2
The OEB aggregation file (for plotting in OEB) will only contain the [signature metrics](#signature-metrics) described above. Those have to be explicitly specified for visualization in the [aggregation template][aggr_temp].

### Output 3

The OEB consolidation file contains all validation, assessment and aggregation objects for all challenges a participant has been run on.


[//]: # (References)
 
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[in3]: ./example_files/input3.gtf
[assessment_out]: ./example_files/assessment_out.json
[aggregation_out]: ./example_files/aggregation_out.json
[aggr_temp]: ../quantification_dockers/q_consolidation/aggregation_template_Q.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[spec-gtf]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format4>
[matchpas]: ../../../utils/matchPAS/
[oeb]: <https://openebench.bsc.es/>