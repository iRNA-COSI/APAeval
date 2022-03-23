# Identification benchmarks specification

## Synopsis

Benchmarks to assess the performance of RNAseq-based poly(A) site (PAS) identification compared with ground truth (poly(A)-sites identified from orthogonal 3'end seq data or existing PAS databases).

### Input data:

Comparison of PAS predicted from RNA-Seq data with 3'end sequencing data

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites identified based on orthogonal 3'end seq dataset

Comparison of PAS predicted from simulated RNA-Seq data with dataset used for simulation

1. Poly(A) sites identified based on simulated RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites of transcripts used to simulate RNAseq reads

### Metrics:

Based on the input data the following metrics are computed:

1. Sensitivity (TPR, True Positive Rate) = TP/(TP+FN)
2. FPR (False Positive Rate) = FP/(FP+TN)
3. Precision = TP/(TP+FP)
4. Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for a range of distance thresholds

TP - true positives - PAS identified by the tool and present in the orthogonal dataset  
FP - false positives - PAS identified by the tool and not present in the orthogonal dataset  
FN - false negatives - PAS not identified by the tool but present in the orthogonal dataset

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the orthogonal 3'end seq dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.

### OpenEBench challenges

The metrics are visualised using 2D scatter plots and barplots, as described in _Plots_ section, which are then used for ranking the participating tools.
A plot for any given input dataset constitutes a benchmarking challenge as understood in OpenEBench schema.
Not all plots have to be prepared for each dataset, as described in _Plots_ section.
All identification challenges belong to the same benchmarking event.

## Inputs

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :----------|
| 1 | BED | [Specification][spec-bed] | [Link][in1] | BED file with PAS identification by the benchmarked tool |
| 2 | BED | [Specification][spec-bed] | [Link][in2] | BED file with PAS identification from ground truth dataset |

### Additional info
  
#### Input 1

This BED file contains single-nucleotide position of poly(A) sites identified by the benchmarked method.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Input 2

This BED file contains single-nucleotide position of poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Fields are the same as in Input 1.
	
## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

1. **bar plot** visualizing **Precision** of poly(A) site identification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - benchmarked tool  
**Y axis** - precision

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation


Ranking: The best performing tool is the one with the highest precision value.

2. **2D scatter plot** visualizing **TPR and FPR** of poly(A) site identification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - FPR(d)  
**Y axis** - TPR(d)  
where _d - distance threshold_

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

Ranking: The best performing tool is the one with the highest TPR combined with lowest FPR (top left part of the plot) and the worst performing tool is the one with the lowest TPR combined with highest FPR (bottom left part of the plot. The plot should be divided into diagonal quartiles based on the distance from optimal performance. Alternatively, if the plot is divided into square quartiles, the following ranking order should be applied: top-left, top-right, bottom-left, bottom-right.

3. **bar plot** visualizing **AUC** of ROC plot with single AUC value for each tool.

**X axis** - benchmarked tool  
**Y axis** - AUC

Input data:

- RNA-Seq data compared with 3'end sequencing data


Ranking: The best performing tool is the one with the highest AUC value.

Optional plots:

4. **2D line plot** visualising **TPR and FPR** in the form of [ROC curve][roc-curve], i.e. TPR and FPR calculated for a range of distance thresholds are plotted together. 

**X axis** - FPR(d)  
**Y axis** - TPR(d)  
where _d - distance threshold_ in range 0-100 nt with 10 nt increments.

Input datasets:

- RNA-Seq data compared with 3'end sequencing data

Ranking: For each tool, the area under the curve is calculated. The best performing tool is the one with the highest AUC.

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
- **challenge_id** - describing the combination of input dataset with indication whether it is based on real or simulated data and any parameters used for metric calculation; e.g. datasetA\_simulated\_10nt
- **participant_id** - benchmarked tool
- **stderr** - standard error, set to zero
- **value** - metric value
- **metric_id** - metric name; metric names used in this benchmark are specified in the table below
 
The following tables list the metric names, value types and units, and a description of each attribute-value pair:

| Attribute | Type | Unit | Description |
| :--- | :--- | :--- | :----------------- |
| `Precision_10nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 10 nt distance threshold |
| `Precision_50nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 50 nt distance threshold |
| `Precision_100nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 100 nt distance threshold |
| `Sensitivity_10nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 10 nt distance threshold |
| `Sensitivity_50nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 50 nt distance threshold |
| `Sensitivity_100nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 100 nt distance threshold |
| `FPR_10nt` | `float` | N/A | False Positive Rate of PAS identification compared with orthogonal dataset; FPR = FP/(FP+TN); calculated for 10 nt distance threshold |
| `FPR_50nt` | `float` | N/A | False Positive Rate of PAS identification compared with orthogonal dataset; FPR = FP/(FP+TN); calculated for 50 nt distance threshold |
| `FPR_100nt` | `float` | N/A | False Positive Rate of PAS identification compared with orthogonal dataset; FPR = FP/(FP+TN); calculated for 100 nt distance threshold |
| `AUC` | `float` | N/A | Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for all distance thresholds between 0 nt and 100 nt with 10 nt intervals. |

#### Output 2

The OpenEBench consolidation file specifies visualization types (2D scatter plot, barplot), descriptions of metrics used for X and Y axis and (X,Y) value pairs for each challenge.

[//]: # (References)
  
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[assessment_out]: ./example_files/assessment_out.json
[consolidation_out]: ./example_files/consolidation_out.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[roc-curve]: <https://en.wikipedia.org/wiki/Receiver_operating_characteristic#ROC_curves_beyond_binary_classification>
  
