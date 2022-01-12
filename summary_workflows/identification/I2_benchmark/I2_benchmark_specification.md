# I2 benchmark specification

## Synopsis

Benchmark to test sensitivity and False Discovery Rate (FDR) of RNAseq-based poly(A) site identification compared with ground truth (poly(A)-sites identified from orthogonal 3'end seq data).

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites identified based on orthogonal 3'end seq dataset

Based on the input data the following metrics are computed:

1. Sensitivity = TP/(TP+FN)
2. FDR = FP/(TP+FP)
3. FPR = FP/(FP+TN)
4. Precision = TP/(TP+FP)

TP - true positives - PAS identified by the tool and present in the orthogonal dataset  
FP - false positives - PAS identified by the tool and not present in the orthogonal dataset  
FN - false negatives - PAS not identified by the tool but present in the orthogonal dataset

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the orthogonal 3'end seq dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.  
Distance thresholds should be between 0 nt and 100 nt with 10 nt increments.

Optional metrics:

5. Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for all distance thresholds


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
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 2

This BED file contains single-nucleotide position of poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Fields are the same as in Format 1.
	
## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

1. **bar plot** visualizing **Precision** of poly(A) site identification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - benchmarked tool  
**Y axis** - precision

Ranking: The best performing tool is the one with the highest precision value.

2. **bar plot** visualizing **Sensitivity** of poly(A) site identification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - benchmarked tool  
**Y axis** - sensitivity

Ranking: The best performing tool is the one with the highest sensitivity value.

3. **2D scatter plot** visualizing **TPR and FPR** of poly(A) site identification. Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

**X axis** - FPR(d)  
**Y axis** - TPR(d)  
where _d - distance threshold_

Ranking: The best performing tool is the one with the highest TPR combined with lowest FPR (top left part of the plot) and the worst performing tool is the one with the lowest TPR combined with highest FPR (bottom left part of the plot. The plot should be divided into diagonal quartiles based on the distance from optimal performance. Alternatively, if the plot is divided into square quartiles, the following ranking order should be applied: top-left, top-right, bottom-left, bottom-right.

4. **2D line plot** visualising **TPR and FPR** in the form of [ROC curve][roc-curve], i.e. TPR and FPR calculated for a range of distance thresholds are plotted together.

**X axis** - FPR(d)  
**Y axis** - TPR(d)  
where _d - distance threshold_

Ranking: For each tool, the area under the curve is calculated. The best performing tool is the one with the highest AUC.

Note: 2D line plot is not supported in OpenEBench yet. If it's not implemented, the data should be visualised outside of OpenEBench.

Optional plots:

5. **bar plot** visualizing **AUC** of ROC plot with single AUC value for each tool.

**X axis** - benchmarked tool  
**Y axis** - AUC

Ranking: The best performing tool is the one with the highest AUC value.


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
  | `FPR` | `vector` | N/A | A vector of length=11 containing False Positive Rate values of PAS identification compared with orthogonal dataset; FPR = FP/(FP+TN); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `AUC` | `float` | N/A | Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for all distance thresholds between 0 nt and 100 nt with 10 nt intervals. |
 
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Sensitivity | N/A | Output 1 | Read file, parse JSON and extract attribute `sensitivity` that has type `vector` of `float`.  | `array` | N/A |
  | 2 | FDR | N/A | Output 1 | Read file, parse JSON and extract attribute `FDR` that has type `vector` of `float`.  | `array` | N/A |
  | 3 | FDR | N/A | Output 1 | Read file, parse JSON and extract attribute `FPR` that has type `vector` of `float`. | `array` | N/A |
  | 4 | Precision | N/A | Output 1 | Read file, parse JSON and extract attribute `FDR` that has type `vector` of `float`. | `array` | N/A |
  | 5 | AUC | N/A | Output 1 | Read file, parse JSON and extract attribute `FDR` that has type `float`. | `vector` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[roc-curve]: <https://en.wikipedia.org/wiki/Receiver_operating_characteristic#ROC_curves_beyond_binary_classification>
  
