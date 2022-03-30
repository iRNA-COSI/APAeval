# I3 benchmark specification

## Synopsis

Benchmark to compare output of tool's RNAseq based site identification to an established atlas of polyA-sites (PolyASite, PolyA_DB). The goal is to calculate precision with respect to data from existing databases for several different distance thresholds.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites identified in PolyASite database
3. Poly(A) sites identified in PolyA_DB database

Based on the input data the following metrics are computed:

1. Sensitivity = TP/(TP+FN)
2. FDR = FP/(TP+FP)
3. FPR = FP/(FP+TN)
4. Precision = TP/(TP+FP)

TP - true positives - PAS identified by the tool and present within X nucleotides of PAS identified in the database  
FP - false positives  - PAS identified by the tool and not present within X nucleotides of PAS identified in the database

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified in the database i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified in the database for the prediction to be considered true.  
Distance thresholds should be between 0 nt and 100 nt with 10 nt increments.

Optional metrics:

5. Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for all distance thresholds


## General info

* **Challenge:** Identification
* **Identifier:** I3

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | BED | [Specification][spec-bed-polyAsite] | [Link][in2] |
  | 3 | custom TSV | [Specification][spec-custom-polyAdb] | [Link][in3] |

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

This BED file contains PAS from PolyASite database

#### Format 3

This TSV file contains PAS from PolyA_DB database

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
  | 2 | JSON | [Specification][spec-json] | [Link][out2] |
  
### Additional info outputs

#### Format 1
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `sensitivity` | `vector` | N/A | A vector of length=11 containing sensitivity values of PAS identification compared with database; Sensitivity = (TP/(TP+FN)); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `FDR` | `vector` | N/A | A vector of length=11 containing False Discovery Rate values of PAS identification compared with database; FDR = (FP/(TP+FP)); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `FPR` | `vector` | N/A | A vector of length=11 containing False Positive Rate values of PAS identification compared with database; FPR = FP/(FP+TN); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with database; Precision = TP/(TP+FP); calculated for distance between 0 nt and 100 nt with 10 nt intervals; Each value in the vector is of type `float` |
  | `AUC` | `float` | N/A | Area under the curve (AUC) of receiver operating characteristic (ROC) calculated from TPR and FPR values for all distance thresholds between 0 nt and 100 nt with 10 nt intervals. |



#### Format 2
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `Precision` | `vector` | N/A | A vector of length=11 containing precision values of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float` |

A vector of length=11 containing precision values of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); calculated for distance between 0 nt and 200 nt with 20 nt intervals; Each value in the vector is of type `float`
  
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
[in3]: ./example_files/input3.tsv

[out1]: ./example_files/output1.json
[out2]: ./example_files/output2.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[spec-bed-polyAsite]: <https://polyasite.unibas.ch/atlas>
[spec-custom-polyAdb]: <https://exon.apps.wistar.org/polya_db/v3/download/3.2/readme.txt>
[roc-curve]: <https://en.wikipedia.org/wiki/Receiver_operating_characteristic#ROC_curves_beyond_binary_classification>
 
