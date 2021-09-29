# D2 benchmark specification

## Synopsis

Benchmark to test sensitivity and False Discovery Rate (FDR) of RNAseq-based identification of genes with differentially used poly(A) sites compared with ground truth (genes identified from orthogonal 3'end seq data).

Input data:

1. List of genes with information on differentially used poly(A) sites identified from RNA-Seq data using the benchmarked tool
2. List of genes with information on differentially used poly(A) sites identified from orthogonal 3'end seq dataset

Based on the input data the following metrics are computed:

1. Sensitivity = (TP/(TP+FN))
2. FDR = (FP/(TP+FP))

TP - true positives - genes identified by the tool and present in the orthogonal dataset  
FP - false positives - genes identified by the tool and not present in the orthogonal dataset  
FN - false negatives - genes not identified by the tool but present in the orthogonal dataset


## General info

* **Challenge:** Differential usage
* **Identifier:** D2

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | CSV | [Specification][spec-csv] | [Link][in1] |
  | 2 | CSV | [Specification][spec-csv] | [Link][in2] |

### Additional info inputs
  
#### Format 1

This CSV file contains a list of genes with differentially used poly(A) sites identified from RNA-Seq data by the benchmarked method.  
Columns:

- gene name
- is_de: information whether **any** PAS in the gene was differentially expressed; possible values: "0" if no differential expression detected, "1" if differentially expressed
- is_lengthened: information whether shortening/lengthening events were detected for the gene; possible values: "0" if no shortening/lengthening events detected, "1" if shortening/lengthening events detected

#### Format 2

This CSV file contains a list of genes with differentially used poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Columns:

- gene name
- is_de: information whether **any** PAS in the gene was differentially expressed; possible values: "0" if no differential expression detected, "1" if differentially expressed
- is_lengthened: information whether shortening/lengthening events were detected for the gene; possible values: "0" if no shortening/lengthening events detected, "1" if shortening/lengthening events detected


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
  | `sensitivity` | `float` | N/A | Sensitivity of detection of differentially expressed PAS compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)) |
  | `FDR` | `float` | N/A | False Discovery Rate of detection of differentially expressed PAS compared with orthogonal dataset; FDR = (FP/(TP+FP)) |
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Sensitivity | N/A | Output 1 | Read file, parse JSON and extract attribute `sensitivity`. | `vector` | N/A |
  | 2 | FDR | N/A | Output 1 | Read file, parse JSON and extract attribute `FDR`. | `vector` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
  [in1]: ./example_files/input1.csv
[in2]: ./example_files/input2.csv
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-csv]: <https://en.wikipedia.org/wiki/Comma-separated_values>

  
