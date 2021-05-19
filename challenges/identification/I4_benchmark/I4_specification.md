# I4 benchmark specification

## Synopsis

Benchmark to annotate the identified poly(A) sites to specific features (terminal exon, intron, 3'utr etc) and count the polyA sites assigned to those features.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. GTF/GFF file with genome annotation

Based on the input data the following metrics are computed:

1. Proportion of poly(A) sites annotated to 3'-UTR
2. Proportion of poly(A) sites annotated to introns
3. Proportion of poly(A) sites annotated to terminal exon
4. Proportion of poly(A) sites annotated to other features


## General info

* **Challenge:** Identification
* **Identifier:** I4

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | GTF | [Specification][spec-gtf] | [Link][in2] |

### Additional info inputs
  
#### Format 1

This BED file contains PAS identified by the benchmarked tool from RNAseq data

#### Format 2

This GTF file contains genome annotation for the RNAseq dataset

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
  | `annotation_3utr` | `float` | % | Proportion of poly(A) sites which were annotated to 3'-UTR |
  | `annotation_intron` | `float` | % | Proportion of poly(A) sites which were annotated to intron |
  | `annotation_terminal_exon` | `float` | % | Proportion of poly(A) sites which were annotated to terminal exon |
  | `annotation_other` | `float` | % | Proportion of poly(A) sites which were annotated to other features |
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Annotated_to_3utr | % | Output 1 | Read file, parse JSON and extract attribute `annotation_3utr` | `vector` | N/A |
  | 2 | Annotated_to_introns | % | Output 1 | Read file, parse JSON and extract attribute `annotation_intron` | `vector` | N/A |
  | 3 | Annotated_to_terminal_exon | % | Output 1 | Read file, parse JSON and extract attribute `annotation_terminal_exon` | `vector` | N/A |
  | 4 | Annotated_to_other | % | Output 1 | Read file, parse JSON and extract attribute `annotation_other` | `vector` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
  [in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.gtf
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
  [spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
  [spec-gtf]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format4>
  