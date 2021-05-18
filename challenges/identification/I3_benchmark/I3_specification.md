# I3 benchmark specification

## Synopsis

Benchmark to compare output of tool's RNAseq based site identification to an established atlas of polyA-sites (PolyASite, PolyA_DB, etc.). The goal is to calculate precision with respect to data from existing databases for several different distance thresholds.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Poly(A) sites identified in PolyASite database
3. Poly(A) sites identified in PolyA_DB database

Based on the input data the following metrics are computed:

1. Precision = (TP/(TP+FP))

TP - true positives - PAS identified by the tool and present within X nucleotides of PAS identified in the database  
FP - false positives  - PAS identified by the tool and not present within X nucleotides of PAS identified in the database

The metrics should be calculated for the following thresholds of distance between PAS identified by the tool and PAS present in databases:

- 50 nt
- 100 nt
- 200 nt


## General info

* **Challenge:** Identification
* **Identifier:** I3

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | BED | [Specification][spec-bed-polyAsite] | [Link][in2] |
  | 3 | custom TXT | [Specification][spec-custom-polyAdb] | [Link][in3] |

### Additional info inputs
  
#### Format 1

This BED file contains PAS identified by the benchmarked tool from RNAseq data

#### Format 2

This BED file contains PAS from PolyASite database

#### Format 3

This BED file contains PAS from PolyA_DB database


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
  | `Precision_50nt_PolyASite` | `float` | N/A | Precision of PAS identification compared with PolyASite database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 50 nt |
  | `Precision_100nt_PolyASite` | `float` | N/A | Precision of PAS identification compared with PolyASite database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 100 nt |
  | `Precision_200nt_PolyASite` | `float` | N/A | Precision of PAS identification compared with PolyASite database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 200 nt |

#### Format 2
  
  The following table lists the attribute names, value types and units, and a
description of each attribute-value pair:
  
  | Attribute | Type | Unit | Description |
  | --- | --- | --- | --- |
  | `Precision_50nt_PolyAdb` | `float` | N/A | Precision of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 50 nt |
  | `Precision_100nt_PolyAdb` | `float` | N/A | Precision of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 100 nt |
  | `Precision_200nt_PolyAdb` | `float` | N/A | Precision of PAS identification compared with PolyA_DB database; Precision = (TP/(TP+FP)); maximum distance between identified PAS and PAS in database: 200 nt |
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Precision_50nt_PolyASite | N/A | Output 1 | Read file, parse JSON and extract attribute `Precision_50nt_PolyASite` | `vector` | N/A |
  | 2 | Precision_100nt_PolyASite | N/A | Output 1 | Read file, parse JSON and extract attribute `Precision_100nt_PolyASite` | `vector` | N/A |
  | 3 | Precision_200nt_PolyASite | N/A | Output 1 | Read file, parse JSON and extract attribute `Precision_200nt_PolyASite` | `vector` | N/A |
  | 4 | Precision_50nt_PolyAdb | N/A | Output 2 | Read file, parse JSON and extract attribute `Precision_50nt_PolyAdb` | `vector` | N/A |
  | 5 | Precision_100nt_PolyAdb | N/A | Output 2 | Read file, parse JSON and extract attribute `Precision_100nt_PolyAdb` | `vector` | N/A |
  | 6 | Precision_200nt_PolyAdb | N/A | Output 2 | Read file, parse JSON and extract attribute `Precision_200nt_PolyAdb` | `vector` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
  [in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[in3]: ./example_files/input3.bed

[out1]: ./example_files/output1.json
[out2]: ./example_files/output2.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
  [spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
  [spec-bed-polyAsite]: <https://polyasite.unibas.ch/atlas>
  [spec-custom-polyAdb]: <https://exon.apps.wistar.org/polya_db/v3/download/3.2/readme.txt>
  