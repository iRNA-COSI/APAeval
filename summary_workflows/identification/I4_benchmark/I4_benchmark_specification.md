# I4 benchmark specification

## Synopsis

Benchmark to assign the identified poly(A) sites to specific features (5'-UTR, 3'-UTR, CDS, introns and intergenic regions) and count the proportion of polyA sites assigned to those features.

Input data:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. GFF file with genome annotation

Based on the input data the proportion of poly(A) sites assigned to the following genomic features should be computed:

1. 3'-UTRs
2. 5'-UTRs
3. introns
4. CDS
5. terminal exons
6. intergenic regions

Additionally the proportion of poly(A) sites assigned to regions downstream of CDS should be computed:

7. 1kb downstream of CDS
8. 5kb downstream of CDS
9. 10kb downstream of CDS


## General info

* **Challenge:** Identification
* **Identifier:** I4

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | BED | [Specification][spec-bed] | [Link][in1] |
  | 2 | GFF | [Specification][spec-gtf] | [Link][in2] |

### Additional info inputs
  
#### Format 1

This BED file contains single-nucleotide position of poly(A) sites identified by the benchmarked method.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 2

This GFF file contains reference genome annotation for the RNAseq dataset.

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
  | `annotation_3utr` | `float` | % | Proportion of poly(A) sites which were annotated to 3'-UTRs |
  | `annotation_5utr` | `float` | % | Proportion of poly(A) sites which were annotated to 5'-UTRs |
  | `annotation_intron` | `float` | % | Proportion of poly(A) sites which were annotated to introns |
  | `annotation_cds` | `float` | % | Proportion of poly(A) sites which were annotated to CDS |
  | `annotation_terminal_exon` | `float` | % | Proportion of poly(A) sites which were annotated to terminal exon |
  | `annotation_intergenic` | `float` | % | Proportion of poly(A) sites which were annotated to intergenic regions |
  | `annotation_ds_1kb` | `float` | % | Proportion of poly(A) sites which were annotated to region 1kb downstream from CDS |
  | `annotation_ds_5kb` | `float` | % | Proportion of poly(A) sites which were annotated to region 5kb downstream from CDS |
  | `annotation_ds_10kb` | `float` | % | Proportion of poly(A) sites which were annotated to region 10kb downstream from CDS |
  
## Metrics
  
  | # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | PAS annotated to 3'-UTRs | % | Output 1 | Read file, parse JSON and extract attribute `annotation_3utr` | `vector` | N/A |
  | 2 | PAS annotated to 5'-UTRs | % | Output 1 | Read file, parse JSON and extract attribute `annotation_5utr` | `vector` | N/A |
  | 3 | PAS annotated to introns | % | Output 1 | Read file, parse JSON and extract attribute `annotation_intron` | `vector` | N/A |
  | 4 | PAS annotated to CDS | % | Output 1 | Read file, parse JSON and extract attribute `annotation_cds` | `vector` | N/A |
  | 5 | PAS annotated to terminal exon | % | Output 1 | Read file, parse JSON and extract attribute `annotation_terminal_exon` | `vector` | N/A |
  | 6 | PAS annotated to intergenic region | % | Output 1 | Read file, parse JSON and extract attribute `annotation_intergenic` | `vector` | N/A |
  | 7 | PAS annotated to region 1kb downstream from CDS | % | Output 1 | Read file, parse JSON and extract attribute `annotation_ds_1kb` | `vector` | N/A |
  | 8 | PAS annotated to region 5kb downstream from CDS | % | Output 1 | Read file, parse JSON and extract attribute `annotation_ds_5kb` | `vector` | N/A |
  | 9 | PAS annotated to region 10kb downstream from CDS | % | Output 1 | Read file, parse JSON and extract attribute `annotation_ds_10kb` | `vector` | N/A |
  
### Additional info metrics
  
  N/A

[//]: # (References)
  
  [in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.gff
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
  [spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
  [spec-gtf]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format3>
  