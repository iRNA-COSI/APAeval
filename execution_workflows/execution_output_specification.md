# Execution workflows output specification

## Synopsis

This specification describes required output of `execution_workflows` to ensure that the `summary_workflows` have access to files of the same format from each participant.

Outputs depend on features available for the participant (i.e. not all participants perform *de novo* identification of PAS) and can be grouped into three categories:

- identification benchmarking event:
  - BED file with identified poly(A) sites with single nucleotide resolution
- quantification benchmarking event:
  - BED file with identified unique poly(A) sites and TPM values for each site
  - BED file with identified unique poly(A) sites and relative usage values for each site
  - BED file with coordinates of regions (e.g. terminal exons, genes) containing poly(A) sites
  - BED file with coordinates of regions and the relative usage value for each region
- differential usage benchmarking event:
  - TSV file with gene ID and significance of differential PAS usage

## Inputs

Inputs to execution workflows are provided by APAeval.
>Detailed specification currently missing. Please refer to the [execution workflow README][ex-readme-in]
## Outputs


| OUTCODE | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 01 | BED | [Specification][spec-bed] | [Link][out1] |
  | 02 | BED | [Specification][spec-bed] | [Link][out2] |
  | 03 | TSV | [Wikipedia][wiki-tsv] | [Link][out3] |
  | 04 | BED | [Specification][spec-bed] | [Link][out4] |
  | 05 | BED | [Specification][spec-bed] | [Link][out5] |
  | 06 | BED | [Specification][spec-bed] | [Link][out6] |
### Additional info inputs

#### Format 01

This BED file contains single-nucleotide position of poly(A) sites identified by the participant.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 02

This BED file contains positions of unique poly(A) sites with TPM values for each identified site in the **score** column.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - TPM value for the identified site
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 03

This TSV file contains two columns:

1. gene ID
2. significance of differential PAS usage

Column names should not be added to the file.

#### Format 04

This BED file contains positions of unique poly(A) sites with relative usage values for each identified site in the **score** column.

Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - relative usage value for the identified site
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 05

This BED file contains positions of regions (e.g. terminal exons of genes, whole genes) that contain identified poly(A) sites.

Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the first nucleotide of the region (e.g. terminal exon, gene); the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; this corresponds to the last nucleotide of the region.
- **name** - defines the name of the identified region. It's recommended to use a conventional identifier (e.g. Ensembl transcript ID, gene ID)
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 06

This BED file contains positions of regions (e.g. terminal exons of genes, whole genes) and the relative usage values for each identified region in the **score** column.

Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the first nucleotide of the region (e.g. terminal exon, gene); the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; this corresponds to the last nucleotide of the region.
- **name** - defines the name of the identified region. It's recommended to use a conventional identifier (e.g. Ensembl transcript ID, gene ID)
- **score** - relative usage value for the identified site
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

## Naming conventions
For naming conventions please refer to the [execution workflow README][ex-readme]

[//]: # (References)

[out1]: ./example_output_files/output1.bed
[out2]: ./example_output_files/output2.bed
[out3]: ./example_output_files/output3.tsv
[out4]: ./example_output_files/output4.bed
[out5]: ./example_output_files/output5.bed
[out6]: ./example_output_files/output6.bed
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[wiki-tsv]: <https://en.wikipedia.org/wiki/Tab-separated_values>
[ex-readme]: ./README.md
[ex-readme-in]: ./README.md#more-details
