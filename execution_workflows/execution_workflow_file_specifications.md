# Execution workflows output specification

## Synopsis

This specification describes allowed input and required output of `execution_workflows` to ensure that the `summary_workflows` have access to files of the same format from each participant.

Outputs depend on features available for the participant (i.e. not all participants perform *de novo* identification of PAS) and can be grouped into four categories:

- identification benchmarking event:
  - BED file with identified poly(A) sites with single nucleotide resolution
- absolute quantification benchmarking event:
  - BED file with identified unique poly(A) sites and TPM values for each site
- relative quantification benchmarking event:
  - BED file with identified unique poly(A) sites and relative usage values for each site
- differential usage benchmarking event (not yet implemented)
  - TSV file with gene ID and significance of differential PAS usage

## Inputs

Inputs to execution workflows are provided by APAeval.

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :----------|
| 1 | BAM | [Specification][spec-bam] | - | BAM file with mapped RNA seq reads<br/>Pre-processing done with [nfcore-rnaseq][nfcore-rnaseq]
| 2 | GTF | [Specification][spec-gtf] | [Link][in2] | GTF file with genome annotation including 3' UTRs |

## Outputs


| OUTCODE | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 01 | BED | [Specification][spec-bed] | [Link][out-i] |
  | 02 | BED | [Specification][spec-bed] | [Link][out-aq] |
  | 03 | TSV | [Wikipedia][wiki-tsv] | Not yet implemented |
  | 04 | BED | [Specification][spec-bed] | [Link][out-rq] |

### Additional info outputs

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

This BED file contains positions of unique poly(A) sites with relative usage values for each identified site in the **score** column. The relative usage values must be fractional for each poly(A) site (e.g. fraction used vs other sites in the same terminal exon).

Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - relative usage value for the identified site (between 0-1)
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".


## Naming conventions
For naming conventions please refer to the [execution workflow README][ex-readme]

[//]: # (References)

[in2]: ./example_files/input2.gtf
[nfcore-rnaseq]: <https://nf-co.re/rnaseq>
[out-aq]: ./example_files/output_absQuant.bed
[out-i]: ./example_files/output_identification.bed
[out-rq]: ./example_files/output_relQuant.bed
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[spec-bam]: <https://samtools.github.io/hts-specs/SAMv1.pdf>
[spec-gtf]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format4>
[wiki-tsv]: <https://en.wikipedia.org/wiki/Tab-separated_values>
[ex-readme]: ./README.md
[ex-readme-in]: ./README.md#more-details
