# Execution pipelines output specification

## Synopsis

This specification describes required output of `execution_workflows` to ensure that the `summary_workflows` have access to files in the same format from each tool.

Outputs depend on features available for the tool (i.e. not all tools perform *de novo* identification of PAS) and can be grouped into three categories:

- identification challenge:
  - BED file with identified poly(A) sites with single nucleotide resolution
- quantification challenge:
  - BED file with identified unique poly(A) sites and TPM values for each site
- differential usage challenge:
  - TSV file with gene ID and significance of differential PAS usage

## Inputs

Inputs to the execution pipeline depend on the method/tool.

## Outputs


| OUTCODE | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 01 | BED | [Specification][spec-bed] | [Link][out1] |
  | 02 | BED | [Specification][spec-bed] | [Link][out2] |
  | 03 | TSV | [Wikipedia][wiki-tsv] | [Link][out3] |
  
### Additional info inputs
  
#### Format 01

This BED file contains single-nucleotide position of poly(A) sites identified by the tool.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 02

This BED file contains positions of unique poly(A) sites with TPM values for each identified site in the **score** column.

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the name of the identified poly(A) site
- **score** - TPM value for the identified site
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Format 03

This TSV file contains two columns:

1. gene ID
2. significance of differential PAS usage

Column names should not be added to the file.


## Naming conventions
For naming conventions please refer to the [execution workflow README][ex_readme]

[//]: # (References)
  
[out1]: ./example_output_files/output1.bed
[out2]: ./example_output_files/output2.bed
[out3]: ./example_output_files/output3.tsv
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[wiki-tsv]: <https://en.wikipedia.org/wiki/Tab-separated_values>
[ex_readme]: /execution_workflows/README.md

