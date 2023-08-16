# Pilot benchmark

## Synopsis

Test benchmark to:

1. Figure out how to specify individual benchmarks, and write corresponding
   method and summary workflows.
2. Write templates for benchmark specification documents, method workflows
   (Nextflow and Snakemake) and summary workflows (only Nextflow).

For this test, the run/execution times and maximum memory usage of two tools
are computed and compared.

## General info

* **Challenge:** Pilot
* **Identifier:** P1

## Inputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | BAM | [Specification][spec-sam-bam] | [Link][in1] |
| 2 | BAM.BAI | [Specification][spec-sam-bam] | [Link][in2] |
| 3 | GTF | [Wikipedia][wiki-gtf] | [Link][in3] |
| 4 | GFF | [Wikipedia][wiki-gff] | [Link][in4] |

### Additional info inputs

N/A

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
| `run_time_sec` | `float` | sec | Total execution time of the tool to be benchmarked, excluding any custom pre-/postprocessing steps |
| `max_mem_mib` | `float` | MiB | Maximum memory used during the execution of the tool to be benchmarked, excluding any custom pre-/postprocessing steps |

## Metrics

| # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Run time | sec | Output 1 | Read file, parse JSON and extract attribute `run_time_sec` | `vector` | N/A |
| 2 | Maximum memory usage | MiB | Output 1 | Read file, parse JSON and extract attribute `max_mem_mib` | `vector` | N/A |

### Additional info metrics

N/A

[//]: # (References)

[in1]: ./example_files/input1.bam
[in2]: ./example_files/input1.bam.bai
[in3]: ./example_files/input3.gtf
[in4]: ./example_files/input4.gff3
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-sam-bam]: <https://samtools.github.io/hts-specs/SAMv1.pdf>
[wiki-gtf]: <https://en.wikipedia.org/wiki/Gene_transfer_format>
[wiki-gff]: <https://en.wikipedia.org/wiki/General_feature_format>
