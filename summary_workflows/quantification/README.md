# APAeval Quantification

*OpenEBench compatible Nextflow workflow for assessment of a bioinformatics tool's performance in **quantifying poly(A) site usage** from RNA-seq data*
___

This README describes the APAeval **(absolute) quantification** summary workflow. For a more general introduction to summary workflows see [the main summary workflow `README.md`][readme-swf]. For the specification of metrics, in- and output file formats, see [the quantification benchmarks specification][q-spec].

- [Description of steps](#description-of-steps)
  - [1. Validation](#1-validation)
  - [2. Metrics Computation](#2-metrics-computation)
  - [3. Results Consolidation](#3-results-consolidation)
- [Usage](#usage)
## Description of steps
### 1. Validation
- `input_file`: output file from execution workflow in bed6 format
- Validation checks performed in [`quantification_dockers/q_validation/validation.py`][validation-py]:
   - input file has to be tab separated file with 6 columns
   - start and end coordinates (col 2,3) have to be int64
   - strand (col 6) has to be one of [+,-]
   - chromosome (col 1) has to match the ones from the genome annotation (see below `genome_dir`)
  
- The `validated_[participant].[challenge].[event].json` file is used in the consolidation step, but not in the compute metrics one. However, the workflow exits after the validation step if *one or more* of the input files don't comply to the specifications of the current benchmarking event
  
### 2. Metrics Computation
- "input file" and "gold standard file" will be compared in order to calculate the metrics
- `input_file`: output file from execution workflow in bed6 format
- `gold standard`: bed6 file derived from 3'end sequencing on the same sample(s) as the RNA-seq data used in the challenge
>NOTE: the gold standard file MUST be named in the format `[challenge].bed`, where `[challenge]` is specified in `challenges_ids` in [`[tool]_[event].config`][tool-event-config]. The extension `.bed` is hardcoded within [`compute_metrics.py`][metrics-py].
- `windows` parameter is used to compute metrics for a list of window sizes.
    - For running on OEB: the parameter is read from `nextflow.config`.
- `genome_dir`: Directory to genome annotation in gtf format with 9 fields as specified [here](https://www.gencodegenes.org/pages/data_format.html). The gtf is used for the relative PAS usage metric computation.
  - For running on OEB: The genome directory is specified in `nextflow.config`
  - For the test data, challenge `challenge_1.mm10` with ground truth file `challenge_1.mm10.bed` will use genome file `gencode.test.mm10.gtf`, because both contain `mm10` within two dots in the filename.
> NOTE: the genome file needs to contain the same substring as the challenge. That is, challenge `[partone].[organism].[partwo].bed` requires a genome annotation file like `[partone].[organism].[partwo].gtf`, where `[organism]` starts with *mm* or *hg* (only these two currently supported). And `[partone]` and `[parttwo]` can be an aribitrary string (or empty string).
- APAeval custom functions called in [`quantification_dockers/q_metrics/compute_metrics.py`][metrics-py] are defined in `utils/apaeval`
- The `assessments_[participant].[challenge].[event].json` file is used in the consolidation step


### 3. Results Consolidation
- Gathers *all* `validated_[participant].[challenge].[event].json` files from the validation step, *all* `assessments_[participant].[challenge].[event].json` files from the metrics computation step, and - if available - existing aggregation data (currently imported from the `data/` directory; in `nextflow.config`: `aggregation_dir`)
- Outputs OEB compatible `consolidated_result.json` file for the tested participant
- "aggregation" objects in the `consolidated_result.json` determine which metrics are to be plotted against each other on the OEB website
- In order to specify which of the metrics present in the assessment objects should be plotted on OEB, the file `quantification_dockers/q_consolidation/aggregation_template.json` has to be modified.

## Usage
Please check out the sections on [building docker images][build-images] and [running the summary workflow][run-workflow] in the main [APAeval summary workflow README][readme-swf]



[//]: # (References)
[readme-swf]: ../README.md
[build-images]: ../README.md#7-build-images
[run-workflow]: ../README.md#8-test-run
[q-spec]: ./specification/
[validation-py]:./quantification_dockers/q_validation/validation.py
[metrics-py]:./quantification_dockers/q_metrics/compute_metrics.py
[nextflow-config]: ./nextflow.config
[tool-event-config]: ./tool_event.config
