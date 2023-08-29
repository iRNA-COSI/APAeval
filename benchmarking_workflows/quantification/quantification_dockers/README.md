# APAeval Quantification: Docker declarations

*OpenEBench/iRNA COSI APAeval benchmarking **Docker declarations** for the **quantification** benchmarking event; to be used with the [OEB compatible Nextflow benchmarking workflow][main-nf]*.
___

This README is a brief recap of the three OEB benchmarking workflow steps that are computed in their individual containers. For a general introduction to **APAeval benchmarking workflows** see [the main benchmarking workflow `README.md`][readme-bwf]. For a more detailed description of the APAeval **quantification** benchmarking workflow see [the quantification benchmarking workflow `README.md`][readme-qbwf].


>**NOTE for developers:**    
In order to make the workflow containers reproducible and stable in the long-term, make sure to use specific versions in the container base image (e.g.*ubuntu:20.04*, NOT *ubuntu:latest*).

- [Structure](#structure)
- [Usage](#usage)
## Structure
Our benchmarking workflow structure is composed of three docker images / steps:
1. [**Validation**](./q_validation):
The input file format is checked and, if required, the content of the file is validated. The validation generates a
participant dataset. In order to create datasets with structure compatible with the [Elixir Benchmarking Data Model][elixir-data-model], please use the following [python module and JSON schema][oeb-json].

2. [**Metrics_computation**](./q_metrics):
Predictions are compared with the 'Gold Standards' provided by the community, resulting in a set of quality metrics that reflect the performance of the benchmarking participants. Those metrics are written into assessment datasets. In order to create datasets with structure compatible with the [Elixir Benchmarking Data Model][elixir-data-model], please use the following [python module and JSON schema][oeb-json].
Functions for APAeval metrics calculation should be defined in the module [`utils/apaeval`][apa-module] and imported in `compute_metrics.py`


3. [**Consolidation**](./q_consolidation):
Validation results and assessment metrics of all completed challenges are gathered and written to `consolidated_result.json`, which is compatible with the [Elixir Benchmarking Data Model][elixir-data-model]). If present, existing aggregations are imported from a specified directory and included in the consolidation. The metrics to be plotted are collected in "aggregation" objects, which are specified in `aggregation_template.json`.


## Usage
Please check out the sections on [building docker images][build-images] and [running the benchmarking workflow][run-workflow] in the main [APAeval benchmarking workflow README][readme-bwf]


[//]: # (References)
[main-nf]: ../main.nf
[readme-qbwf]: ../README.md
[readme-bwf]: ../../README.md
[apa-module]: ../../../utils/apaeval/src/apaeval/main.py
[build-images]: ../../README.md#7-build-images
[run-workflow]: ../../README.md#8-test-run
[elixir-data-model]: https://github.com/inab/benchmarking-data-model
[oeb-json]: ../../JSON_templates/src/JSON_templates/

