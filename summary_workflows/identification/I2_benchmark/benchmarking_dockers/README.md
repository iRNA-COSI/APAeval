# TCGA benchmarking Docker declarations

OpenEBench TCGA benchmarking Docker declarations, which define the architecture of benchmarking workflows to be implemented in OpenEBench.

**NOTE for developers.** In order to make the workflow containers reproducible and stable in the long-term, make sure to use specific versions in the container base image (e.g.*ubuntu:16.04*, NOT *ubuntu:latest*).

## Structure
Our benchmarking workflow structure is composed by three docker images / steps:
1. [**Validation**](./tcga_validation):
the input file format is checked and, if required, the content of the file is validated. The validation generates a
participant dataset. In order to create datasets with structure compatible with the [Elixir
    Benchmarking Data Model](https://github.com/inab/benchmarking-data-model), please use the following [python module and JSON schema](./tcga_validation/JSON_templates)
2. [**Metrics_computation**](./tcga_metrics):
the predictions are compared with the 'Gold Standards' provided by the community, which, in this case, results in two
performance metrics - precision (Positive Predictive Value) and recall(True Positive Rate). Those metrics are written
into assessment datasets. In order to create datasets with structure compatible with the [Elixir
    Benchmarking Data Model](https://github.com/inab/benchmarking-data-model), please use the following [python module and JSON schema](./tcga_metrics/JSON_templates)
3. [**Consolidation**](./tcga_consolidation):
the benchmark itself is performed by merging the assessment metrics with the rest of TCGA data. The results are provided
SVG format - scatter plot, and JSON format - aggregation/summary datasets, which are also compatible with the [Elixir
    Benchmarking Data Model](https://github.com/inab/benchmarking-data-model).

Find more information about the TCGA Cancer Drivers Pipeline [here](https://github.com/inab/TCGA_benchmarking_workflow).

## TCGA sample files
* [example_input.txt](./example_input.txt) is a gene predictions file which can be used as input to test the containers.
* [TCGA_full_data](./TCGA_full_data) folder contains all the reference data required by the containers. It is derived from the manuscript:
[Comprehensive Characterization of Cancer Driver Genes and Mutations](https://www.cell.com/cell/fulltext/S0092-8674%2818%2930237-X?code=cell-site), Bailey et al, 2018, Cell [![doi:10.1016/j.cell.2018.02.060](https://img.shields.io/badge/doi-10.1016%2Fj.cell.2018.02.060-green.svg)](https://doi.org/10.1016/j.cell.2018.02.060)
* [sample_results](./sample_results) folder contains an example output for each of the containers run, separated in subfolders, with two cancer types / challenges selected (ACC, BLCA). Results found in [sample_results/consolidation_results](./sample_results/consolidation_results) can be visualized in the browser using [this javascript library](https://github.com/inab/benchmarking_workflows_results_visualizer).


## Usage
In order to build the Docker images locally, please run `./build.sh`
