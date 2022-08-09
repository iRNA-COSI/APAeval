## JSON_templates

OEB specific module, script and minimal specification copied from https://github.com/inab/TCGA_benchmarking_dockers.

Excerpt from module doc:

> This module contains functions that generate JSON objects with structure compatible with the Elixir Benchmarking Data Model (https://github.com/inab/benchmarking-data-model). It should be used in the docker declarations to generate the output files in any benchmarking workflow which might be implemented in the OpenEBench infrastructure. Benchmarking workflows architecture can be found in https://github.com/inab/TCGA_benchmarking_workflow
Docker declarations for each step: https://github.com/inab/TCGA_benchmarking_dockers 

### Usage

Locally build the module for debugging: `pip install -e .` into the active APAeval conda environment.
