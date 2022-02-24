# APAeval benchmarking Docker declarations

OpenEBench/iRNA COSI APAeval benchmarking Docker declarations, which define the architecture of benchmarking workflows to be implemented in OpenEBench.

**NOTE for developers.** In order to make the workflow containers reproducible and stable in the long-term, make sure to use specific versions in the container base image (e.g.*ubuntu:16.04*, NOT *ubuntu:latest*).

## Structure
Our benchmarking workflow structure is composed by three docker images / steps:
1. [**Validation**](./q_validation):
the input file format is checked and, if required, the content of the file is validated. The validation generates a
participant dataset. In order to create datasets with structure compatible with the [Elixir
    Benchmarking Data Model](https://github.com/inab/benchmarking-data-model), please use the following [python module and JSON schema](./apaeval_validation/JSON_templates)
2. [**Metrics_computation**](./q_metrics):
the predictions are compared with the 'Gold Standards' provided by the community, resulting in a set of quality metrics that reflect the performance of the benchmarking participants. Those metrics are written into assessment datasets. In order to create datasets with structure compatible with the [Elixir
    Benchmarking Data Model](https://github.com/inab/benchmarking-data-model), please use the following [python module and JSON schema](./apaeval_metrics/JSON_templates)
3. [**Consolidation**](./q_consolidation):
the benchmark itself is performed by merging the assessment metrics of all APAeval participants. The results are provided in both SVG format (scatter plot) and JSON format (aggregation/summary datasets that are compatible with the [Elixir Benchmarking Data Model](https://github.com/inab/benchmarking-data-model))

For the participant execution workflows please refer to the [iRNA-COSI/APAeval/execution_workflows](https://github.com/iRNA-COSI/APAeval/tree/main/execution_workflows) repository.

## Usage
After making the necessary changes for your specific challenge, you will have to build the docker image locally by either of the following two methods:

1. Go to the `[Y]_dockers/` directory and run the following (note: the `tag_id` should match the one in ../`nextflow.config`):
```
run `./build.sh <tag_id>`
```
2. Go to the specific docker directory for each step in `[Y]_dockers/`:
 - `[X]_validation/`, `[X]_metrics/`, and `[X]_consolidation/`
and run the following
```
docker build . -t apaeval/[X]_[validation/metrics/consolidation]:1.0
```
If you want to update the docker container, please remove your original image first:
```
docker image ls #look for the IMAGE_ID of your docker image
docker rmi [IMAGE_ID]
```
Then, you can rebuild the docker image locally (see above).
 - Note: please don't push it up to docker hub because that may use quite a bit `AWS` rates

