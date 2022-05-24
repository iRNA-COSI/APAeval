# APAeval Identification: Docker declarations

*OpenEBench/iRNA COSI APAeval benchmarking **Docker declarations** for the **identification** benchmarking event; to be used with the [OEB compatible Nextflow benchmarking workflow][main-nf]*.
___

This README is a brief recap of the three OEB benchmarking workflow steps that are computed in their individual containers. For a general introduction to **APAeval summary workflows** see [the main summary workflow `README.md`][readme-swf]. For a more detailed description of the APAeval **identification** summary workflow see [the identification summary workflow `README.md`][readme-iswf].


>**NOTE for developers:**    
In order to make the workflow containers reproducible and stable in the long-term, make sure to use specific versions in the container base image (e.g.*ubuntu:20.04*, NOT *ubuntu:latest*).

- [Structure](#structure)
- [Usage](#usage)
## Structure
Our benchmarking workflow structure is composed of three docker images / steps:
1. [**Validation**](./i_validation):
the input file format is checked and, if required, the content of the file is validated. The validation generates a
participant dataset. In order to create datasets with structure compatible with the [Elixir Benchmarking Data Model][elixir-data-model], please use the following [python module and JSON schema][json-val].
2. [**Metrics_computation**](./i_metrics):
the predictions are compared with the 'Gold Standards' provided by the community, resulting in a set of quality metrics that reflect the performance of the benchmarking participants. Those metrics are written into assessment datasets. In order to create datasets with structure compatible with the [Elixir Benchmarking Data Model][elixir-data-model], please use the following [python module and JSON schema][json-met].
3. [**Consolidation**](./i_consolidation):
the benchmark itself is performed by merging the assessment metrics of all APAeval participants. The results are provided in both SVG format (scatter plot) and JSON format (aggregation/summary datasets that are compatible with the [Elixir Benchmarking Data Model][elixir-data-model]).


## Usage
After making the necessary changes for your specific challenge, you will have to build the docker image(s) locally by either of the following two methods:

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
If you want to update the docker containers, please remove your original images first:
```
docker image ls #look for the IMAGE_ID of your docker image
docker rmi [IMAGE_ID]
```
Then, you can rebuild the docker image locally (see above).
> NOTE: please don't push to APAeval docker hub unless you created a new version for production that was sanctioned by APAeval admins (ask us on Slack).


[//]: # (References)
[main-nf]: ../main.nf
[readme-swf]: ../../README.md
[readme-iswf]: ../README.md
[elixir-data-model]: https://github.com/inab/benchmarking-data-model
[json-val]: ./i_validation/JSON_templates
[json-met]: ./i_metrics/JSON_templates
