# OpenEBench assessment workflow in Nextflow

Example pipeline with Nextflow used to assess results, comparing the metrics being computed with this workflow with APAeval pilot benchmark results.

#### DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! 😉
There are three steps in the summary workflow:
 - Validation
    - `input_file`: tab-separated output file from execution workflow
    - Change the validation in `[Y]_dockers/[x]_validation/validation.py` (around line 50) for the kind of input_file requested by your benchmarking event. 
      > NOTE: Currently, the following restrictions are in place:
      >- tab separated file with 6 columns
      >- start and end coordinates (col 2,3) have to be int64
      >- strand (col 6) has to be one of [+,-]
      >- chromosome (col 1) has to be one of [1..22,X,Y]

    - The `validated-participant-data.json` file is not used in the subsequent steps, but the workflow exits if the input file doesn't comply to the specifications of the current benchmarking event
 - Metrics Computation
    - `input_file`: tab-separated output file from execution workflow
    - Change the `[Y]_dockers/[x]_metrics/compute_metrics.py` for the specific metric calculation
    - the gold standard from `metrics_ref_datasets/[challenge].[extension]` and input_file values are used for computing the metrics
      >NOTE: the gold standard file MUST be named in the format `[challenge].[extension]`, `[challenge]` being specified in `challenges_ids` in `nextflow.config`, and `[extension]` being hardcoded within `compute_metrics.py`
    - The `Assessment_datasets.json` file is used in the following step
 - Results Consolidation
    - Inputs the `Assessment_datasets.json` file from the metrics computation step and the `data/` directory, which stores files with benchmark values
    - The current python scripts are as they are in https://github.com/inab/TCGA_benchmarking_dockers, and only supports 2D plots with x and y axes

## Usage
After making the necessary changes for your specific challenge, you will have to build the docker image locally by either of the following two methods:

1. Go to the `[Y]_dockers/` directory and run the following
(note: the `tag_id` should match the one in ../`nextflow.config`)
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

#### Running the summary workflow
One can use the following command to run the summary workflow from command line:
```
nextflow run main.nf -profile docker
```
This reads the parameters from the [nextflow.config](nextflow.config) file.
