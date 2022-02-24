# OpenEBench assessment workflow in Nextflow

Example pipeline with Nextflow used to assess results, comparing the metrics being computed with this workflow with APAeval pilot benchmark results.

#### DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! 😉
There are three steps in the summary workflow:
 - Validation
    - `input_file`: tab-separated output file from execution workflow
    - Change the validation in `benchmarking_dockers/apaeval_validation/validation.py` (around line 50)for the kind of input_file requested by your benchmarking event
    - The `[output].json` file is not used in the subsequent steps, but the workflow exits if the input file doesn't comply to the specifications of the current benchmarking event
 - Metrics Computation
    - `input_file`: tab-separated output file from execution workflow
    - Change the `benchmarking_dockers/apaeval_metrics/compute_metrics.py` for the specific input_file and the specific metric calculation
    - the gold standard from `metrics_ref_dataset/[challenge].txt` and input_file values are used for computing the metrics
    - The `output.json` file is used in the following step
 - Results Consolidation
    - Inputs the `output.json` file from the metrics computation step and the `data/` directory, which stores files with benchmark values
    - The current python scripts are as they are in https://github.com/inab/TCGA_benchmarking_dockers, and only supports 2D plots with x and y axes

#### After making the necessary changes for your specific challenge, you will have to build the docker image locally by either of the following two methods:
1. Go to the `benchmarking_dockers/` directory and run the following
```
run `./build.sh <tag_id>`
```
2. Go to the specific docker directory for each step in `benchmarking_dockers/`:
 - `apaeval_validation/`, `apaeval_metrics/`, and `apaeval_consolidation/`
and run the following
```
docker build . -t apaeval/[challenge]_[validation/metrics/consolidation]:1.0
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
