# OpenEBench assessment workflow in Nextflow

Example pipeline with Nextflow used to assess results, comparing the metrics being computed with this workflow with APAeval pilot benchmark results.

#### DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! 😉
There are three steps in the summary workflow:
 - Validation
    - `input_file`: tab-separated output file from execution workflow
    - Change the `benchmarking_dockers/apaeval_validation/validation.py` for the specific input_file
    - Each input_file may have different fields from different execution workflows
    - `public_ref/[validation_ref].txt` stores the values required to be in the input_files for validating the input_file 
    - The `[output].json` file is not used in the subsequent steps
 - Metrics Computation
    - `input_file`: tab-separated output file from execution workflow
    - Change the `benchmarking_dockers/apaeval_metrics/compute_metrics.py` for the specific input_file and the specific metric calculation
    - the gold standard from `metrics_ref_dataset/[challenge].txt` and input_file values are used for computing the metrics
    - The `output.json` file is used in the following step
 - Results Consolidation
    - Inputs the `output.json` file from the metrics computation step and the `data/` directory, which stores files with benchmark values
    - The current python scripts are as they are in https://github.com/inab/TCGA_benchmarking_dockers, and only supports 2D plots with x and y axes

#### After making the necessary changes for your specific challenge, you will have to build the docker image locally
Go to the specific docker directory for each step in `benchmarking_dockers/`:
 - `apaeval_validation/`, `apaeval_metrics/`, and `apaeval_consolidation/`
and type the following
```
docker build . -t apaeval_[challenge]_[validation/metrics/consolidation]:1.0
```
If you want to update the docker container, please remove your original image first:
```
docker image ls #look for the IMAGE_ID of your docker image
docker rmi [IMAGE_ID]
```
Then, you can rebuild the docker image locally (see above).
 - Note: please don't push it up to docker hub because that may use quite a bit `AWS` rates