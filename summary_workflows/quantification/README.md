# APAeval Quantification
*OpenEBench compatible Nextflow workflow for assessment of a bioinformatics tool's performance in **quantifying poly(A) site usage** from RNA-seq data*
___

This README describes the APAeval **quantification** summary workflow. For a more general introduction to summary workflows see [the main summary workflow `README.md`][readme-swf]

- [Description of steps](#description-of-steps)
  - [1. Validation](#1-validation)
  - [2. Metrics Computation](#2-metrics-computation)
  - [3. Results Consolidation](#3-results-consolidation)
- [Usage](#usage)
  - [Building docker images](#building-docker-images)
  - [Running the summary workflow](#running-the-summary-workflow)
## Description of steps
### 1. Validation
- `input_file`: output file from execution workflow in bed6 format
- Validation checks performed in `quantification_dockers/q_validation/validation.py`:
   - input file has to be tab separated file with 6 columns
   - start and end coordinates (col 2,3) have to be int64
   - strand (col 6) has to be one of [+,-]
   - chromosome (col 1) has to be one of [1..22,X,Y]   
  
- The `validated-participant-data.json` file is not used in the subsequent steps, but the workflow exits if the input file doesn't comply to the specifications of the current benchmarking event
  
### 2. Metrics Computation
- "input file" and "gold standard file" will be compared in order to calculate the metrics
- `input_file`: output file from execution workflow in bed6 format
- `gold standard`: bed6 file derived from 3'end sequencing on the same sample(s) as the RNA-seq data used in the challenge
>NOTE: the gold standard file MUST be named in the format `[challenge].bed`, where `[challenge]` is specified in `challenges_ids` in `nextflow.config`. The extension `.bed` is hardcoded within `compute_metrics.py`
- APAeval custom functions called in `quantification_dockers/q_metrics/compute_metrics.py` are defined in `quantification_dockers/q_metrics/matchPAS`
- The `Assessment_datasets.json` file is used in the following step

### 3. Results Consolidation
- Gathers the `Assessment_datasets.json` file from the metrics computation step and existing assessment data (currently imported from the `data/` directory; in `nextflow.config`: `assess_dir`)
- Python scripts have been copied from https://github.com/inab/TCGA_benchmarking_dockers, and only support 2D plots with x and y axes

## Usage
### Building docker images

1. Go to the `[Y]_dockers/` directory and run the following
(note: the `tag_id` should match the one in ../`nextflow.config`)
```
run `./build.sh <tag_id>`
```
OR   
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
> NOTE: please don't push to APAeval docker hub unless sanctioned by APAeval admins (ask us on Slack)

### Running the summary workflow
One can use the following command to run the summary workflow from command line:
```
nextflow run main.nf -profile docker
```
This reads the parameters from the [nextflow.config](nextflow.config) file.


[//]: # (References)
[readme-swf]: ../README.md