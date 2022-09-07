# Summary workflows
This is where the summary workflows for APAeval live. The dedicated code for each ***benchmarking event*** resides in its own subdirectory. 

- [Overview](#overview)
- [Summary workflow general description](#summary-workflow-general-description)
  - [Directory structure for summary workflow code](#directory-structure-for-summary-workflow-code)
- [HOW TO: DEVELOP](#how-to-develop)
  - [1. Copy template](#1-copy-template)
  - [2. Establish proper validation](#2-establish-proper-validation)
  - [3. Calculate the metrics](#3-calculate-the-metrics)
  - [4. Consolidation](#4-consolidation)
  - [5. Adapt `nextflow.config` and `main.nf`](#5-adapt-nextflowconfig-and-mainnf)
  - [6. Don't forget to update `README.md`](#6-dont-forget-to-update-readmemd)
  - [7. Build images](#7-build-images)
  - [8. Test run](#8-test-run)
- [HOW TO: "PRODUCTION"](#how-to-production)
  - [1. Get input data (= participant output)](#1-get-input-data--participant-output)
  - [2. Adapt configs](#2-adapt-configs)
  - [3. Containers & images](#3-containers--images)
- [Origin](#origin)
## Overview
![apaeval-swfs][apaeval-swfs]

APAeval consists of a number of ***benchmarking events*** to evaluate the performance of different tasks that the methods of interest (=***participants***) might be able to perform: **poly(A) site identification, absolute quantification, relative quantification and assessment of their differential usage**. A method can participate in one or several events, depending on its functions.   

Within a benchmarking event, one or more ***challenges*** will be performed. A **challenge is primarily defined by the input dataset used for performance assessment**. A challenge is evaluated within a ***summary workflow***, which can be run with either docker or singularity, locally or on an HPC infrastructure (currently a profile for Slurm is included). The summary workflow will **compute all metrics relevant for the benchmarking event**. A list of challenge IDs and input files (= output files of one participant for all specified challenges) is passed to the workflow.    

In order to compare the performance of participants within a challenge/event, the respective summary workflow will be run on output files from all eligible participant execution workflows (For the participant execution workflows please refer to the [iRNA-COSI/APAeval/execution_workflows][apaeval-ewfs] repository!). The calculated metrics will be written to `.json` files that can be either submitted to OEB for database storage and online vizualisation, or transformed into a table format that can be used for creating custom plots with the help of scripts from the [APAeval utils directory][apaeval-utils].

## Summary workflow general description
In a first step the provided input files are **validated**. Subsequently, all specified **metrics** are computed, using the matched ground truth files, if applicable. Finally, the **results are gathered** in OEB specific `.json` files per participant.
Based on the created .json files, results can be vizualized on OEB per challenge, such that performance of participants can be compared for each metric.   

In order to eventually be compatible with the OEB infrastructure, summary workflows are written in Nextflow and are structured in a predefined manner that will be described in the following sections.

> DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! 😉

### Directory structure for summary workflow code
```bash
summmary_workflows/
    |- JSON_templates/
    |- [benchmarking_event]/
        |- main.nf
        |- nextflow.config
        |- [participant]_[event].config
        |- specification/
            |- example_files/
            |- specification.md
        |- [benchmarking_event]_dockers/
            |- validation/
                |- Dockerfile
                |- requirements.txt
                |- validation.py
                |- ...
            |- metrics/
                |- Dockerfile
                |- requirements.txt
                |- compute_metrics.py
                |- ...
            |- consolidation/
                |- Dockerfile
                |- requirements.txt
                |- aggregation.py
                |- merge_data_model_files.py
                |- ...    
...
utils/
    |- matchPAS/
        |-  src/matchPAS/main.py       
```

Within such a directory we find the `main.nf` and `nextflow.config` files, which specify the workflow and all its event-specific parameters, respectively, as well as a `[participant]_[event].config `, which contains the input file- and challenge names for a particular participant. `main.nf` ideally does NOT have to be changed (at least not much) between benchmarking events, as it simply connects the three steps `validation`, `metrics_computation` and `consolidation` inherent to the OEB workflow structure. In contrast, file and tool names have to be adapted in `[participant]_[event].config ` for dedicated workflow runs. The name of `[participant]_[event].config` also has to be specified in `nextflow.config` at the top under `includeConfig`, and finally `nextflow.config` is the place to change additional parameters if necessary.

> ATTENTION: Keep `nextflow.config` unchanged (apart from the above mentioned `includeConfig`) within an event, in order to be able to directly compare the different participant runs.    

Within the benchmarking event's directory resides a subdirectory `specification` with a detailed description of required input and output file formats, as well as of the metrics to be calculated for the respective benchmarking event. The *actual code* is hidden in the directory `[benchmarking_event]_dockers`; For each of the three summary workflow steps required by OEB, a separate docker container will be built:

1. Validation
2. Metrics calculation
3. Consolidation


The "dockers" directories contain Dockerfiles, requirements, and dedicated python scripts. In order to create datasets that are compatible with the [Elixir Benchmarking Data Model][elixir-data-model], the JSON templates in the main `summary_workflows` directory are imported in the respective docker containers. The provided *python scripts*, as well as the module `utils/matchPAS` they import, are where the action happens: **These scripts are where you most likely will have to make adjustments for different benchmarking events**. 


## HOW TO: DEVELOP
For an example of a summary workflow and further instructions, refer to the [quantification summary workflow][q-swf].   
### 1. Copy template
If not done so already, copy the whole contents of the `quantification` directory into the directory for your new benchmarking event. Specify the objectives of your event by adapting the contents of `specification/` .

### 2. Establish proper validation
OEB requires all inputs to be validated. To check for correct input file formats for your benchmarking event, adapt the validation in [`validation.py`][validation-py] (around line 50). Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if necessary.

### 3. Calculate the metrics
Adapt [`compute_metrics.py`][metrics-py] to compare the participant output to the community provided gold standard file(s). You can define custom functions in modules within the [`utils` directory][apaeval-utils] (in the example of the quantification summary workflow, we created a module [`matchPAS`][matchpas]).

>NOTE: the extension of the gold standard file is currently hardcoded in [`compute_metrics.py`][metrics-py] in line 56. Change this according to your gold standard file format.

Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if applicable.
### 4. Consolidation
The json outputs from the first two steps will be gathered here, and "aggregation objects" for OEB vizualisation will be created based on the `aggregation_template.json`. Thus, this is the file you want to adapt in order to control which metrics are plotted in OEB. You can set visualization types for local plotting in [`manage_assessment_data.py`][assess-py]. The current python scripts have been copied from https://github.com/inab/TCGA_benchmarking_dockers, and only support 2D plots with x and y axes.

### 5. Adapt `nextflow.config` and `main.nf`
In the former you'll have to adjust the docker container names and general workflow parameters, whereas in the latter you'll only have to make changes if you have introduced new workflow parameters (or want to change the wiring of steps, which is not recommended for the sake of attempted OEB compatibility).

### 6. Don't forget to update `README.md`
Describe the type of validation and metric calculation you perform in the `README.md` in your benchmarking event directory (see [example from quantification summary workflow][q-swf]).

### 7. Build images
After making the necessary changes for your specific event, you will have to build the docker images locally by either of the following two methods:

1. Go to the `[X]_dockers/` directory and run the following
(note: the `tag_id` should match the one in your `nextflow.config`)
```
run `./build.sh <tag_id>`
```
2. Go to the specific docker directory for each step in `[X]_dockers/`:
 - `[X]_validation/`, `[X]_metrics/`, or `[X]_consolidation/`
and run the following
```
docker build . -t apaeval/[X]_[validation OR metrics OR consolidation]:<tag_id>
```
If you want to update the docker containers, please remove your original images first:
```
docker image ls #look for the IMAGE_ID of your docker image
docker rmi [IMAGE_ID]
```
Then, you can rebuild the docker image locally (see above).
> NOTE: please don't push to APAeval docker hub unless sanctioned by APAeval admins (ask us on Slack)

### 8. Test run
One can use the following command to run the quantification summary workflow with the provided test files from command line:
```
nextflow run main.nf -profile docker -c tool_event.config

# Or for running with singularity and slurm:
nextflow run main.nf -profile slurm -c tool_event.config

```
> NOTE: Parameters from the [nextflow.config][nextflow-config] file are read **in addition** to the ones specified with the `-c` flag, but the latter will override any parameters of the same name in the nextflow.config. Don't forget to `includeConfig` the `tool_event.config` in the `nextflow.config`

## HOW TO: "PRODUCTION"
When you have completed the steps described above you can finally run the summary workflow on real data. Below are some hints to help you get going.

### 1. Get input data (= participant output)
Place the participant output into a directory like `DATA/PARTICIPANT_NAME/` and make sure the files are named like `PARTICIPANT_NAME.CHALLENGE_ID.EVENT.EXT`

### 2. Adapt configs
You're going to run the workflow for one participant at a time, but you can specify multiple challenges for that participant. To do so, create a participant specific `[participant]_[event].config` (copy `tool_event.config`). There you'll specify input files and challenge names. ***Don't forget to set your new config's name at the top of `nextflow.config` (directive `includeConfig`).***

### 3. Containers & images
Make sure you have the images appropriate for your system ready. If you're running docker you can use the images you built locally in the [HOW TO: DEVELOP](#7-build-containers) section. If you want to use singularity you'll first have to push those images to a publicly accessible repo, e.g. dockerhub. Make sure to rename the images (see bash command below) and adjust the paths in the `nextflow.config` accordingly.

```bash
docker tag apaeval/q_consolidation:1.0 your_docker_repo/q_consolidation:1.0
docker push your_docker_repo/q_consolidation:1.0
```

> ATTENTION: always make sure you're using up to date versions of the images. More specifically: DO make sure you have removed old local images and/or cleared singularity caches on your system and checked your `nextflow.config`

## Origin
The APAeval OEB summary workflow is an adaptation of the [TCGA_benchmarking_workflow][tcga-wf] with the [corresponding docker declarations][tcga-docker]. The structure of output files is compatible with the [ELIXIR Benchmarking Data Model][elixir-data-model]. The current version of the workflow is not compatible with the OEB VRE setup, however, only minor changes should be needed to re-establish compatibility.

[//]: # (References)
[apaeval-swfs]: ../images/SWFs.png
[apaeval-ewfs]: https://github.com/iRNA-COSI/APAeval/tree/main/execution_workflows
[apaeval-utils]: https://github.com/iRNA-COSI/APAeval/tree/main/utils
[assess-py]:quantification/quantification_dockers/q_consolidation/manage_assessment_data.py
[elixir-data-model]: https://github.com/inab/benchmarking-data-model
[metrics-py]:quantification/quantification_dockers/q_metrics/compute_metrics.py
[matchpas]: https://github.com/iRNA-COSI/APAeval/tree/main/utils/matchPAS
[oeb]: <https://openebench.bsc.es/>
[nextflow-config]: quantification/nextflow.config
[q-swf]: quantification/README.md
[tcga-wf]: https://github.com/inab/TCGA_benchmarking_workflow
[tcga-docker]: https://github.com/inab/TCGA_benchmarking_dockers
[validation-py]:quantification/quantification_dockers/q_validation/validation.py


