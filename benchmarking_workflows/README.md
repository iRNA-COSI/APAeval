# Benchmarking workflows
This is where the benchmarking workflows for APAeval live. The dedicated code for each ***benchmarking event*** resides in its own subdirectory. 

- [Overview](#overview)
- [Benchmarking workflow general description](#benchmarking-workflow-general-description)
  - [Directory structure for benchmarking workflow code](#directory-structure-for-benchmarking-workflow-code)
- [HOW TO: (File) naming requirements](#how-to-file-naming-requirements)
  - [Challenges](#challenges)
  - [Ground truth files](#ground-truth-files)
  - [Participant output (=input) files](#participant-output-input-files)
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

Within a benchmarking event, one or more ***challenges*** will be performed. A **challenge is primarily defined by the ground truth dataset used for performance assessment**. A challenge is evaluated within a ***benchmarking workflow***, which can be run with either docker or singularity, locally or on an HPC infrastructure (currently a profile for Slurm is included). The benchmarking workflow will **compute all metrics relevant for the benchmarking event**. A list of challenge IDs and input files (= output files of one participant for all specified challenges) is passed to the workflow.    

In order to compare the performance of participants within a challenge/event, the respective benchmarking workflow will be run on output files from all eligible participant [method workflows][apaeval-mwfs]. The calculated metrics will be written to `.json` files that can be either submitted to OEB for database storage and online vizualisation, or transformed into a table format that can be used for creating custom plots with the help of scripts from the [APAeval utils directory][apaeval-utils].

## Benchmarking workflow general description
In a first step the provided input files are **validated**. Subsequently, all specified **metrics** are computed, using the matched ground truth files, if applicable. Finally, the **results are gathered** in OEB specific `.json` files per participant.
Based on the created .json files, results can be vizualized on OEB per challenge, such that performance of participants can be compared for each metric.   

In order to eventually be compatible with the OEB infrastructure, benchmarking workflows are written in Nextflow and are structured in a predefined manner that will be described in the following sections.

> DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! 😉

### Directory structure for benchmarking workflow code
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
    |- apaeval/
        |-  src/apaeval/main.py       
```

Within such a directory we find the `main.nf` and `nextflow.config` files, which specify the workflow and all its event-specific parameters, respectively, as well as a `[participant]_[event].config `, which contains the input file- and challenge names for a particular participant. `main.nf` ideally does NOT have to be changed (at least not much) between benchmarking events, as it simply connects the three steps `validation`, `metrics_computation` and `consolidation` inherent to the OEB workflow structure. In contrast, file and tool names have to be adapted in `[participant]_[event].config ` for dedicated workflow runs. The name of `[participant]_[event].config` also has to be specified in `nextflow.config` at the top under `includeConfig`, and finally `nextflow.config` is the place to change additional parameters if necessary.

> ATTENTION: Keep `nextflow.config` unchanged (apart from the above mentioned `includeConfig`) within an event, in order to be able to directly compare the different participant runs.    

Within the benchmarking event's directory resides a subdirectory `specification` with a detailed description of required input and output file formats, as well as of the metrics to be calculated for the respective benchmarking event. The *actual code* is hidden in the directory `[benchmarking_event]_dockers`; For each of the three benchmarking workflow steps required by OEB, a separate docker container will be built:

1. Validation
2. Metrics calculation
3. Consolidation


The "dockers" directories contain Dockerfiles, requirements, and dedicated python scripts. In order to create datasets that are compatible with the [Elixir Benchmarking Data Model][elixir-data-model], the JSON templates in the main `benchmarking_workflows` directory are imported in the respective docker containers. The provided *python scripts*, as well as the module `utils/apaeval` they import, are where the action happens: **These scripts are where you most likely will have to make adjustments for different benchmarking events**. 

## HOW TO: (File) naming requirements
### Challenges
Challenge IDs have to be of the form 
```
[SAMPLE_NAME].([ADDITIONAL_INFO].)[GENOME]
```
where `[SAMPLE_NAME]` is a unique id of the condition represented in the ground truth and assessed by the participant. `[ADDITIONAL_INFO]` is optional; this can be used if several ground truths are obtained from the same condition but differ otherwise, e.g. one is a subset of the other.`[GENOME]` is the genome version used for creating the ground truth. MUST contain either "mm" or "hg", e.g. `mm10` or `hg38_v26`

#### Examples:    
> MmusCortex_adult_R1.TE.mm10    
> GTEXsim_R19.hg38_v26
### Ground truth files
The gold standard file MUST be named in the format 
```
[CHALLENGE].[EXT]
```
where `[CHALLENGE]` is specified in `challenges_ids` in [`[tool]_[event].config`][tool-event-config]. The extension `.bed` is hardcoded within [`compute_metrics.py`][metrics-py] and `[CHALLENGE]` itself has to be of the format described above.   

### Participant output (=input) files
Participant outputs MUST contain the exact `[SAMPLE_NAME]` part of the challenge(s) (see requirements above) they want to participate in. They have to be of the format
```
[PARTICIPANT].[SAMPLE_NAME].[EVENT_ID].[EXT]
```
where `[PARTICIPANT]` is the unique name of the participant to be tested. If a tool is for example run in two different modes, that should be reflected here (like MYTOOL and MYTOOL_SPECIAL). `[EVENT_ID]` is a two-digit-code as follows:   
01 - Identification   
02 - absolute quantification   
03 - differential expression   
04 - relative quantification   

### Metrics
Metric names MUST be **exactly** the same in the respective `compute_metrics.py` and `aggregation_template_X.json` files of a benchmarking workflow. These metric names will then appear in the result `.json` files of the workflow, and will be appearing on OEB plots after uploading the results there.

#### Examples:  
```
Jaccard_index:10nt
percentage_genes_w_correct_nPAS
```


## HOW TO: DEVELOP
For an example of a benchmarking workflow and further instructions, refer to the [quantification benchmarking workflow][q-swf].   
### 1. Copy template
If not done so already, copy the whole contents of the `quantification` directory into the directory for your new benchmarking event. Specify the objectives of your event by adapting the contents of `specification/` .

### 2. Establish proper validation
OEB requires all inputs to be validated. To check for correct input file formats for your benchmarking event, adapt the validation in [`validation.py`][validation-py] (around line 50). Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if necessary.

### 3. Calculate the metrics
Adapt [`compute_metrics.py`][metrics-py] to compare the participant output to the community provided gold standard file(s). You can define custom functions in the [`utils/apaeval`][apa-module] module.

>NOTE: the extension of the gold standard file is currently hardcoded in [`compute_metrics.py`][metrics-py] in line 56. Change this according to your gold standard file format.

Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if applicable.
### 4. Consolidation
The json outputs from the first two steps will be gathered here, and "aggregation objects" for OEB vizualisation will be created based on the `aggregation_template.json`. Thus, this is the file you want to adapt in order to control which metrics are plotted in OEB. You can set visualization types for local plotting in [`manage_assessment_data.py`][assess-py]. The current python scripts have been copied from https://github.com/inab/TCGA_benchmarking_dockers, and only support 2D plots with x and y axes.

### 5. Adapt `nextflow.config` and `main.nf`
In the former you'll have to adjust the docker container names and general workflow parameters, whereas in the latter you'll only have to make changes if you have introduced new workflow parameters (or want to change the wiring of steps, which is not recommended for the sake of attempted OEB compatibility).

### 6. Don't forget to update `README.md`
Describe the type of validation and metric calculation you perform in the `README.md` in your benchmarking event directory (see [example from quantification benchmarking workflow][q-swf]).

### 7. Build images
> ATTENTION: the [apaeval module][apa-module] is installed inside the containers via a git url specified in the respective `requirements.txt` (for [q_validation](quantification/quantification_dockers/q_validation/requirements.txt) and [q_metrics](quantification/quantification_dockers/q_metrics/requirements.txt)). If you made changes to the module, don't forget to push your branch and adjust those urls accordingly.


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


### 8. Test run
One can use the following command to run the quantification benchmarking workflow with the provided test files from command line:
```
nextflow run main.nf -profile docker -c tool_event.config

# Or for running with singularity and slurm:
nextflow run main.nf -profile slurm -c tool_event.config

```
> NOTE: Parameters from the [nextflow.config][nextflow-config] file are read **in addition** to the ones specified with the `-c` flag, but the latter will override any parameters of the same name in the nextflow.config. Don't forget to `includeConfig` the `tool_event.config` in the `nextflow.config`

## HOW TO: "PRODUCTION"
When you have completed the steps described above you can finally run the benchmarking workflow on real data. Below are some hints to help you get going.

### 1. Get input data (= participant output)
Place the participant output into a directory like `DATA/PARTICIPANT_NAME/` and make sure the files are named like `PARTICIPANT_NAME.CHALLENGE_ID.EVENT.EXT`

### 2. Adapt configs
You're going to run the workflow for one participant at a time, but you can specify multiple challenges for that participant. To do so, create a participant specific `[participant]_[event].config` (copy `tool_event.config`). There you'll specify input files and challenge names. ***Don't forget to set your new config's name at the top of `nextflow.config` (directive `includeConfig`).***

### 3. Containers & images
Make sure you have the images appropriate for your system ready. If you're running docker you can use the images you built locally in the [HOW TO: DEVELOP](#7-build-containers) section. If you want to use singularity you'll first have to push those images to a publicly accessible repo, ideally [biocontainers][biocontainers]. Make sure to rename the images (see bash command below) and adjust the paths in the `nextflow.config` accordingly.

```bash
docker tag apaeval/q_consolidation:1.0 your_docker_repo/q_consolidation:1.0
docker push your_docker_repo/q_consolidation:1.0
```

> ATTENTION: always make sure you're using up to date versions of the images. More specifically: DO make sure you have removed old local images and/or cleared singularity caches on your system and checked your `nextflow.config`

## Origin
The APAeval OEB benchmarking workflow is an adaptation of the [TCGA_benchmarking_workflow][tcga-wf] with the [corresponding docker declarations][tcga-docker]. The structure of output files is compatible with the [ELIXIR Benchmarking Data Model][elixir-data-model]. The current version of the workflow is not compatible with the OEB VRE setup, however, only minor changes should be needed to re-establish compatibility.

[//]: # (References)
[apaeval-swfs]: ../images/SWFs.png
[apaeval-mwfs]: ../method_workflows
[apaeval-utils]: https://github.com/iRNA-COSI/APAeval/tree/main/utils
[assess-py]:quantification/quantification_dockers/q_consolidation/manage_assessment_data.py
[biocontainers]: <https://biocontainers-edu.readthedocs.io/en/latest/index.html>
[elixir-data-model]: https://github.com/inab/benchmarking-data-model
[metrics-py]:quantification/quantification_dockers/q_metrics/compute_metrics.py
[apa-module]: https://github.com/iRNA-COSI/APAeval/tree/main/utils/apaeval
[oeb]: <https://openebench.bsc.es/>
[nextflow-config]: quantification/nextflow.config
[q-swf]: quantification/README.md
[tcga-wf]: https://github.com/inab/TCGA_benchmarking_workflow
[tcga-docker]: https://github.com/inab/TCGA_benchmarking_dockers
[tool-event-config]:quantification/tool_event.config
[validation-py]:quantification/quantification_dockers/q_validation/validation.py


