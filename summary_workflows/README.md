# Summary workflows
This is where the summary workflows for APAeval live.

- [Overview](#overview)
- [Summary workflow general description](#summary-workflow-general-description)
- [HOW TO](#how-to)
  - [1. Copy template](#1-copy-template)
  - [2. Establish proper validation](#2-establish-proper-validation)
  - [3. Calculate the metrics](#3-calculate-the-metrics)
  - [4. Consolidation](#4-consolidation)
  - [5. Adapt `nextflow.config`](#5-adapt-nextflowconfig)
  - [6. Don't forget to update `README.md`](#6-dont-forget-to-update-readmemd)
  - [7.Build containers](#7build-containers)
  - [8. Run the summary workflow](#8-run-the-summary-workflow)
- [Origin](#origin)
## Overview
![apaeval-swfs][apaeval-swfs]

APAeval consists of **three benchmarking** events to evaluate the performance of different tasks that the methods of interest (=*participants*) might be able to perform: **poly(A) site identification, quantification, and assessment of their differential usage**. A method can participate in one, two or all three events, depending on its functions.   

Within a benchmarking event, one or more *challenges* will be performed. A **challenge is primarily defined by the input dataset used for performance assessment**. A challenge is evaluated within a *summary workflow*, which is run on the OEB infrastructure. The summary workflow will **compute all metrics relevant for the challenge**. Summary workflows can be re-used between challenges, however, depending on the input dataset, different metrics might be calculated, and summary workflows might thus be adapted to individual challenges (Example here: in challenge Ix, metrics I1 and I2 are computed, whereas in challenge Iy, an additional metric I5 is assessed. Apart from calculating metric I5 however, the summary workflows for the challenges Ix and Iy are the same.)    

In order to compare the performance of participants within a challenge, the respective summary workflow will be run on output files from all eligible participant execution workflows (For the participant execution workflows please refer to the [iRNA-COSI/APAeval/execution_workflows][apaeval-ewfs] repository!).   
In a first step the provided files are **validated**. Subsequently, all required **metrics** (scripts In for Identification, Qn for quantification, Dn for differential usage) are computed, using the matched ground truth files, if applicable. Finally, the **results will be gathered in OEB** specific .json files per participant.
Based on the created .json files, OEB will visualize all results per challenge, such that performance of participants can be compared for each metric.

## Summary workflow general description
In order to be compatible with the OEB infrastructure, summary workflows are written in Nextflow and are structured in a predefined manner that will be described in the following sections.

> DON'T FREAK OUT IF YOU'RE UNFAMILIAR WITH `NEXTFLOW`! MOST CHANGES YOU'LL MAKE ARE IN `PYTHON`! ????

The summary workflow code for each *benchmarking event* resides in its own directory (`identification`, `quantification`, `differential_usage`). 


```bash
[benchmarking_event]/
    |- main.nf
    |- nextflow.config
    |- specification/
        |- example_files/
        |- specification.md
    |- [benchmarking_event]_dockers/
        |- validation/
            |- JSON_templates/
            |- Dockerfile
            |- requirements.txt
            |- validation.py
            |- ...
        |- metrics/
            |- ...
            |- compute_metrics.py
            |- [python_module]/
            |- ...
        |- consolidation/
            |- ...
            |- manage_assessment_data.py
            |- ...            
```

Within such a directory we find the `main.nf` and `nextflow.config` files, which specify the workflow and all its parameters, respectively. `main.nf` ideally does NOT have to be changed between benchmarking events, as it simply connects the three steps `validation`, `metrics_computation` and `consolidation` inherent to the OEB infrastructure. In contrast, file and tool names, and parameters, have to be adapted in `nextflow.config` for *each individual tool and challenge*.   
Within the benchmarking event's directory resides a subdirectory `specification` with a detailed description of required input and output file formats, as well as of the metrics to be calculated for the respective benchmarking event. The *actual code* is hidden in the directory `[benchmarking_event]_dockers`; For each of the three summary workflow steps required by OEB, a separate docker container will be built:

1. Validation
2. Metrics calculation
3. Consolidation


The "dockers" directories contain Dockerfiles, requirements, and OEB `.json` templates. In order to create datasets that are compatible with the [Elixir Benchmarking Data Model][elixir-data-model], **the provided JSON templates have to be used**. The provided *python scripts* are where the action happens: **These scripts are where you most likely will have to make adjustments for different benchmarking events**. 


## HOW TO
For an example of a summary workflow and further instructions, refer to the [quantification summary workflow][q-swf].   
### 1. Copy template
If not done so already, copy the whole contents of the `quantification` directory into the respective benchmarking event directory (`identification` or `differential_usage`). Adapt the contents of `specification/`.

### 2. Establish proper validation
OEB requires all inputs to be validated. To check for correct input file formats for your benchmarking event, adapt the validation in [`validation.py`][validation-py] (around line 50). Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if applicable.

### 3. Calculate the metrics
Adapt [`compute_metrics.py`][metrics-py] to compare the participant output to the community provided gold standard file(s). You can define custom functions in modules within the `metrics` docker directory (in the example of the quantification summary workflow, we created a module [`matchPAS`][matchpas])

>NOTE: the extension of the gold standard file is currently hardcoded in [`compute_metrics.py`][metrics-py] in line 50. Change this according to your gold standard file format.

Update the corresponding `requirements.txt`, `constraints.txt` and `Dockerfile` for installation of additional packages, if applicable.
### 4. Consolidation
You can set visualization types in [`manage_assessment_data.py`][assess-py]. The current python scripts have been copied from https://github.com/inab/TCGA_benchmarking_dockers, and only support 2D plots with x and y axes.

### 5. Adapt `nextflow.config`
Here you'll adjust the docker container names, as well as input, participant and challenge names for *each tool and each challenge*.

### 6. Don't forget to update `README.md`
Describe the type of validation and metric calculation you perform in the `README.md` in your benchmarking event directory (see [example from quantification summary workflow][q-swf]).

### 7.Build containers
After making the necessary changes for your specific challenge, you will have to build the docker images locally by either of the following two methods:

1. Go to the `[Y]_dockers/` directory and run the following
(note: the `tag_id` should match the one in `nextflow.config`)
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
> NOTE: please don't push to APAeval docker hub unless sanctioned by APAeval admins (ask us on Slack)

### 8. Run the summary workflow
One can use the following command to run the summary workflow from command line:
```
nextflow run main.nf -profile docker
```
This reads the parameters from the [nextflow.config][nextflow-config] file.

## Origin
The APAeval OEB summary workflow is an adaptation of the [TCGA_benchmarking_workflow][tcga-wf] with the [corresponding docker declarations][tcga-docker]. The output structure is compatible with the [ELIXIR Benchmarking Data Model][elixir-data-model].

[//]: # (References)
[apaeval-swfs]: ../images/SWFs.png
[apaeval-ewfs]: https://github.com/iRNA-COSI/APAeval/tree/main/execution_workflows
[elixir-data-model]: https://github.com/inab/benchmarking-data-model
[q-swf]: quantification/README.md
[validation-py]:quantification/quantification_dockers/q_validation/validation.py
[metrics-py]:quantification/quantification_dockers/q_metrics/compute_metrics.py
[matchpas]: quantification/quantification_dockers/q_metrics/matchPAS
[assess-py]:quantification/quantification_dockers/q_consolidation/manage_assessment_data.py
[nextflow-config]: quantification/nextflow.config
[tcga-wf]: https://github.com/inab/TCGA_benchmarking_workflow
[tcga-docker]: https://github.com/inab/TCGA_benchmarking_dockers


