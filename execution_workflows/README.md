# Execution workflows
This is where the execution workflows for APAeval live. 

## Overview
_Execution workflows_ contain all steps that need to be run _per method_ (in OEB terms: per _participant_). Depending on the participant, an execution workflow will have to perform **pre-processing** steps to convert the APAeval sanctioned input files into a format that the participant can consume. This does not include e.g. adapter trimming or mapping of reads, as those steps are already performed in our general pre-processing pipeline. After pre-processing, the actual **execution of the method** has to be implemented, and subsequently **post-processing** steps might be required to convert the obtained output into the format defined by the APAeval specifications.

![EWFs][apaeval-ewfs]

## More details
1. **Sanctioned input files**: Each of the processed input data APAeval is using for their challenges is provided as .bam file (For ease of use, as some participants require .fastq input, we also provide a .fastq file reconstructed from the processed .bam files). If participants need other file formats, these HAVE TO be created as part of the pre-processing **within** execution workflows (see 2.). Similarly, for each dataset we provide a gencode annotation in .gtf format, as well as a reference PAS atlas in .bed format for participants that depend on pre-defined PAS. All other annotation formats that might be needed HAVE TO be created from those. Non-sanctioned annotation- or similar auxiliary files MUST NOT be downloaded as part of the execution workflows, in order to ensure comparability of all participants’ performance.  

> As several execution workflows might have to do the same pre-processing tasks, we created a [utils][utils] directory, where scripts (which have their corresponding docker images uploaded to the APAeval dockerhub) are stored. Please check the [utils][utils] directory before writing your own conversion scripts, and/or add your pre-processing scripts to the utils directory if you think others might be able to re-use them.



2. **Method execution:** For each method to be benchmarked (“participant”) one execution workflow has to be written. The workflow MUST include all necessary pre- and post-processing steps that are needed to get from the input formats provided by APAeval (see 1.), to the output specified by APAeval in their metrics specifications (see 3.). Those specifications might require that more than one output file is created by the workflow. If a method has distinct run modes, the calls to those should be parameterized, so that the workflow can easily be run again with a different mode. If those runmodes significantly alter the behaviour of the method, please discuss with the APAeval community whether the distinct runmodes should actually be treated as distinct participants in APAeval (see [section on parameters](#parameters)).   
In general, all relevant participant parameters should be configurable in the workflow config files. Parameters, file names, run modes, etc. MUST NOT be hardcoded within the workflow. 

> IMPORTANT: Do not download any other annotation files because the docs of your participant say so. Instead, create all files the participant needs from the ones provided by APAeval. If you don't know how, please don't hesitate to start discussions within the APAeval community! Chances are high that somebody already encountered a similar problem and will be able to help.

3. **Post-processing:** To ensure compatibility with the OEB benchmarking events, [specifications for file formats][spec-doc] (output of execution workflows = input for summary workflows) are provided by APAeval. There is one specification per metric (=statistical parameter to assess performance of a participant), but calculation of several metrics can require a common input file format (thus, the file has to be created only once by the execution workflow). The three required execution workflow outputs are a) a bed file containing coordinates of identified PAS and their respective expression in tpm, b) a tsv file containing information on differential expression (if applicable), c) a json file containing information about compute resource and time requirements (see specifications for detailed description of the file formats). These files have to be created within the execution workflows as post-processing steps.


_Execution workflows_ should be implemented in either [Nexflow][nf] or
[Snakemake][snakemake], and individual steps should be isolated through the
use of either [Conda][conda] virtual environments (deprecated; to run on AWS we need containerized workflows) or
[Docker][docker]/[Singularity][singularity] containers. For more information on how to create these containers, see section [containers](#containers).

## Templates
To implement an execution workflow for a participant, copy either the [snakemake template][snakemake-template] 
or the [nextflow template dsl1][nextflow-template-dsl1]/[nextflow template dsl2][nextflow-template-dsl2] into 
the participant's directory and adapt the workflow directory names as described in the template's `README`. 
Don't forget to adapt the `README` itself as well.  


Example:    
```bash
execution_workflows/
 |--QAPA/
     |--QAPA_snakemake/
          |--workflow/Snakefile
          |--config/config.QAPA.yaml
          |--envs/QAPA.yaml
          |--envs/QAPA.Dockerfile
          |-- ...
 |--MISO/
      |--MISO_nextflow/
          |-- ...
```

## Containers
For the sake of reproducibility and interoperability, we require the use of docker containers in our execution workflows. The participants to be benchmarked have to be available in a container, but also any other tools that are used for pre- or post-processing in an execution workflow should be containerized. Whether you get individual containers for all the tools of your workflow, or combine them inside one container is up to you (The former being the more flexible option of course).

> IMPORTANT: Do check out the [utils directory][utils] before you work on containers for pre- or post-processing tools, maybe someone already did the same thing. If not, and you're gonna build useful containers, don't forget to add them there as well.

Here are some pointers on how to best approach the containerization:

1. Check if your participant (or other tool) is already available as a Docker container, e.g. at
   - [dockerhub][dockerhub]
   - [biocontainers][biocontainers]
   - google for `[TOOL_NAME] Docker`  or `[TOOL_NAME] Dockerfile`


2. If no Docker image is availabe for your tool
   - build one from a Dockerfile. There are thousands of Docker tutorials, for example [here][docker-tutorial]. An example of a Dockerfile in APAeval you can find [here][dockerfile].
   - once you've successfully built your image locally, you should push the corresponding (*tested*) `Dockerfile` to your branch on the APAeval repo. In order to then get your image to our dockerhub account, please make yourself heard in our slack space and [Alex][docker-contact-alex] or [Yuk Kei][docker-contact-yukkei] will help you.
   - naming conventions: 
        - if your container only contains one tool:
        `apaeval/{tool_name}:{tool_version}`, e.g. `apaeval/my_tool:v1.0.0`
        - if you combine all tools required for your workflow: `apaeval/exwf_{participant_name}:{commit_hash}`, where `commit_hash` is the short SHA of the Git commit in the APAeval repo that last modified the corresponding Dockerfile, e.g., 65132f2

3. Now you just have to specify the docker image(s) in your execution workflow:
    - For [nextflow][nf-docker], the individual containers can be specified in the processes.
    - For [Snakemake][smk-docker], the individual containers can be specified per rule.
## Input
### Test data
For more information about input files, see ["sanctioned input files"](#more-details) above. For development and debugging you can use the small [test input dataset][test-data] we provide with this repository. You should use the `.bam`, `.fastq.gz` and/or`.gtf` files as input to your workflow. The `.bed` file serves as an example for a ground truth file.

### Parameters 
Both [snakemake template][snakemake-template] and [nextflow template][nextflow-template-dsl2] contain example `sample.csv` files. Here you'd fill in the paths to the samples you'd be running, and any other *sample specific* information required by the workflow you're implementing. Thus, you can/must adapt the fields of this `samples.csv` according to your workflow's requirements.   

Moreover, both workflow languages require additional information in `config` files. This is the place to specify *run- or participant-specific* parameters

>**Important notes:**   
>* Describe in your README extensively where parameters (sample info, participant specific parameters) have to be specified for a new run of the pipeline.
>* Describe in the README if your participant has different run modes, or parameter settings that might alter the participant's performance considerably. In such a case you should suggest that the different modes should be treated in APAeval as entirely *distinct participants*. Please raise such considerations in our slack space.
>* Parameterize your code as much as possible, so that the user will only have to change the sample sheet and config file, and *not the code*. E.g. output file paths should be built from information the user has filled into the sample sheet or config file.
>* For information on how files need to be named see [below](#output)!

## Output 
In principle you are free to store output files how it best suits you (or the participant). 
However, the "real" and final outputs for each run of the benchmarking will need to be *copied* to a directory in the format   
`PATH/TO/APAEVAL/CHALLENGECODE/PARTICIPANT/`

This directory *must* contain:
- Output files (check [formats](#formats) and [filenames](#filenames))
- Configuration files (with parameter settings), e.g. `config.yaml` and `samples.csv`.
- `logs/` directory with all log files created by the workflow exeuction.

### Formats
File formats for the 3 benchmarking events are described in the [output specification][spec-doc] which also contains the `OUTCODE` (01 - Identification, 02 - Quantification, 03 - Differential expression) needed for correct naming.
### Filenames
> As mentioned [above](#parameters) it is best to parameterize filenames, such that for each run the names and codes can be set by changing only the sample sheet and config file!

File names **must** adhere to the following schema: `CHALLENGECODE_PARTICIPANT_OUTCODE.ext`   
For the codes please refer to the following documents:   
- CHALLENGECODE: in [`challenge_codes.md`][challenge-code]
- PARTICIPANT: same as directory name in [`execution_workflows`][participant]
- OUTCODE: in [`execution_output_specification.md`][spec-doc]

**Example:**   
 `AA/MISO/AA_MISO_01.bed` would be the output of MISO (your participant) for the identification benchmarking event (OUTCODE 01, we know that from [`execution_output_specification.md`][spec-doc]), run on dataset "P19" using 4 cores (CHALLENGECODE AA, we know that from) [`summary_input_specification.md`][param-code])


## PR reviews
At least 2 independent reviews are required before your code can be merged into the main APAeval branch. Why not review some other PR while you wait for yours to be accepted? You can find some instructions in [Sam's PR review guide][pr-review-guide].

## Benchmarking Participants (bioinformatic methods)
List of bioinformatic methods benchmarked in APAeval. Please update columns as the execution workflows progress.

| Method | Citation | Type | Status in APAeval | OpenEBench link | Installation | Can take your :poodle: for a walk? |
|-|-|-|-|-|-|-|
| [APA-Scan](https://github.com/compbiolabucf/APA-Scan) | [Fahmi et al. 2020](https://www.biorxiv.org/content/10.1101/2020.02.16.951657v2) | Identification, Quantification | [Issue #26][issue-26] | https://dev-openebench.bsc.es/tool/apa-scan | github | Maybe |
| [APAlyzer](https://bioconductor.org/packages/release/bioc/html/APAlyzer.html) | [Wang & Tian 2020](https://pubmed.ncbi.nlm.nih.gov/32321166/) | Quantification, Differential Quantitation | [Issue #19][issue-19] | https://dev-openebench.bsc.es/tool/apalyzer | Bioconductor | No, but it can take *you* for a ride |
| [diffUTR](https://github.com/ETHZ-INS/diffUTR) | [Gerber et al. 2021](https://doi.org/10.1186/s12859-021-04114-7) | Differential Quantitation | [Issue #27][issue-27] | https://dev-openebench.bsc.es/tool/diffutr | Bioconductor | Only if you ask nicely |
| [QAPA](https://github.com/morrislab/qapa) | [Ha et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1414-4) | Quantification, Differential Quantitation | [Issue #28][issue-28] | https://openebench.bsc.es/tool/qapa | github | Probably not |
| [APAtrap](https://sourceforge.net/projects/apatrap/) | [Ye et al. 2018](https://academic.oup.com/bioinformatics/article/34/11/1841/4816794) | Identification (de novo), Quantification, Differential Quantitation | [Issue #29][issue-29] | NA | SourceForge | Nope |
| [Aptardi](https://github.com/luskry/aptardi) | [Lusk et al. 2021](https://www.nature.com/articles/s41467-021-21894-x) | Identification | [Issue #30][issue-30] | https://openebench.bsc.es/tool/aptardi | Bioconda | No, and neither could my :poodle: walk *it* |
| [CSI-UTR](https://github.com/UofLBioinformatics/CSI-UTR) | [Harrison et al. 2019](https://doi.org/10.3389/fgene.2019.00182) | Quantification, Differential Quantification | [Issue #31][issue-31] | NA | github | Maybe |
| [DaPars2](https://github.com/3UTR/DaPars2) | [Feng et al. 2018](https://academic.oup.com/nar/article/46/D1/D1027/4372484?) | Identification (de novo), Quantification, Differential Quantitation | [Issue #32][issue-32] | NA | github | No way |
| [GETUTR](http://big.hanyang.ac.kr/GETUTR/manual.htm) | [Kim et al. 2015](https://www.sciencedirect.com/science/article/abs/pii/S1046202315001577?via%3Dihub) | Identification, Quantification | [Issue #33][issue-33] | https://openebench.bsc.es/tool/getutr | lab website | Nope |
| [IsoSCM](https://github.com/shenkers/isoscm) | [Shenker et al. 2015](https://rnajournal.cshlp.org/content/21/1/14) | Identification (de novo), Differential Quantitation | [Issue #34][issue-34] | https://dev-openebench.bsc.es/tool/isoscm | github | Heck yes |
| [LABRAT](https://github.com/TaliaferroLab/LABRAT) | [Goering et al. 2020](https://www.biorxiv.org/content/10.1101/2020.10.05.326702v1) | Quantification, Differential Quantitation | [Issue #35][issue-35] | https://openebench.bsc.es/tool/labrat | Bioconda | Only if you ask nicely |
| [MISO](http://hollywood.mit.edu/burgelab/miso/) | [Katz et al. 2010](https://www.nature.com/articles/nmeth.1528) | Quantification, Differential Quantitation | [Issue #36][issue-36] | https://openebench.bsc.es/tool/miso | pypi or bioconda | No, and neither could my :poodle: walk *it* |
| [mountainClimber](https://github.com/gxiaolab/mountainClimber) | [Cass & Xiao 2019](https://www.sciencedirect.com/science/article/pii/S2405471219302686) | Identification (de novo), Quantification, Differential Quantitation | [Issue #37][issue-37] | https://openebench.bsc.es/tool/mountainclimber | github | Maybe |
| [Roar](https://github.com/vodkatad/roar/) | [Grassi et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1254-8) | Identification, Differential Quantitation | [Issue #38][issue-38] | https://openebench.bsc.es/tool/roar | Bioconductor | Probably not |
| [TAPAS](https://github.com/arefeen/TAPAS) | [Arefeen et al. 2018](https://academic.oup.com/bioinformatics/article/34/15/2521/4904269) | Identification (de novo), Quantification, Differential Quantitation | [Issue #39][issue-39] | https://openebench.bsc.es/tool/tapas | github | Nope |
| [DaPars](https://github.com/ZhengXia/dapars) | [Xia et al. 2014](https://www.nature.com/articles/ncomms6274) | Identification (de novo), Quantification, Differential Quantitation | [Issue #40][issue-40] | NA | github | Heck yes |
| [PAQR/KAPAC](https://github.com/zavolanlab/PAQR_KAPAC) | [Gruber et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1415-3) | Quantification, Differential Quantitation | [Issue #41][issue-41] | https://openebench.bsc.es/tool/kapac | miniconda | No way |
| [InPAS](http://www.bioconductor.org/packages/release/bioc/vignettes/InPAS/inst/doc/InPAS.html) | [Sheppard et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789547/) | Identification (de novo), Quantification, Differential Quantitation |  | https://openebench.bsc.es/tool/inpas | Bioconductor | Yep |

*Note: Can take your :poodle: for a walk was assigned randomly for initial table*

[//]: # (References)

[apaeval-ewfs]: ../images/EWFs.png 
[conda]: <https://docs.conda.io/en/latest/>  
[docker]: <https://www.docker.com/>
[docker-contact-alex]: <https://app.slack.com/client/T01PW9SAN7K/D01PP4WK7TL/user_profile/U01PEJ5TW4V>
[docker-contact-yukkei]: <https://app.slack.com/client/T01PW9SAN7K/D01PP4WK7TL/user_profile/U01SFJM5FM5>
[dockerhub]: <https://hub.docker.com/>
[dockerfile]: ../docs/templates/snakemake/workflow/envs/[METHOD].Dockerfile
[docker-tutorial]: <https://stackify.com/docker-build-a-beginners-guide-to-building-docker-images/>
[biocontainers]: <https://biocontainers.pro/registry>
[nf]: <https://www.nextflow.io/>
[nf-docker]: <https://www.nextflow.io/docs/latest/docker.html#multiple-containers.>
[nextflow-template-dsl1]: <https://github.com/iRNA-COSI/APAeval/tree/main/docs/templates/nextflow_dsl1>
[nextflow-template-dsl2]: <https://github.com/iRNA-COSI/APAeval/tree/main/docs/templates/nextflow_dsl2>
[spec-doc]: execution_output_specification.md 
[challenge-code]: ../summary_workflows/challenge_codes.md
[participant]: ../execution_workflows/
[pr-review-guide]: ./execution_workflows/PR_review_guide.md
[singularity]: <https://sylabs.io/singularity/>
[snakemake-template]: <https://github.com/iRNA-COSI/APAeval/docs/templates/snakemake>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[smk-docker]: <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#running-jobs-in-containers>
[test-data]: ../tests/test_data
[issue-26]: <https://github.com/iRNA-COSI/APAeval/issues/26>
[issue-19]: <https://github.com/iRNA-COSI/APAeval/issues/19>
[issue-27]: <https://github.com/iRNA-COSI/APAeval/issues/27>
[issue-28]: <https://github.com/iRNA-COSI/APAeval/issues/28>
[issue-29]: <https://github.com/iRNA-COSI/APAeval/issues/29>
[issue-30]: <https://github.com/iRNA-COSI/APAeval/issues/30>
[issue-31]: <https://github.com/iRNA-COSI/APAeval/issues/31>
[issue-32]: <https://github.com/iRNA-COSI/APAeval/issues/32>
[issue-33]: <https://github.com/iRNA-COSI/APAeval/issues/33>
[issue-34]: <https://github.com/iRNA-COSI/APAeval/issues/34>
[issue-35]: <https://github.com/iRNA-COSI/APAeval/issues/35>
[issue-36]: <https://github.com/iRNA-COSI/APAeval/issues/36>
[issue-37]: <https://github.com/iRNA-COSI/APAeval/issues/37>
[issue-38]: <https://github.com/iRNA-COSI/APAeval/issues/38>
[issue-39]: <https://github.com/iRNA-COSI/APAeval/issues/39>
[issue-40]: <https://github.com/iRNA-COSI/APAeval/issues/40>
[issue-41]: <https://github.com/iRNA-COSI/APAeval/issues/41>
[utils]: ../utils




