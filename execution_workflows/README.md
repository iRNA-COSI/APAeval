# Execution workflows
This is where the execution workflows for APAeval live. 

## Overview
_Execution workflows_ contain all steps that need to be run _per method_:

1. **Pre-processing:** Convert the input files the APAeval team has prepared
  into the input files a given method consumes, if applicable. This does not include e.g. adapter trimming or mapping of reads, as those steps are already performed in our general pre-processing pipeline. Pre-processing here means you have to convert the provided `.bam`, `.fastq.gz`, `.gtf` or `.gff` files to a format that your method can use.
2. **Method execution:** Execute the method in any way necessary to compute the
  output files for _all_ challenges (may require more than one run of the tool
  if, e.g., run in different execution modes).
3. **Post-processing:** Convert the output files of the method into the [formats][spec-doc]
  consumed by the _summary workflows_ as specified by the APAeval team, if
  applicable.

_Execution workflows_ should be implemented in either [Nexflow][nf] or
[Snakemake][snakemake], and individual steps should be isolated through the
use of either [Conda][conda] virtual environments (deprecated; to run on AWS we need containerized workflows) or
[Docker][docker]/[Singularity][singularity] containers.

## Templates
To implement an execution workflow for a method, copy either the [snakemake template][snakemake-template] or the [nextflow template][nextflow-template] into the method's directory and adapt the workflow directory names as described in the template's `README`. Don't forget to adapt the `README` itself as well.  


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

## Input
### Test data
For development and debugging you can use the small [test input dataset][test-data] we provide with this repository. You should use the `.bam`, `.fastq.gz`, `.gtf` and/or `.gff` files as input to your workflow. The `.bed` file serves as an example for a ground truth file.

### Parameters
Both [snakemake template][snakemake-template] and [nextflow template][nextflow-template] contain example `sample.csv` files. Here you'd fill in the paths to the samples you'd be running, and any other *sample specific* information required by the method you're implementing. Thus, you can/must adapt the fields of this `samples.csv` according to your workflow's requirements.   

Moreover, both workflow languages require additional information in `config` files. This is the place to specify *run- or method-specific* parameters

**Important notes:**   
* Describe in your README extensively where parameters (sample info, method specific parameters) have to be specified for a new run of the pipeline.
* Parameterize your code as much as possible, so that the user will only have to change the sample sheet and config file, and *not the code*. E.g. output file paths should be built from information the user has filled into the sample sheet or config file.
* For information on how files need to be named see [below](#output)!

## Output 
In principle you are free to store output files how it best suits you (or the method). 
However, the "real" and final outputs for each run of the benchmarking will need to be *copied* to a directory in the format   
`PATH/TO/s3-BUCKET/PARAMCODE/METHOD/`

This directory *must* contain:
- Output files (check [formats](#formats) and [filenames](#filenames))
- Configuration files (with parameter settings), e.g. `config.yaml` and `samples.csv`.
- `logs/` directory with all log files created by the workflow exeuction.

### Formats
File formats for the 3 challenges are described in the [output specification][spec-doc] which also contains the `OUTCODE` needed for correct naming.
### Filenames
> As mentioned [above](#parameters) it is best to parameterize filenames, such that for each run the names and codes can be set by changing only the sample sheet and config file!

File names **must** adhere to the following schema: `PARAMCODE_METHOD_OUTCODE.ext`   
For the codes please refer to the following documents:   
- PARAMCODE: in [`summary_input_specification.md`][param-code]
- METHOD: same as directory name in [`execution_workflows`][method]
- OUTCODE: in [`execution_output_specification.md`][outcode]

**Example:**   
 `AA/MISO/AA_MISO_01.bed` would be the output of MISO (your method) for the identification challenge (OUTCODE 01, we know that from [`execution_output_specification.md`][outcode]), run on dataset "P19" using 4 cores (PARAMCODE AA, we know that from) [`summary_input_specification.md`][param-code])


## Tools
List of tools used in APAeval. Please update columns as the execution workflows progress.

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
 
[conda]: <https://docs.conda.io/en/latest/>  
[docker]: <https://www.docker.com/>
[nf]: <https://www.nextflow.io/>
[nextflow-template]: <https://github.com/iRNA-COSI/APAeval/docs/templates/nextflow>
[spec-doc]: /execution_workflows/execution_output_specification.md 
[param-code]: /summary_workflows/parameter_codes.md
[method]: /execution_workflows/
[outcode]: /execution_workflows/execution_output_specification.md
[singularity]: <https://sylabs.io/singularity/>
[snakemake-template]: <https://github.com/iRNA-COSI/APAeval/docs/templates/snakemake>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[test-data]: /tests/test_data
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





