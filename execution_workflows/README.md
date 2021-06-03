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

| Method | Citation | Repository | Type | Status in APAeval | Present in OpenEBench | Installation | Can take your dog for a walk? |
|-|-|-|-|-|-|-|-|
| APA-Scan | [Fahmi et al. 2020](https://www.biorxiv.org/content/10.1101/2020.02.16.951657v2) | [APA-Scan](https://github.com/compbiolabucf/APA-Scan) | Identification, Quantification | In progress | Yes | github | Maybe |
| APAlyzer | [Wang & Tian 2020](https://pubmed.ncbi.nlm.nih.gov/32321166/) | [APAlyzer](https://bioconductor.org/packages/release/bioc/html/APAlyzer.html) | Quantification, Differential Quantitation | On hold | Yes | Bioconductor | No, but it can take *you* for a ride |
| diffUTR | [Gerber et al. 2021](https://doi.org/10.1186/s12859-021-04114-7) | [diffUTR](https://github.com/ETHZ-INS/diffUTR) | Differential Quantitation | To do | Yes | Bioconductor | Only if you ask nicely |
| QAPA | [Ha et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1414-4) | [QAPA](https://github.com/morrislab/qapa) | Quantification, Differential Quantitation | In progress | Yes | github | Probably not |
| APAtrap | [Ye et al. 2018](https://academic.oup.com/bioinformatics/article/34/11/1841/4816794) | [APAtrap](https://sourceforge.net/projects/apatrap/) | Identification (de novo), Quantification, Differential Quantitation | In progress | No | SourceForge | Nope |
| Aptardi  | [Lusk et al. 2021](https://www.nature.com/articles/s41467-021-21894-x) | [Aptardi](https://github.com/luskry/aptardi) | Identification | To do | Yes | Bioconda | No, and neither could my dog walk *it* |
| CSI-UTR | [Harrison et al. 2019](https://doi.org/10.3389/fgene.2019.00182) | [CSI-UTR](https://github.com/UofLBioinformatics/CSI-UTR) | Quantification, Differential Quantification | In progress | No | github | Maybe |
| DaPars2 | [Feng et al. 2018](https://academic.oup.com/nar/article/46/D1/D1027/4372484?) | [DaPars2](https://github.com/3UTR/DaPars2) | Identification (de novo), Quantification, Differential Quantitation | In progress | No | github | No way |
| GETUTR | [Kim et al. 2015](https://www.sciencedirect.com/science/article/abs/pii/S1046202315001577?via%3Dihub) | [GETUTR](http://big.hanyang.ac.kr/GETUTR/manual.htm) | Identification, Quantification | To do | Yes | lab website | Nope |
| IsoSCM | [Shenker et al. 2015](https://rnajournal.cshlp.org/content/21/1/14) | [IsoSCM](https://github.com/shenkers/isoscm) | Identification (de novo), Differential Quantitation | To do | Yes | github | Heck yes |
| LABRAT | [Goering et al. 2020](https://www.biorxiv.org/content/10.1101/2020.10.05.326702v1) | [LABRAT](https://github.com/TaliaferroLab/LABRAT) | Quantification, Differential Quantitation | In progress | Yes | Bioconda | Only if you ask nicely |
| MISO | [Katz et al. 2010](https://www.nature.com/articles/nmeth.1528) | [MISO](http://hollywood.mit.edu/burgelab/miso/) | Quantification, Differential Quantitation | On hold | Yes | pypi or bioconda | No, and neither could my dog walk *it* |
| mountainClimber | [Cass & Xiao 2019](https://www.sciencedirect.com/science/article/pii/S2405471219302686) | [mountainClimber](https://github.com/gxiaolab/mountainClimber) | Identification (de novo), Quantification, Differential Quantitation | To do | Yes | github | Maybe |
| Roar | [Grassi et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1254-8) | [Roar](https://github.com/vodkatad/roar/) | Identification, Differential Quantitation | In progress | Yes | Bioconductor | Probably not |
| TAPAS | [Arefeen et al. 2018](https://academic.oup.com/bioinformatics/article/34/15/2521/4904269) | [TAPAS](https://github.com/arefeen/TAPAS) | Identification (de novo), Quantification, Differential Quantitation | In progress | Yes | github | Nope |
| DaPars | [Xia et al. 2014](https://www.nature.com/articles/ncomms6274) | [DaPars](https://github.com/ZhengXia/dapars) | Identification (de novo), Quantification, Differential Quantitation | On hold | No | github | Heck yes |
| PAQR/KAPAC | [Gruber et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1415-3) | [PAQR/KAPAC](https://github.com/zavolanlab/PAQR_KAPAC) | Quantification, Differential Quantitation | On hold | Yes | miniconda | No way |
| InPAS | [Sheppard et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789547/) | [InPAS](http://www.bioconductor.org/packages/release/bioc/vignettes/InPAS/inst/doc/InPAS.html) | Identification (de novo), Quantification, Differential Quantitation | To do | Yes | Bioconductor | Maybe |

*Note: Can take your dog for a walk was assigned randomly for initial table*

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
