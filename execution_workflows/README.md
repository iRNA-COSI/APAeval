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
