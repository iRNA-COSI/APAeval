# Methods directory
This is where the execution workflows for APAeval live. 

Besides the method execution, an execution workflow has to contain the pre-processing steps to convert the given input files [LINK TO INPUTS] to formats that the method consumes, if any, and the post-processing steps required to create all output files [LINK TO OUTPUTS] needed to run all summary workflows of challenges that are applicable to the method (e.g. not every method performs PAS identification, etc.).


To implement an execution workflow for a method, copy either the [snakemake template](https://github.com/iRNA-COSI/APAeval/docs/templates/snakemake) or the [nextflow template](https://github.com/iRNA-COSI/APAeval/docs/templates/nextflow) into the method's directory and adapt the workflow directory names as described in the template's `README`.   


Example:    
```bash
methods/
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

## Output 

You are free to store output files how it best suits you (or the method). 
However, final outputs need to be copied to a directory in this format:
"PATH/TO/s3-BUCKET/PARAM_CODE/METHOD/"

This directory must contain:
- Output files in the format described in the [output specification][spec-doc]. File names must adhere to the following scheme: PARAM_CODE_METHOD_OUTCODE.ext.
- Configuration files (includes parameter settings), e.g. `config.yaml` and `samples.csv`.
- `logs/` directory with all log files created by the workflow exeuction.

### Where to find the codes

- PARAM_CODE: in [`summary_input_specification.md`][param-code]
- METHOD: same as directory name in [`execution_workflows`][method]
- OUTCODE: in [`execution_output_specification.md`][outcode]

[//]: # (References)
  
[spec-doc]: /execution_workflows/execution_output_specification.md 
[param-code]: /summary_workflows/summary_input_specification.md
[method]: /execution_workflows/
[outcode]: /execution_workflows/execution_output_specification.md