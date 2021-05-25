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

