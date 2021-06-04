# Snakemake pipeline execution

Methods are saved in `workflows`, and corresponding conda environments in `envs`.
Configuration file paths and samples table paths are saved in `configs`. 

## Installation with conda

Create a conda environment with

``` bash
conda env create -f envs/snakemake.yaml
conda activate snakemake
```

# Benchmark Q1 on Miso

Adjust file paths for dataset run in `configs/config.miso.yaml`. 
* Miso also requires an own configuration file and is specified in `configs/miso_settings.txt`.
* The samples table `configs/samples.csv` must include the following columns: **sample** and **bam**, for the sample name and absolute path to BAM file respectively. 

Run the MISO pipeline with conda environment:

``` bash
bash run_miso.sh
```

Run the MISO pipeline with singularity containers:

``` bash
bash run_miso_singularity.sh
```

> Note: Running this requires that singularity is installed.

The benchmark output file is `AA_MISO_04.json` and follows the benchmark specification. 
It contains the sum of run times for all samples and the index generation and the maximum Proportional Set Size (PSS) from all rules.

## Current challenges

* CPU time seems to be not correctly recorded with snakemake's `benchmark` directive when using multiple threads. 
* The pipeline is tailored to the pilot dataset, which is unstranded in single-end. It might be necessary to create new rules to incorporate this.