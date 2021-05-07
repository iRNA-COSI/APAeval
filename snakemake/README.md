# Snakemake pipeline execution

Methods are saved in `workflows`, and corresponding conda environments in `envs`.
Configuration file paths and samples table paths are saved in `configs`. 

## Installation with conda

Create a conda environment with

``` bash
conda env create -f envs/snakemake.yaml
conda activate snakemake
```

# Miso

Adjust file paths for dataset run in `configs/config.miso.yaml`. 
* Miso also requires an own configuration file and is specified in `configs/miso_settings.txt`.
* The samples table `configs/samples.csv` must include the following columns: **sample** and **bam**, for the sample name and absolute path to BAM file respectively. 

Run the MISO pipeline with

``` bash
bash run_miso.sh
```

The benchmark output file is `benchmark.Q1_miso.json`