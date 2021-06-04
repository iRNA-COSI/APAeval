# Snakemake pipeline execution

Methods are saved in `workflow`, and corresponding conda environments in `envs`. `Snakefile` is on the working directory.
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

Run the MISO pipeline with

``` bash
bash run_miso.sh
```

The benchmark output file is `benchmark.Q1_miso.json` and follows the benchmark specification. 
It contains the sum of run times for all samples and the index generation and the maximum Proportional Set Size (PSS) from all rules.

# Running on AWS

## Prerequesites

AWS is set up and two IAM users are available, one with administration permissions and one with working permissions. Please refer to the respective documentation on how to do this.

AWS command line tools and tibanna are installed. If not, append it in your `snakemake` environment:

``` bash
conda activate snakemake
conda install awscli=1.19.78 tibanna
```

Configure your AWS profile for the admin (necessary to deploy tibanna):

``` bash
aws configure --profile admin
```

It will guide you through configuration, with default region and output format. 
> Note: The access keys you can download as csv file from the AWS console.

Do the same for the work profile:

``` bash
aws configure --profile work
```

Check if you can access the s3 bucket with

``` bash
aws s3 ls s3://<bucket>
```

## Set-up and deploy tibanna

This section is based on [executing a snakemake workflow via tibanna on AWS](https://snakemake.readthedocs.io/en/stable/executing/cloud.html?#executing-a-snakemake-workflow-via-tibanna-on-amazon-web-services).

Activate the AWS profile with admin permissions

``` bash
export AWS_PROFILE=admin
```

Deploy tibanna with usergroup name and the s3 bucket and add user **work** to the usergroup.

``` bash
tibanna deploy_unicorn -g <name> -b <bucket>
tibanna add_user -u work -g <name>
```

Next, the function that tibanna should execute must be defined. You find the following line in `~/.bashrc`

``` bash
export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=tibanna_unicorn_<name>
```

## Run

Copy `test_data` into this directory. This is necessary, because all source files need to be in the same working directory as the `Snakefile`.
Replace bucketname and subdir in `run_miso_tibanna.sh` and run the workflow with

``` bash
bash run_miso_tibanna.sh
```

# Current challenges

* CPU time seems to be not correctly recorded with snakemake's `benchmark` directive when using multiple threads. 
* The pipeline is tailored to the pilot dataset, which is unstranded in single-end. It might be necessary to create new rules to incorporate this.