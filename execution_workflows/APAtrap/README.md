# APAtrap
The APAtrap R package does the following three steps: 
1. Refines annotated 3'UTR and identifies novel 3' UTRs and 3' UTR extensions
2. Aims to identify all potential APA (alternative polyadenylation) sites
3. Detects genes with differential APA site usage between conditions by leveraging 
   the resolution of RNA-seq data

The [paper](https://academic.oup.com/bioinformatics/article/34/11/1841/4816794) is titled APAtrap: identification and quantification of 
alternative polyadenylation sites from RNA-seq data <br>
The [application download](https://sourceforge.net/projects/apatrap/files/) is free to download
and [user manual](https://sourceforge.net/p/apatrap/wiki/User%20Manual/) was used as a referrence
to create the nextflow pipeline flow of this module

## Running APAtrap workflow

## Input & pre-processing
Required files are to be specified in the input `samplesheet_example_files.csv`. Each row in the sample sheet has four
columns:

- condition: name of the condition (e.g Control)
- sample: unique name of the sample for logs (e.g control_replicate1)
- bam: BAM input file for the sample 
- bai: BAI index file for sample's bam input

It is important to name samples of the same condition with the exact condition name under the condition
column in the sample sheet since samples are grouped per condition to be processed in the differential step. In addition,
sample names should be unique.

To run APAtrap with test data provided for APAeval, check the path to APAtrap with `pwd` and replace 
the `path_to` in samplesheet_example_files.csv with the path from the `pwd` command. 

When using your own data and input file instead of the provided test data and `samplesheet_example_files.csv`, make sure to include in the 
input file you are using the absolute path to the file, with the four column names following the column
names in `samplesheet_example_files.csv`.

### Running with Docker or Singularity
## Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running 
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
Additionally, make sure that line 49 in Apatrap/nextflow.config file `docker.enabled=true` is uncommented while line
51 `singularity.enabled=true` is commented out.

## Singularity
To run with singularity, comment out line 49 in Apatrap/nextflow.config file `docker.enabled=true` and make sure that line
51 `singularity.enabled=true` is uncommented.

### Parameters
Parameters used to run APAtrap are specified in conf/modules.config file. 
Parameters relevant to the workflow itself are:
- `run_identification` - set to true to obtain identification challenge output. Specifying any other value will throw an error.
- `run_quantification` - set to true to obtain quantification challenge output. Specifying any other value will throw an error.
- `run_differential` - set to true to obtain differential challenge output. Specifying any other value will throw an error.
- `output_dir` - name of the folder that the final output files are going to be in, located under Apatrap/results/apatrap/
- `identification_out_suffix` - suffix of the output file(s) for the current run ending with .bed when running identification,
                                the prefix will be the different sample names obtained from the sample column in the sample sheet
- `quantification_out_suffix` - suffix of the output file(s) for the current run ending with .bed when running quantification,
                                the prefix will be the different sample names obtained from the sample column in the sample sheet
- `differential_out` - name of the output file for the current run ending with .tsv when running differential
- `genome_file` - absolute path from the APAtrap folder to the input GTF annotation file can be obtained by replacing `path_to`
   with the path to APAtrap, and if using your own genome file, make sure to use the absolute path to your genome file

### Running the differential workflow
- Set the 'run_differential' parameter in conf/modules.config to true
- Change 'differential_out' parameter in conf/modules.config to the desired file name that ends with '.tsv'
- Ensure the sample sheet contains exactly two distinct conditions in the condition column

### Running the identification and quantification workflow
- Set the 'run_identification' and/or 'run_quantification' parameters in conf/modules.config to true to 
  run identification and/or quantification workflows
- Change 'identification_out_suffix' and/or 'quantification_out_suffix' parameters in conf/modules.config to 
  the desired file suffix that ends with '.bed'
- In the sample sheet, at least one row is required to be provided 


### Running the workflow
Once parameters have been set in conf/modules.config file, run the pilot benchmark nextflow pipeline with 
`nextflow main.nf --input samplesheet_example_files.csv`. 

Note that it's recommended to delete the existing `results` folder to remove results from the previous run, if any. 
This is to avoid polluting the next run with the results from the previous run when we are using the same sample name(s),
which would likely result in an error.

## Output & post-processing
When using the default `output_dir` parameter value in conf/modules.config, APAtrap run results in files of the outputs 
of the challenges located under APAtrap/results/apatrap/challenges_outputs folder.
For identification outputs, the files have sample names as prefixes to differentiate the different runs.
The differential output file will stay as the name specified in modules.config file.


## Notes
- When running differential, as clarified with APAtrap author, Dr. Congting Ye, predictAPA step requires that the input files are placed
  by group and in the same order of values specified by parameter -n. So for example, if we have two replicates for each condition, we'll then do
  `predictAPA -i condition1_replicate1.bedgraph condition1_replicate2.bedgraph condition2_replicate1.bedgraph condition2_replicate2.bedgraph -g 2 -n 2 2 -u  hg19.utr.bed -o output.txt`. 
  This means that we have to loop through the sample files by condition, where each condition has a folder in the results folder.
  As such, the results folder needs to be clean of files and folders from the previous run as to not pollute the predictAPA step of the next run. This is
  the reason we require the results folder to be deleted before running the workflow.

## Author contact
If you have any question or comment about APAtrap, please contact with Dr. Congting Ye (yec@xmu.edu.cn).

