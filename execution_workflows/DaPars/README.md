# DaPars
DaPars is a python-based model that directly infers the dynamic alternative polyadenylation (APA)
usage by comparing standard RNA-seq. Given the annotated gene model, DaPars can infer the de novo proximal APA sites 
as well as the long and short 3’UTR expression levels. Finally, it identifies the dynamic APA usages between two 
conditions. 

The [paper](https://www.nature.com/articles/ncomms6274) is titled Dynamic analyses of alternative polyadenylation from 
RNA-seq reveal a 3′-UTR landscape across seven tumour types <br>
The [application download](https://github.com/ZhengXia/dapars) is free to download
and [documentation](http://xiazlab.org/dapars_tutorial/html/DaPars.html) was used as a reference
to create the nextflow pipeline flow of this module

## Running DaPars workflow

### Modes:
This workflow has two modes. The first mode is for identification challenge, and the second mode is
for differential challenge.

### Parameters
Parameters used to run the two steps of DaPars are specified in conf/modules.config file. This file includes the 
parameter called 'mode' to determine whether to obtain identification or differential challenge output.
When 'mode' parameter is set to 'differential', the workflow to obtain differential challenge
output will be run. When 'mode' parameter is set to 'identification', the workflow to obtain identification
challenge will be run. When mode is set to anything other than 'differential' or 'identification'. the workflow
will throw an error. This file also has a parameter called 'genome_file' which holds the relative path from DaPars 
to the gtf genome file to be used. The relative path is the path that leads to the gtf genome file starting from
the DaPars folder.

### Steps to run this:
 - When 'mode' parameter is set to 'differential', DaPars requires exactly two distinct conditions to be provided. An 
   example input file is samplesheet_example_files.csv, which contains paths to all the input test files needed for 
   each differential workflow run for two distinct conditions: Control and Srsf
 - When 'mode' parameter is set to 'identification', DaPars requires the same condition to be provided twice. An 
   example input file is samplesheet_example_files_identification.csv, which contains paths to all the input test files needed for 
   each identification workflow run for the condition: Control. The only difference between the two lines
   is the 'condition' column, which has to be distinct. The workflow will then treat the two rows as two different
   conditions, a requirement for DaPars to run successfully.
 - To run DaPars with test data provided for APAeval, check the path to APAeval/execution_workflows/DaPars with `pwd` and replace the `path_to` 
   in samplesheet_example_files.csv or samplesheet_example_files_identification.csv with the 
   path from the `pwd` command. 
 - Note that for the 'identification' mode, exactly two rows have to be provided in the samplesheet. For
   'differential' mode, at least two rows have to be provided containing exactly two distinct conditions.
 - Make sure to specify parameters used for DaPars run in conf/modules.config file before running the workflow.
 - Then, you are good to run the pilot benchmark nextflow pipeline with `DaPars`

```
nextflow main.nf --input samplesheet_example_files.csv
```
when 'mode' is set to 'differential' <br>
OR
```
nextflow main.nf --input samplesheet_example_files_identification.csv
```
when 'mode' is set to 'identification'

### Docker containers
This workflow uses docker container. To run, make sure that docker is installed and running

## Input & pre-processing
Required files are to be specified in the input `samplesheet_example_files.csv`. Each row in the samplesheet has four
columns:

- sample: name of the sample for logs (e.g control_replicate1)
- condition: name of the condition (e.g Control) 
- bam: BAM input file for the sample 
- bai: BAI index file for sample's bam input

Each row in the samplesheet has the files needed for each sample file.

It is important to name samples of the same condition with the exact condition name under the condition
column in the samplesheet since samples are grouped per condition to be processed in the differential step.

## Output & post-processing
Each DaPars run results in differential challenge file located under DaPars/results/dapars/dapars_differential_output.tsv.
Identification challenege file is located under DaPars/results/dapars/dapars_identification_output.tsv.

## Notes
- It is not possible to obtain quantification challenge output data since the TPM value in the output file
  is provided per transcript instead of per site. Hence, this tool is not compatible for quantification
  challenge. 
- Make sure that the input bam files have leading 'chr' in the chromosome column. Otherwise, once 
  converted to input bedgraph file for DaPars, the workflow will exit with an error
![](../../../../Screen Shot 2021-09-21 at 1.08.09 PM.png)
   This is because DaPars checks for the leading 'chr' in the process and would error out otherwise.
   Hence, we've added a step to check for this leading 'chr' early in the workflow to prevent having to
   go through the entire workflow before erroring out for efficiency.

## Author contact
If you have any question or comment about DaPars, please post on DaPars Google Groups (https://groups.google.com/u/1/g/DaPars) or the author, Dr. Zheng Xia (xiaz@ohsu.edu).
