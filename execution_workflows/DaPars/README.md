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

### Steps to run this:
 - DaPars requires exactly two distinct conditions to be provided for the differential step
 - samplesheet_example_files.csv contains paths to all the input test files provided for each workflow run
 - To run DaPars with test data provided for APAeval, check the path to DaPars with `pwd` and replace the `path_to` in samplesheet_example_files.csv with the 
   path from the `pwd` command. Otherwise, replace the bam/bai paths in the samplesheet with the absolute paths
   to the files you would like to use
 - Each row in the samplesheet contains files for one sample to be processed
 - Any number of rows (i.e. any number of replicates per condition) can be provided in the samplesheet in any order
 - The first step of DaPars requires a gene symbol file. Download gene symbol file from 
 [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). Follow the guide from this [thread](https://www.biostars.org/p/92939/)
 that shows how to download a gene symbol file form UCSC table browser
 - Make sure to specify parameters used for DaPars run in conf/modules.config file
 - Then, you are good to run the pilot benchmark nextflow pipeline with `DaPars`

```
nextflow main.nf --input samplesheet_example_files.csv
```

## Params
Parameters used to run the two steps of DaPars are specified in conf/modules.config file.

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

## Notes
It is not possible to obtain identification data since only proximal APA sites shared by
both groups are shown in the final output file. Since only a subset of APA sites both conditions
are shown, it doesn't qualify for identification challenge and it only qualifies for differential
challenge


## Author contact
If you have any question or comment about DaPars, please contact Zheng Xia (zxia@bcm.edu) or Wei Li (wl1@bcm.edu).