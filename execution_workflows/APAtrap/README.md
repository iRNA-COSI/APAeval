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

### Steps to run this:
 - We can either run the workflow to obtain identification and quantification outputs, or differential output. There
   is currently no way to run the workflow once to get outputs for all three challenges
 - samplesheet_example_files.csv contains paths to all the input files provided for each workflow run
 - Replace "path_to" in samplesheet_example_files.csv with the path to the reference input files
 - Check the path to the input files with `pwd` and replace the `path_to` in samplesheet_example_files.csv with the 
   path from the `pwd` command
 - Each row in the samplesheet contains files for one sample to be processed
 - Any number of rows can be provided in the samplesheet
 - If we are running the differential step of APAtrap, make sure to have at least two rows in the samplesheet and all
   the rows are processed at the same time
 - If we are not running the differential step of APAtrap, each row in the samplesheet will be processed individually 
   for identification and quantification challenges
 - Go to `conf/modules.config` to configure the parameters required to run APAtrap workflow. Descriptions of the parameters
   are located in the file. Head to `Params` for more info
 - Then, you are good to run the pilot benchmark nextflow pipeline with `APAtrap`

```
nextflow main.nf --input samplesheet_example_files.csv
```

## Params
Parameters used to run all three steps of APAtrap are specified in conf/modules.config file. In the file, under workflow,
there is a parameter called run_differential. When set to true, only the output file for differential challenge is produced.
When set to false, only the output files for identification and quantification challenges are produced. The separation
between identification and quantification versus differential runs is necessary since the inputs for these two modes are
different. 

### Docker containers
This workflow uses docker container. To run, make sure that docker is installed and running
 
## Input & pre-processing
Required files are to be specified in the input `samplesheet_example_files.csv`. Each row in the samplesheet has four
columns:

- sample: name of the sample
- bam: BAM input file for the sample 
- bai: BAI index file for sample's bam input
- bed: gene model file in bed format. Can be obtained from [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1133780495_lZCAEdlwBd7HbE03thrN4Tsi6lSF&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=wgEncodeGencodeVM18&hgta_table=0&hgta_regionType=genome&position=chr12%3A56%2C694%2C976-56%2C714%2C605&hgta_outputType=bed&hgta_outFileName=)

Each row in the samplesheet has the files needed for each sample file.
When run_differential is set to false, at least one row is required to be provided.
When run_differential is set to true, at least two rows (i.e. two conditions) are required to be provided. 

## Output & post-processing
Each APAtrap run results in files of the outputs of the challenges located under APAtrap/results/apatrap/challenges_outputs folder.
For identification and quantification outputs, the files have sample names as prefixes to differentiate the different runs.
The differential output file will stay as the name specified in modules.config file.


## Notes
The "chr" in the chromosome column in the input bam file needs to be there. Otherwise, identifyDistal3UTR will
product an empty output file. The check for leading 'chr' in the chromosome column is done under 
PREPROCESSING step. If no leading chr is detected, it'll be added.

