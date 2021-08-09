# APAtrap
The APAtrap R package does the following: 
1. Refines annotated 3'UTR and identifies novel 3' UTRs and 3' UTR extensions
2. Aims to identify all potential APA (alternative polyadenylation) sites
3. Detects genes with differential APA site usage between conditions by leveraging 
   the resolution of RNA-seq data

The [paper](https://academic.oup.com/bioinformatics/article/34/11/1841/4816794) is titled APAtrap: identification and quantification of 
alternative polyadenylation sites from RNA-seq data
The [application download](https://sourceforge.net/projects/apatrap/files/) is free to download
and [user manual](https://sourceforge.net/p/apatrap/wiki/User%20Manual/) was used as a referrence
to create the nextflow pipeline flow of this module

### Steps to run this:
 - Replace "path_to" in samplesheet_example_files.csv with the path to the reference input files
 - Check the path to the input files with `pwd` and replace the `path_to` in samplesheet_example_files.csv with the 
   path from the `pwd` command
 - Go to `conf/modules.config` to configure the parameters required to run APAtrap workflow. Descriptions of the parameters
   are located in the file
 - Then, you are good to run the pilot benchmark nextflow pipeline with `APAtrap`
```
nextflow main.nf --input samplesheet_example_files.csv
```

### Docker containers
This workflow uses docker container. To run, make sure that docker is installed and running
 
## Input & pre-processing
Required files to be specified in the input `samplesheet_example_files.csv`:

- sample: name of the run, will be used as the folder name of final output files 
- bam1: BAM input file for sample 1
- bai1: BAI index file for sample1's bam input
- bam2: BAM input file for sample 2
- bai2: BAI index file for sample2's bam input
- bed: gene model file in bed format. Can be obtained from [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1133780495_lZCAEdlwBd7HbE03thrN4Tsi6lSF&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=wgEncodeGencodeVM18&hgta_table=0&hgta_regionType=genome&position=chr12%3A56%2C694%2C976-56%2C714%2C605&hgta_outputType=bed&hgta_outFileName=)

## Params
In this workflow, default parameters are used to run all three steps of APAtrap


## Output & post-processing
Each APAtrap run results in files for all three challenges:
- apatrap_identification_output.bed
- apatrap_quantification_output.bed
- apatrap_differential_output.tsv 

The output files are located under APAtrap/results/<SAMPLE> where <SAMPLE>
is the sample specified on the first column of the input samplesheet_example_files.csv
file. This is done to differentiate the results from the different APAtrap runs.
So sample should be specified as a unique name as to not overwrite
results from different runs

## Notes
The "chr" in the chromosome column in the input bam file needs to be there. Otherwise, identifyDistal3UTR will
product an empty output file. The check for leading 'chr' in the chromosome column is done under 
PREPROCESSING step. If no leading chr is detected, it'll be added

