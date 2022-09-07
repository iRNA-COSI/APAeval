# TAPAS
TAPAS is an R-based software that detects polyadenylation sites within a gene from RNA-Seq data and identifies differentially polyadenylated sites between two samples.

The [paper](https://academic.oup.com/bioinformatics/article/34/15/2521/4904269) is titled TAPAS: tool for alternative polyadenylation site analysis <br>
The [application](https://github.com/arefeen/TAPAS) is free to download,
and the [README documentation](https://github.com/arefeen/TAPAS#tapas-tool-for-alternative-polyadenylation-site-analysis) was used as a reference
to create the nextflow pipeline of this module.

This workflow qualifies for the APAeval **identification and quantification challenges**.

## Running TAPAS workflow

### Input & pre-processing
An example sample sheet is available at `samplesheet.csv`. `samplesheet.csv` can contain multiple entries, where each row in the samplesheet has two
columns:

- sample: name of the sample for logs (e.g control_replicate1)
- bam: BAM input file for the sample
- read_length: read length of the sample

### Docker and Singularity containers
This workflow uses docker containers. To run, make sure that docker is installed and running
(e.g. by running the command `docker --help` and seeing a help message printed).
- If running with Docker, please include the `-profile docker` in the command, which enables Docker.
- If running with Singularity, please include the `-profile singularity` in the command, which enables Singularity.

### Parameters
Parameters used to run the TAPAS are specified in the nextflow.config file.
Parameters relevant to the workflow itself are:
- `input` - samplesheet.csv

### Running the TAPAS execution workflow
- Download the test data [here](to be added). The current dataset is in a genomic region where there are enough reads to test TAPAS's PAS quantification functionality.

- Update the samplesheet.csv with the full path to the downloaded bam file(s).
```
sample,bam,read_length
sample,[path_to]/SRR6795721.bam,[int(read_length)]
```
- Run the pipeline with the samplesheet.csv with the input paths updated using either docker or singularity containers
```
nextflow main.nf --input samplesheet.csv --gtf [/path/to/gtf] -profile <docker/singularity>
```

## Output & post-processing
TAPAS outputs a file containing TAPAS identified entries, which are formatted into a bed file with the following columns:
```
chrom,chromStart,chromEnd,name,score,strand
```
Please do note that the column names are extracted from the transcript_id attribute in the last column of the gtf file.

## Author contact
If you have any question or comment about TAPAS, please post on TAPAS GitHub Issues (https://github.com/arefeen/TAPAS/issues) or the author, [Dr. Ashraful Arefeen](https://scholar.google.com/citations?user=qaJhymQAAAAJ&hl=en).

