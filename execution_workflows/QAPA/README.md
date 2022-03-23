# QAPA
QAPA is a python-and-R-based model for quantifying alternative polyadenylation (APA) from RNA-seq data.

The [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1414-4) is titled QAPA: a new method for the systematic analysis of alternative polyadenylation from RNA-seq data <br>
The [application download](https://github.com/morrislab/qapa) is free to download
and [documentation](https://github.com/morrislab/qapa#rna-seq-quantification-of-alternative-polyadenylation-qapa) was used as a reference
to create the nextflow pipeline flow of this module

## Running the QAPA workflow

### Input & pre-processing
An example sample sheet is available at `samplesheet_example.csv`. Each row in the samplesheet has nine
columns:

- sample: name of the sample for logs (e.g control_replicate1)
- fastq1: FASTQ1 for both single-end and paired-end RNA-seq
- fastq2: FASTQ2 for paired-end RNA-seq. If single-end, please leave this empty
- bam: BAM input file for the sample 
- bai: BAI index file for sample's bam input
- gff: GFF annotation file
- fasta: FASTA reference sequence file
- bed: BED 3'UTR library 
- mart export: Ensembl gene metadata table from Biomart. Please see [QAPA's GitHub page for more details](https://github.com/morrislab/qapa#prepare-annotation-files)

To run QAPA with test data provided for APAeval, check the path to QAPA with `pwd` and replace 
the `path_to` in samplesheet_example.csv.

When using your own data and input file instead of the provided test data and sample sheet, make sure to include in the 
input file you are using the absolute paths to the files, with the column names following the column
names above.

### Running with Docker or Singularity
## Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running 
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
To run with `docker`, please indicate `-profile docker`
```
nextflow main.nf --input samplesheet_example_files.csv` -profile docker
```

## Singularity
To run with `singularity`, please indicate `-profile singularity`
```
nextflow main.nf --input samplesheet_example_files.csv` -profile singlularity
```

### Parameters
Parameters used to run the two steps of DaPars are specified in conf/modules.config file. 
Parameters relevant to the workflow itself are:
- `input` - path to the `samplesheet.csv`
- `outdir` - name of the folder that the final output files are going to be in, located under Dapars/results/dapars/

### Running the workflow
The QAPA workflow only does quantification. Once parameters have been set in conf/modules.config file, run the pilot benchmark nextflow pipeline with 
```
nextflow main.nf --input samplesheet_example_files.csv` -profile [docker/singularity]
```

## Output & post-processing
When using the default output_dir parameter value in conf/modules.config, QAPA store results under
`QAPA/results/` folder.

## Author contact
If you have any question or comment about QAPA, please submit an issue on [QAPA's GitHub repository](https://github.com/morrislab/qapa/issues)
