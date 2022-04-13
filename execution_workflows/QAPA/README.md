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

When using your own data and input file instead of the provided test data and sample sheet, make sure to include in the 
input file you are using the absolute paths to the files, with the column names following the column
names above.

### Parameters
Parameters used to run QAPA are specified in conf/modules.config file. 
Parameters relevant to the workflow itself are:
- `input` - path to the `samplesheet.csv`
- `outdir` - name of the folder that the final output files are going to be in, located under QAPA/
- `gtf`: GTF annotation file
- `fasta`: FASTA reference sequence file
- `polyabed`: BED 3'UTR library
- `run_qapa_build`: run qapa built (default: false)

Notes on `--polyabed` and `--run_qapa_build`:
1. If building annotations from scratch with `qapa build`, please pass the `--run_qapa_build` flag and provide both GTF (`--gtf`) and poly(A) BED file (`--polyabed`)
2. If providing pre-generated QAPA annotations, please pass both GTF (`--gtf`) and poly(A) BED file (`--polyabed`) only (DON"T pass the `--run_qapa_build` flag)

### Running the workflow
The QAPA workflow only does quantification. Once parameters have been set in conf/modules.config file, run the pilot benchmark nextflow pipeline with 
```
nextflow main.nf --input samplesheet_example_files.csv` --gtf <path_to_gtf> --polyabed <path_to_poly(A)_bed> --fasta <path_to_fasta> --run_qapa_build -profile [docker/singularity]
```
#### Running with Docker or Singularity
##### Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running 
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
To run with `docker`, please indicate `-profile docker`
```
nextflow main.nf --input samplesheet_example_files.csv` --gtf <path_to_gtf> --polyabed <path_to_poly(A)_bed> --fasta <path_to_fasta> --run_qapa_build -profile docker
```

##### Singularity
To run with `singularity`, please indicate `-profile singularity`
```
nextflow main.nf --input samplesheet_example_files.csv --gtf <path_to_gtf> --polyabed <path_to_poly(A)_bed> --fasta <path_to_fasta> --run_qapa_build -profile singlularity
```


## Output & post-processing
When using the default output_dir parameter value in conf/modules.config, QAPA store results under
`results/qapa` folder, and the quantification output BED files will be stored in `results/qapa/<sample_name>/qapa_quant.bed`.

## Author contact
If you have any question or comment about QAPA, please submit an issue on [QAPA's GitHub repository](https://github.com/morrislab/qapa/issues)
