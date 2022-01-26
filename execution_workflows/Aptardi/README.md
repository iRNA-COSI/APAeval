# Aptardi
Aptardi is a python-based software that uses gtf/gff-formatted transcript features, which can be constructed from RNA-seq data using StringTie2, and RNA-seq data mapped to the genome (DNA) sequences to train a machine-learning model that identifies the 3' ends of transcripts.

The [paper](https://www.nature.com/articles/s41467-021-21894-x) is titled Aptardi predicts polyadenylation sites in sample-specific transcriptomes using high-throughput RNA sequencing and DNA sequence <br>
The [application](https://github.com/luskry/aptardi) is free to download,
and the [README documentation](https://github.com/luskry/aptardi#aptardi) was used as a reference
to create the nextflow pipeline of this module.

This workflow qualifies for the APAeval **identification challenge only**.

## Running Aptardi workflow

### Input & pre-processing
An example sample sheet is available at `samplesheet.csv`. `samplesheet.csv` can contain multiple entries, where each row in the samplesheet has two
columns:

- sample: name of the sample for logs (e.g control_replicate1)
- bam: BAM input file for the sample

### Docker and Singularity containers
This workflow uses docker containers. To run, make sure that docker is installed and running
(e.g. by running the command `docker --help` and seeing a help message printed)
If running with Singularity, please include the `-profile singularity` in the command, which enables Singularity.

### Parameters
Parameters used to run the Aptardi are specified in the nextflow.config file.
Parameters relevant to the workflow itself are:
- `input` - samplesheet.csv
- `fasta` - reference_genome_sequence.fasta file
- `gtf` - reference_annotation.gtf file
- `outdir` - name of the folder that the final output files are going to be in, located under ./results/aptardi/
- `output_bed` - name of the output file for the current run ending with .bed
- `use_stringtie2_gtf` - use the gtf output from stringtie2 (default: false)
- `aptardi_model` - for running Aptardi using a pre-built model (default: ./aptardi_default_model_scale/model.hdf5 provided by Aptardi.)
- `aptardi_scale` - for running Aptardi using a pre-built model (default: ./aptardi_default_model_scale/scale.pk provided by Aptardi)
Note: Aptardi also allows users to build their own model during the Aptardi run, but this pipeline does not support building one's own model (one would have to build his/her own model by running Aptardi outside of this pipeline). Please refer to [Aptardi's documentation](https://github.com/luskry/aptardi#options)

### Running the Aptardi execution workflow
- Download the test data [here](https://drive.google.com/drive/folders/1tsDu7TzxoVvnD-0UbVRd-pu-ZL36F190?usp=sharing). The current dataset is in a genomic region where there are enough reads to test Aptardi's PAS identification functionality. If you wish to test the StringTie2 workflow (`use_stringtie2_gtf` set to true), download the test data linked in the [notes section](#notes).

- Update the samplesheet.csv with the full path to the downloaded bam, gtf, and fasta files.
```
sample,bam
sample,[path_to]/test.bam
```
- Run the pipeline with the samplesheet.csv with the input paths updated using either docker or singularity containers
```
nextflow main.nf --input samplesheet.csv --fasta [fasta file] --gtf [gtf file] -profile <docker/singularity>
```

## Output & post-processing
Aptardi outputs a gtf file containing both the Aptardi identified entries and the entries from the inputted gtf files. Hence, the post-processing of the Aptardi results filters out the gtf entries from the reference gtf file and only selects entries that indicates aptardi as their source. The Aptardi entries are formatted into a bed file with the following columns:
```
chrom,chromStart,chromEnd,name,score,strand
```
Please do note that the column names are extracted from the transcript_id attribute in the last column of the gtf file.

## Notes

- Test data for the StringTie2 workflow can be downloaded from [here](https://drive.google.com/drive/folders/1KAWKU_7saa8fHQQ3smARxL6Py-UXac7G?usp=sharing) (Aptardi predicts no new poly(A) sites with the test data provided in the main section). You will need to unzip the FASTA & GTF file before inputting them to the workflow.

## Author contact
If you have any question or comment about Aptardi, please post on Aptardi GitHub Issues (https://github.com/luskry/aptardi/issues) or the author, Dr. Ryan Lusk (ryan.lusk@cuanschutz.edu).
