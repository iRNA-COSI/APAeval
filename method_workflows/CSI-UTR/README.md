# CSI-UTR
CSI-UTR is a perl-based model that detects differential expression of 3-prime UTRs.

The [paper](https://www.frontiersin.org/articles/10.3389/fgene.2019.00182/full) is titled Detection of Differentially
Expressed Cleavage Site Intervals Within 3' Untranslated Regions Using CSI-UTR Reveals Regulated Interaction Motifs <br>
The [tool](https://github.com/UofLBioinformatics/CSI-UTR) is free to download
and [documentation](https://github.com/UofLBioinformatics/CSI-UTR/blob/master/EXAMPLES/TestCases.md) was used as a reference
to create the nextflow pipeline flow of this module

## Running CSI-UTR workflow


### Data download
CSI-UTR requires exactly two distinct conditions to be provided as inputs. Furthermore,  at least two replicates need to be provided per condition otherwise the tool errors out. The test data provided in APAeval only has one replciate per condition. As such, additional test datasets of two distinct conditions and two replicates per condition are provided [here](https://drive.google.com/drive/folders/16BXeJIencg14K8XU5tiYh_XUQhneq-zg?usp=sharing).

Additionally, CSI-UTR annotation and bed files have to be downloaded. In the same folder above, there are folders called annotations and locations that contain CSI annotation and bed files needed to run the tool for Mm10 genome, the genome used for APAeval test data.

Hence, to run CSI-UTR with APAeval test data, download all the folders and files in the folder above and place them under a data/ folder in CSI-UTR working directory.

CSI-UTR is compatible with the following genomes: Hg38, Mm10, Rn6, and Rn6_extended. For each genome, CSI annotation and bed files are required to be downloaded from [CSI-UTR github repo](https://github.com/UofLBioinformatics/CSI-UTR/tree/master/CSI-UTR_v1.1.0/data).

### Input & pre-processing
An example sample sheet is available at `samplesheet_example_files.csv`. Each row in the samplesheet has four
columns:

- sample: name of the sample to be included as prefix in identification challenge output file name(e.g control_replicate1)
- condition: name of the condition (e.g control)
- replicate: replicate number of this condition (e.g 1)
- bam: BAM input file for the sample
- bai: BAI index file for sample's bam input

It is important to name samples of the same condition with the exact condition name under the condition
column in the samplesheet since samples are grouped per condition to be processed in the differential step.
To run CSI-UTR with test data provided for APAeval, check the path to CSI-UTR with `pwd` and replace
the `path_to` in samplesheet_example_files.csv with the path
from the `pwd` command.

When using your own data and input file instead of the provided test data and sample sheet, make sure to include in the
input file you are using the absolute paths to the files, with the column names following the column
names above.

### Running with Docker or Singularity
## Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
Additionally, make sure that line 49 in CSI-UTR/nextflow.config file `docker.enabled=true` is uncommented while line
51 `singularity.enabled=true` is commented out

## Singularity
To run with singularity, comment out line 49 in Dapars/nextflow.config file `docker.enabled=true` and make sure that line
51 `singularity.enabled=true` is uncommented

### Parameters
Parameters used to run the CSI-UTR are specified in conf/modules.config file.
Parameters relevant to the workflow itself are:
- `run_differential` - true/false, whether to run to obtain differential challenge output.
- `output_dir` - name of the folder that the final output files are going to be in, located under CSI-UTR/results/csi-utr/
- `identification_out_suffx` - suffix of the output file for identification challenge ending with .bed, final output file name is prefixed with sample name obtained from the input sample sheet
- `differential_out` - name of differential challenge output file ending in .tsv
- `differential_type` - which differential expression approach to use and report in the differential challenge output. Must be one of 'DEXSeq', 'PAIRWISE' or 'WITHIN_UTR' (case sensitive). Consult publication for definitions of these approaches
- `genome_file` - absolute path from the CSI-UTR folder to the input GTF annotation file can be obtained by replacing `path_to`
   with the path to CSI-UTR, and if using your own genome file, make sure to use the absolute path to your genome file
- `CSI_bed_file` - absolute path from the CSI-UTR folder to the downloaded CSI bed file can be obtained by replacing `path_to`
   with the path to CSI-UTR
- `CSI_annotation_file` - absolute path from the CSI-UTR folder to the downloaded CSI annotation file can be obtained by replacing `path_to` with the path to CSI-UTR
- `genome` - the genome of the samples, value could be Hg38, Mm10, Rn6, or Rn6_extended
- `r` - read length
- `coverage_cut` - coverage cutoff
- `usage_cut` - usage cutoff
- `p` - p-value significance cutoff
- `q` - FDR significance cutoff


### Running the differential workflow
- Set the 'run_differential' parameter in conf/modules.config to true
- Change 'differential_out' parameter in conf/modules.config to the desired file name that ends with '.tsv'
- Ensure the sample sheet contains exactly two distinct conditions in the condition column and at least two replicates per condition. An example input file
  is samplesheet_example_files.csv
- Run the nextflow pipeline with nextflow main.nf --input samplesheet_example_files.csv

## Output & post-processing
CSI-UTR identification run results in a file located under CSi-UTR/results/csi_utr/[output_dir]/[sample_name]_[identification_out_suffix]
CSI-UTR differential run results in a file located under CSi-UTR/results/csi_utr/[output_dir]/[differential_out]

## Notes
- CSI-UTR uses a pre-provided set of polyA site annotations and does not perform de-novo polyA site identification. As such it does not qualify for the identification challenge.
- It is not possible to obtain quantification challenge output data since the TPM value in the output file
  is provided per CSI region and not per APA site. Hence, this tool is not compatible for quantification
  challenge.
- CSI-UTR needs exactly two distinct conditions with at least 2 replicates per condition. If less than 2 replicates are provided per condition, the tool will error out.

## Author contact
If you have any question or comment about CSI-UTR, please contact the corresponding author, Dr. Eric Rouchka (ecrouc01@louisville.edu).
