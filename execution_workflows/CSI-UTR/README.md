# CSI-UTR
CSI-UTR is a perl-based model that detects differential expression of 3-prime UTRs.

The [paper](https://www.frontiersin.org/articles/10.3389/fgene.2019.00182/full) is titled Detection of Differentially
Expressed Cleavage Site Intervals Within 3' Untranslated Regions Using CSI-UTR Reveals Regulated Interaction Motifs <br>
The [tool](https://github.com/UofLBioinformatics/CSI-UTR) is free to download
and [documentation](https://github.com/UofLBioinformatics/CSI-UTR/blob/master/EXAMPLES/TestCases.md) was used as a reference
to create the nextflow pipeline flow of this module

## Running CSI-UTR workflow


### Data download
CSI-UTR requires provided in [CSI-UTR github repo](https://github.com/UofLBioinformatics/CSI-UTR/tree/master/CSI-UTR_v1.1.0/data)

CSI-UTR also requires exactly two distinct conditions to be provided as inputs. Furthermore,  at least two replicates need to be provided per condition otherwise the tool errors out. The test data provided in APAeval only has one replciate per condition. As such, additional test datasets of two distinct conditions and two replicates per condition are provided [here]()
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
Additionally, make sure that line 49 in Dapars/nextflow.config file `docker.enabled=true` is uncommented while line
51 `singularity.enabled=true` is commented out

## Singularity
To run with singularity, comment out line 49 in Dapars/nextflow.config file `docker.enabled=true` and make sure that line
51 `singularity.enabled=true` is uncommented

### Parameters
Parameters used to run the two steps of DaPars are specified in conf/modules.config file. 
Parameters relevant to the workflow itself are:
- `mode` - whether to run to obtain identification ("identification") or differential ("differential") challenge output.
   Specifying any other value will throw an error.
- `output_dir` - name of the folder that the final output files are going to be in, located under Dapars/results/dapars/
- `output_file` - name of the output file for the current run ending with .bed if running identification and .tsv if running differential
- `genome_file` - absolute path from the DaPars folder to the input GTF annotation file can be obtained by replacing `path_to`
   with the path to DaPars, and if using your own genome file, make sure to use the absolute path to your genome file

### Running the differential workflow
- Set the 'run_differential' parameter in conf/modules.config to true
- Change 'differential_out' parameter in conf/modules.config to the desired file name that ends with '.tsv'
- Ensure the sample sheet contains exactly two distinct conditions in the condition column and at least two replicates per condition. An example input file 
  is samplesheet_example_files.csv
- Run the pilot benchmark nextflow pipeline with nextflow main.nf --input samplesheet_example_files.csv

### Running the identification workflow
- Set the 'run_identification' parameter in conf/modules.config to true
- Change 'identification_out_suffix' parameter in conf/modules.config to the desired file name that ends with '.bed'
- Ensure the sample sheet contains exactly two distinct conditions in the condition column and at least two replicates per condition. An example input file
  is samplesheet_example_files.csv
- Run the pilot benchmark nextflow pipeline with nextflow main.nf --input samplesheet_example_files_identification.csv

## Output & post-processing
Each DaPars run results in differential challenge file located under DaPars/results/dapars/dapars_differential_output.tsv.
Identification challenege file is located under DaPars/results/dapars/dapars_identification_output.tsv.

## Notes
- It is not possible to obtain quantification challenge output data since the TPM value in the output file
  is provided per CSI region and not per APA site. Hence, this tool is not compatible for quantification
  challenge. 
- CSI-UTR needs exactly two distinct conditions with at least 2 replicates per condition. If less than 2 replicates are provided per condition, the tool will error out.

## Author contact
If you have any question or comment about CSI-UTR, please contact the corresponding author, Dr. Eric Rouchka (ecrouc01@louisville.edu).
