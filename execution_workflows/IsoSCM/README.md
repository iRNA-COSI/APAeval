# IsoSCM
IsoSCM is a Java application. IsoSCM can be used to estimate the relative usage of 3'UTR isoforms, and identify sites with differential usage between two conditions.
is a python-based model that directly infers the dynamic alternative polyadenylation (APA)
usage by comparing standard RNA-seq. Given the annotated gene model, DaPars can infer the de novo proximal APA sites 
as well as the long and short 3’UTR expression levels. Finally, it identifies the dynamic APA usages between two 
conditions. 

Isoform Structural Change Model (IsoSCM) is a method for transcript assembly that incorporates change-point analysis to improve the 3' UTR annotation process.

The [paper](https://rnajournal.cshlp.org/content/21/1/14) is titled IsoSCM: improved and alternative 3′ UTR annotation 
using multiple change-point inference
The [application download](https://github.com/shenkers/isoscm/releases/tag/IsoSCM-2.0.12) is free to download
The [documentation](https://github.com/shenkers/isoscm) was used as a reference to 
create the nextflow pipeline

## Running IsoSCM workflow

### Input & pre-processing
An example sample sheet is available at `samplesheet_example_files.csv`. Each row in the samplesheet has four
columns:

- sample: name of the sample for logs (e.g control_replicate1)
- condition: name of the condition (e.g Control) 
- bam: BAM input file for the sample 
- bai: BAI index file for sample's bam input

It is important to name samples of the same condition with the exact condition name under the condition
column in the samplesheet since samples are grouped per condition to be processed in the differential step.
To run DaPars with test data provided for APAeval, check the path to DaPars with `pwd` and replace 
the `path_to` in samplesheet_example_files.csv or samplesheet_example_files_identification.csv with the path 
from the `pwd` command. 

When using your own data and input file instead of the provided test data and sample sheet, make sure to include in the 
input file you are using the absolute paths to the four files, with the four column names following the column
names above.

### Running with Docker or Singularity
## Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running 
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
Additionally, make sure that line 49 in IsoSCM/nextflow.config file `docker.enabled=true` is uncommented while line
51 `singularity.enabled=true` is commented out

## Singularity
To run with singularity, comment out line 49 in IsoSCM/nextflow.config file `docker.enabled=true` and make sure that line
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
- Set the 'mode' parameter in conf/modules.config to "differential"
- Change 'output_file' parameter in conf/modules.config to the desired file name that ends with '.tsv'
- Ensure the sample sheet contains exactly two distinct conditions in the condition column. An example input file 
  is samplesheet_example_files.csv
- Run the pilot benchmark nextflow pipeline with nextflow main.nf --input samplesheet_example_files.csv

### Running the identification workflow
- Set the 'mode' parameter in conf/modules.config to "identification"
- Change 'output_file' parameter in conf/modules.config to the desired file name that ends with '.bed'
- In the sample sheet, the same sample should be provided twice, but each row should have a distinct condition in 
  the `condition` column. Exactly two rows must be present in the sample sheet. The workflow will then treat the 
  two rows as two different conditions, a requirement for DaPars to run successfully. An example sample sheet is 
  samplesheet_example_files_identification.csv.
- Run the pilot benchmark nextflow pipeline with nextflow main.nf --input samplesheet_example_files_identification.csv

## Output & post-processing
Each DaPars run results in differential challenge file located under DaPars/results/dapars/dapars_differential_output.tsv.
Identification challenege file is located under DaPars/results/dapars/dapars_identification_output.tsv.

## Notes
- Running STAR requires docker to have 6GB allocated memory

## Author contact
If you have any question or comment about DaPars, please post on DaPars Google Groups (https://groups.google.com/u/1/g/DaPars) or the author, Dr. Zheng Xia (xiaz@ohsu.edu).