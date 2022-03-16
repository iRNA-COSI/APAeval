# IsoSCM
Isoform Structural Change Model (IsoSCM) is a Java application for transcript assembly that incorporates change-point analysis to improve the 3' UTR annotation process.

The [paper](https://rnajournal.cshlp.org/content/21/1/14) is titled IsoSCM: improved and alternative 3′ UTR annotation 
using multiple change-point inference
The [application download](https://github.com/shenkers/isoscm/releases/tag/IsoSCM-2.0.12) is free to download
The [documentation](https://github.com/shenkers/isoscm) was used as a reference to 
create the nextflow pipeline

## Running IsoSCM workflow


### Download genome fasta file
IsoSCM workflow runs STAR genome generate for STAR alignment. This step requires a genome fasta file. To run with test data, run the following to download GRCm38 genome fasta file:
`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.p6.genome.fa.gz`

### Input & pre-processing
An example sample sheet is available at `samplesheet_example_files.csv`. Each row in the samplesheet has four
columns:

- sample: name of the sample (e.g control_replicate1)
- bam: BAM input file for the sample 
- strand: the strandedness of the data

Make sure each sample name is unique.

To run IsoSCM with test data provided for APAeval, check the path to IsoSCM with `pwd` and replace 
the `path_to` in samplesheet_example_files.csv with the path 
from the `pwd` command. 

When using your own data and input file instead of the provided test data and sample sheet, make sure to include in the 
input file you are using the absolute paths to the bam files.

### Running with Docker or Singularity
## Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running 
(e.g. to ensure docker is running, run the command `docker --help` and a help message should be printed).
Additionally, make sure that line 49 in IsoSCM/nextflow.config file `docker.enabled=true` is uncommented while line
52 `singularity.enabled=true` is commented out

## Singularity
To run with singularity, comment out line 49 in IsoSCM/nextflow.config file `docker.enabled=true` and make sure that line
52 `singularity.enabled=true` is uncommented

### Parameters
Parameters used to run the two steps of DaPars are specified in conf/modules.config file. 
Parameters relevant to the workflow itself are:
- `run_star_genome_generate` - if true, the workflow will run STAR genome generate to obtain genome index. This process takes around 20-25 minutes and only needs to run once per genome. Set to false if the genome index folder already exists.
- `output_dir` - name of the folder that the final output files are going to be in, located under Dapars/results/dapars/
- `identification_out_suffix` - suffix of the output file for the current run ending with .bed. This will be prefixed with the sample name specified in the sample column of the input sample sheet e.g. samplesheet_example_files.csv
- `gtf_genome_file` - absolute path to the input GTF annotation file can be obtained by replacing `path_to` with the path to IsoSCM by doing `pwd`, and if using your own genome file, make sure to use the absolute path to your genome file
- `fasta_genome_file` - absolute path to the input fasta annotation file. If running with APAeval test data, download the fasta genome file by running `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.p6.genome.fa.gz` inside of tests/test_data/ folder.
- `star_genome_index` - absolute path to the star genome index folder

### Running the identification workflow
- Change 'identification_out_suffix' parameter in conf/modules.config to the desired file name that ends with '.bed'
- Change the 'output_dir' to the desired directory name under IsoSCM/results/isoscm the identification challenge output will be in  
- An example sample sheet is samplesheet_example_files_identification.csv.
- Run the pilot benchmark nextflow pipeline with nextflow main.nf --input samplesheet_example_files.csv

## Output & post-processing
By default, each IsoSCM run results in an identification challenge file located under IsoSCM/results/isoscm/challenges_outputs/[sample]_identification_output.bed

## Notes
- Running STAR genome generate requires 30-40GB of memory
- Running STAR genome generate takes 20-25 minutes and is only required once per genome
- IsoSCM doesn't qualify for quantification challenge because IsoSCM assemble step outputs mean read density in the segments upstream and downstream of the polyA site but doesn't output the TPM value for each identified site that is required for quantification challengenge output  
- IsoSCM doesn't qualify for differential challnege because the compare step provides a confidence score representing the likelihood that a change-point occurs at the given position, compared to a null model that no change-point at this position. However, this doesn't exactly equal to the significance of differential polyadenylation site usage, which is needed for differential challenge output

## Author contact
If you have any question or comment about IsoSCM, contact the corresponding author, Dr. Eric Lai(laie@mskcc.org).
