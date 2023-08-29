# GETUTR Nextflow method workflow
This is the initial version of the [GETUTR](http://big.hanyang.ac.kr/GETUTR/manual.htm) method workflow.

## Input
The input is a samplesheet.csv file:

```
sample,bam
sample1,/full/path/sample1.bam
sample2,/full/path/sample2.bam
...
```

For each sample line, the results are stored in `results/$sample.PAVA.cps.2.0.0.bed` and `results/$sample.PAVA.smoothed.2.0.0.bed`. For the Mayr samples, the memory requirements are high, around 128GB for the larger samples (bam file of 8GB).

## Output
GETUTR outputs identification and relative usage quantification bed files.

## Parameters
Parameters used to run GETUTR are specified in conf/modules.config file. Parameters relevant to the workflow are:
- `run_identification` - set to true to obtain identification challenge output.
- `run_relative_usage_quantification` - set to true to obtain relative usage quantification challenge output.
- `identification_out_suffix` - suffix of the output file(s) for the current run ending with .bed when running identification,
                                the prefix will be the different sample names obtained from the sample column in the sample sheet
- `relative_usage_quantification_out_suffix` - suffix of the output file(s) for the current run ending with .bed when running relative usage quantification,
                                the prefix will be the different sample names obtained from the sample column in the sample sheet

## Running

Running the nextflow once you define the `samplesheet.csv` is a one line command:

`nextflow main.nf --input samplesheet.csv --gtf /full/path/genome.gtf -profile docker` to run with docker
or
`nextflow main.nf --input samplesheet.csv --gtf /full/path/genome.gtf -profile singularity` to run with singularity
Voila
