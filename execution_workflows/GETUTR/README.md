# GETUTR Nextflow execution workflow
This is the initial version of the [GETUTR](http://big.hanyang.ac.kr/GETUTR/manual.htm) execution workflow.

## Input and Output
The input is a samplesheet.csv file, with this structure:

```
sample,bam,gtf
my_sample,/path/to/my_bam.bam,/path/to/my_gtf.gtf
...
```

For each sample line defined above, the results are stored in `results/$sample.PAVA.cps.2.0.0.bed` and `results/$sample.PAVA.smoothed.2.0.0.bed`. For the Mayr samples, the memory requirements are high, around 128GB for the larger samples (bam file of 8GB).

## Running

Running the nextflow once you define the `samplesheet.csv` is a one line command:

`nextflow main.nf --input samplesheet.csv -profile docker`

Voila