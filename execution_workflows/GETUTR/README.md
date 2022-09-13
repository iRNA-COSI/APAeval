# GETUTR Nextflow execution workflow
This is the initial version of the [GETUTR](http://big.hanyang.ac.kr/GETUTR/manual.htm) execution workflow.

## Input and Output
The input is a samplesheet.csv file:

```
sample,bam
sample1,/full/path/sample1.bam
sample2,/full/path/sample2.bam
...
```

For each sample line, the results are stored in `results/$sample.PAVA.cps.2.0.0.bed` and `results/$sample.PAVA.smoothed.2.0.0.bed`. For the Mayr samples, the memory requirements are high, around 128GB for the larger samples (bam file of 8GB).

## Running

Running the nextflow once you define the `samplesheet.csv` is a one line command:

`nextflow main.nf --input samplesheet.csv -profile docker`

Voila