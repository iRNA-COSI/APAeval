#Background
When running APAtrap on the entire sample, we noticed that 
the workflow used a lot of memory and went for over 70 hours 
before we killed it. As such, we have been running 
APAtrap by splitting the samples by chromosome and running them in parallel 
so that the workflow wouldn't exceed the maximum memory and would complete 
in a reasonable amount of time.

Here, we provided example scripts that are used to preprocess and postprocess the samples 
to allow for the splitting of datasets, creation of input sample sheets, 
and finally the combination of results. We hope that the scripts could be 
used as examples on how to run APAtrap on samples per chromosome in 
case users encounter similar issues with huge memory usage and long 
running time.

# Split BAM files by chormosome
- chr_list.txt - List of chromosomes to split BAM files to
- split_bams_by_chr.py - Python script to generate final shell scripts that split a sample by chromosome for all samples
- split_mayr_bams.sh - Example of one of the final shell scripts outputted by split_bams_by_chr.py to split bam files by chromosome for mayr sample

# Prepare samplesheets per chromosome to run APAtrap execution workflow
- generate_samplesheet.py - generate sample sheets for all chromosomes, following the file name format of the sample files that have been split by chromosome

# Combine challenge outputs per chromosome into one file
- combine_outputs.py - for each sample, combine the per chromosome outputs into one final output file
