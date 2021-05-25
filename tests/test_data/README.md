# Test data for APAeval
Use the files provided here for developing, debugging, testing your code.

* `.bam` files or `.fastq.gz` files as inputs for execution workflows ([Created by randomly choosing 100 genes with MACE-seq PAS reads (GSE151724) and subsetting bam files (generated for the pilot_benchmark: siControl_R1, SRR11918577 and siSrsf3_R1, SRR11918579) to reads contain reads falling within +/- 1kb from gene boundaries from the gtf file])
* corresponding `.gtf` ([Created GENCODE release M18 subsetted to the 100 random genes for the test data])
* `.MACEseq.mm10.bed` as a ground truth example files ([BED6 files for clevage and poly(A) sites for two samples (siControl_R1, SRR11918617 and siSrsf3_R2, SRR11918619) where the score column corresponds to total number of MACE-seq reads overlapping the given PAS in that sample])
* [EXTEND THIS LIST WHEN ADDING MORE TEST FILES TO THE DIRECTORY]
