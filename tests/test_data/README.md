# Test data for APAeval
Use the files provided here for developing, debugging, testing your code.

* `.bam` files with index files (`.bam.bai`) as inputs for execution workflows ([Created by choosing 2 genes with MACE-seq PAS reads at at least 2 sites (GSE151724) and subsetting bam files (generated for the pilot_benchmark: siControl_R1, SRR11918577 and siSrsf3_R1, SRR11918579) to reads contain reads falling within +/- 1kb from gene boundaries from the gtf file]).
* `.fastq` files can be generated to test alignments using `samtools bam2fq input.bam > output.fastq`. 
* corresponding `.gtf` ([Created GENCODE release M18 subsetted to the 2 genes for the test data with leading "chr" removed to match bam files])
* corresponding `.gff3` ([Created GENCODE release M18 subsetted to the 2 genes for the test data with leading "chr" removed to match bam files])
* `.MACEseq.mm10.bed` as a ground truth example files ([BED6 files for clevage and poly(A) sites for the 2 genes from the two samples (siControl_R1, SRR11918617 and siSrsf3_R2, SRR11918619) where the score column corresponds to the TPM for each PAS detected by MACE-seq in that sample])
* [EXTEND THIS LIST WHEN ADDING MORE TEST FILES TO THE DIRECTORY]
