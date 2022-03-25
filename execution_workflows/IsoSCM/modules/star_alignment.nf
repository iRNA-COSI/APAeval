// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Perform STAR alignment containing three steps
    1. Get the fastq file of a sample bam file
    2. Run STAR alignment of the sample bam file
    3. Get the bai file of the new aligned bam file
*/
process STAR_ALIGNMENT {
    container "docker.io/apaeval/isoscm:latest"
    label 'process_high'

    input:
    tuple val(sample), val(strand), val(read_type), path(fastq1), path(fastq2), path(star_index_file), val(star_genome_index_indicator)
    
    output:
    tuple val(sample), val(strand), path(star_out_bam), emit: ch_aligned_bam_files

    script:
    star_out_bam = sample + ".Aligned.sortedByCoord.out.bam"
    if (read_type == "single") {
    	"""
   	STAR \
     	--runThreadN $task.cpus \
      	--outSAMtype BAM SortedByCoordinate \
     	--genomeDir $star_index_file \
      	--outSAMstrandField intronMotif \
      	--readFilesIn $fastq1 \
        --readFilesCommand cat \
        --outFileNamePrefix $sample. 
        """
    }
    else {
	"""
        STAR \
        --runThreadN $task.cpus \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir $star_index_file \
        --outSAMstrandField intronMotif \
        --readFilesIn $fastq1 $fastq2 \
        --readFilesCommand cat \
        --outFileNamePrefix $sample.
        """
    }
}
