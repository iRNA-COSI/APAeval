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
    publishDir "${params.outdir}/isoscm/aligned_bam_files/${sample}.${strand}", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"
    label 'process_high'

    input:
    tuple val(sample), path(bam), val(strand), path(star_index_file), val(star_genome_index_indicator)
    
    output:
    tuple val(sample), val(strand), path(star_out_bam), path(star_out_bai), emit: ch_aligned_bam_files

    script:
    fastq_file = sample + ".fastq"
    star_out_bam = sample + ".Aligned.sortedByCoord.out.bam"
    star_out_bai = sample + ".Aligned.sortedByCoord.out.bam.bai"
    """
    samtools bam2fq $bam > $fastq_file
    
    [ -d star_index ] || mkdir star_index 
    
    STAR \\
      --runThreadN $task.cpus \\
      --outSAMtype BAM SortedByCoordinate \\
      --genomeDir $star_index_file \\
      --outSAMstrandField intronMotif \\
      --readFilesIn $fastq_file \\
      --readFilesCommand cat \\
      --outFileNamePrefix bam_files/$sample. 

    cat bam_files/$star_out_bam > $star_out_bam
    samtools index $star_out_bam > $star_out_bai
    """
}
