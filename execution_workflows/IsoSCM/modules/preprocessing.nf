// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Perform STAR alignment
*/
process PREPROCESSING {
    publishDir "${params.outdir}/isoscm/genome_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"

    input:
    tuple val(sample), path(genome_file), path(bam)

    output:
    path converted_genome_file, emit: ch_genome_file

    script:
    converted_genome_file = "genemodel.bed"
    """
    samtools view $bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > '$sample.fastq'
STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir ./star_index/ \
     --outSAMstrandField intronMotif  --readFilesIn siControl_R1_2genes.fastq \
     --readFilesCommand cat --outFileNamePrefix new/new_siControl_R1_2genes.
STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir ./star_index/ \
     --outSAMstrandField intronMotif  --readFilesIn siSrsf3_R1_2genes.fastq \
     --readFilesCommand cat --outFileNamePrefix new/new_siSrsf3_R1_2genes.

samtools index new/new_siControl_R1_2genes.Aligned.sortedByCoord.out.bam
samtools index new/new_siSrsf3_R1_2genes.Aligned.sortedByCoord.out.bam
"""
}