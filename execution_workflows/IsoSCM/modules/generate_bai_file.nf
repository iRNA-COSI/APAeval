// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Take provided bam file to generate the index bai file
*/
process GENERATE_BAI_FILE {
    publishDir "${params.outdir}/isoscm/aligned_bam_files/${sample}.${strand}", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"

    input:
    tuple val(sample), val(strand), path(bam)

    output:
    tuple val(sample), val(strand), path(bam), path(bai), emit: ch_bam_files

    script:
    """
    samtools index $bam > $bai
    """
}        
