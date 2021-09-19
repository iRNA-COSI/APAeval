// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/samtools:1.13--h8c37831_0"

    input:
    tuple val(sample), path(bam), path(gtf)

    output:
    tuple val(sample), path("*.txt"), emit: read_coverage

    script:
    sorted_bam="sorted_"+"$bam"
    accepted_reads_bam="accepted_reads.bam"
    """
    samtools sort $bam -o $sorted_bam
    samtools index -b $sorted_bam
    samtools view -b $sorted_bam > $accepted_reads_bam
    samtools depth $accepted_reads_bam > read_coverage.txt    
    """
}
