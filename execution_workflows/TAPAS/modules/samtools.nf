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
    tuple val(sample), path(bam), val(read_length)

    output:
    tuple val(sample), path("*.txt"), val(read_length), emit: read_coverage

    script:
    """
    samtools depth $bam > read_coverage.txt    
    """
}
