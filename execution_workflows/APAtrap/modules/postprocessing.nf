// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Create files for identification, quantification, and differential challenges
*/
process POSTPROCESSING {
    tag "$sample"
    publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(sample), path(de_apa_output_file)

    output:
    path "*"

    script:
    """
    postprocessing.py $de_apa_output_file
    """

}
