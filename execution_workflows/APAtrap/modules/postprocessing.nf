// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Create files for identification, quantification, and differential challenges
*/
process POSTPROCESSING {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/$sample", mode: params.publish_dir_mode
    container "docker.io/faricazjj/apatrap:latest"

    input:
    tuple val(sample), path(de_apa_output_file)

    output:
    path "*"

    script:
    """
    postprocessing.py $de_apa_output_file
    """

}
