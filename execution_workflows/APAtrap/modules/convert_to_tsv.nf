// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Create files for differential challenge
*/
process CONVERT_TO_TSV {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/challenges_outputs", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    tuple val(sample), path(de_apa_output_file)

    output:
    path "*"

    script:
    differential_out = inputs.differential_out
    """
    convert_to_tsv.py $de_apa_output_file $differential_out
    """

}
