// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def file_names = modules['postprocessing']

/*
    Create files for identification, quantification, and differential challenges
*/
process POSTPROCESSING {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/$sample", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    tuple val(sample), path(de_apa_output_file)

    output:
    path "*"

    script:
    identification_out = file_names.identification_out
    quantification_out = file_names.quantification_out
    differential_out = file_names.differential_out
    """
    postprocessing.py $de_apa_output_file $identification_out $quantification_out $differential_out
    """

}
