// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process

/*
    Convert GETUTR output file to relative quantification challenge file
*/
process POSTPROCESSING_RELATIVE_QUANTIFICATION {
    publishDir "${params.outdir}/getutr", mode: params.publish_dir_mode
    container "docker.io/apaeval/getutr:latest"

    input:
    tuple val(sample), path(getutr_output)

    output:
    path "*"

    script:
    relative_quantification_output = "${sample}_relative_quantification_output.bed"
    """
    postprocess_relative_quantification.py ${getutr_output} ${relative_quantification_output}
    """
}
