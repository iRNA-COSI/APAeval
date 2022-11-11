// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
def inputs = modules['final_outputs']

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
    relative_quantification_output = sample + inputs.relative_usage_quantification_out_suffix
    """
    postprocess_relative_quantification.py ${getutr_output} ${relative_quantification_output}
    """
}
