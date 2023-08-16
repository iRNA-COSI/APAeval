// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def final_output = modules['final_output']

/*
    Publish differential challenge output file in the directory user specified
*/
process PUBLISH_DIFFERENTIAL_OUTPUT{
    publishDir "${params.outdir}/csi_utr/${final_output.output_dir}", mode: params.publish_dir_mode

    input:
    path file

    output:
    path file, emit: ch_differential_out

    script:
    """
    """
}
