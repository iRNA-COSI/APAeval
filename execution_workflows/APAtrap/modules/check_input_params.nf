// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Checks that mode and the extension of the file names set in conf/modules.config are valid
*/
process CHECK_INPUT_PARAMS {
    publishDir "${params.outdir}/apatrap/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    script:
    run_identification = options.run_identification
    run_quantification = options.run_quantification
    run_differential = options.run_differential
    identification_out_suffix = options.identification_out_suffix
    quantification_out_suffix = options.quantification_out_suffix
    differential_out = options.differential_out
    """
    check_input_params.py $identification_out_suffix \
    $quantification_out_suffix \
    $differential_out \
    $run_identification \
    $run_quantification \
    $run_differential
    """
}