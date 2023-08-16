// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Convert DaPars output file to differential challenge file
*/
process POSTPROCESSING_DIFFERENTIAL {
    publishDir "${params.outdir}/dapars/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    path dapars_output_file
    val sample

    output:
    path "*"

    script:
    run_mode = "differential"
    differential_out = options.differential_out
    """
    convert_output.py $dapars_output_file $differential_out $run_mode
    """
}
