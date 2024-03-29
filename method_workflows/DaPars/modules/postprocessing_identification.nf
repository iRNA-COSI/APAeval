// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Convert DaPars output file to identification challenge file
*/
process POSTPROCESSING_IDENTIFICATION {
    publishDir "${params.outdir}/dapars/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    val sample
    path dapars_output_file

    output:
    path "*"

    script:
    run_mode = "identification"
    identification_out = "${sample}_${options.identification_out_suffix}"
    """
    convert_output.py $dapars_output_file $identification_out $run_mode
    """
}
