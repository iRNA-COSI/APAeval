// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Convert DaPars output file to differential challenge file
*/
process POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION {
    publishDir "${params.outdir}/dapars/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    val sample
    path dapars_output_file

    output:
    path "*"

    script:
    run_mode = "relative_usage_quantification"
    relative_usage_quantification_out = "${sample}_${options.relative_usage_quantification_out_suffix}"
    """
    convert_output.py $dapars_output_file $relative_usage_quantification_out $run_mode
    """
}
