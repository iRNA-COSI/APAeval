// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Convert DaPars output file to differential challenge file
*/
process POSTPROCESSING {
    publishDir "${params.outdir}/dapars/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    tuple path(config_file)

    output:
    path "*"

    script:
    mode = options.mode
    output_file_name = options.output_file
    file = "$PWD/${params.outdir}/dapars/dapars_output_All_Prediction_Results.txt"
    """
    convert_output.py $file $output_file_name $mode
    """
}