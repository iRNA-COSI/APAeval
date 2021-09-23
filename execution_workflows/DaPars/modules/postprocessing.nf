// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert DaPars output file to differential challenge file
*/
process POSTPROCESSING {
    publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    tuple path(config_file), val(mode)

    output:
    path "*"

    script:
    file = "$PWD/${params.outdir}/dapars/dapars_output_All_Prediction_Results.txt"
    """
    convert_output.py $file $mode
    """
}