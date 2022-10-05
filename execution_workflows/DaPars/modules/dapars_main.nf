// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process DAPARS_MAIN {
        tag "$sample"
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
        container "docker.io/apaeval/dapars:latest"

        input:
        path config_file
        val sample

        output:
        path dapars_output_file, emit: ch_dapars_output

        script:
        dapars_output_file = "dapars_output_All_Prediction_Results.txt"
        """
        python /dapars/src/DaPars_main.py $config_file
        """
 }



