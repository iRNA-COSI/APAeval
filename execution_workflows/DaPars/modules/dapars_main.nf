// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()

process DAPARS_MAIN {
        tag "$sample"
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
        container "docker.io/apaeval/dapars:latest"

        input:
        path config_file
        val sample
        val run_mode

        output:
        path output_file, emit: ch_dapars_output

        script:
        if ( run_mode == "identification" ) {
            output_file = sample + "_dapars_main_out.txt"
        }
        else {
            output_file = "dapars_main_out.txt"
        }

        """
        python /dapars/src/DaPars_main.py $config_file > $output_file
        """
 }



