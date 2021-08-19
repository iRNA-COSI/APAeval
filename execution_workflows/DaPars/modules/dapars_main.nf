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
        path(config_file)

        output:
        path "*", emit: ch_dapars_output

        script:
        """
        DaPars_main.py $config_file
        """
 }



