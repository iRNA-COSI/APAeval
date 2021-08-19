// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
// get the configs for this process
def config_file_options    = modules['create_config_file']

process CREATE_CONFIG_FILE {
        tag "$sample"
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
        container "docker.io/apaeval/dapars:latest"

        output:
        path "*", emit: ch_dapars_input

        script:

        """
        create_config_file.py \
        $annotated_3utr \
        $bedgraphs_dir \
        $output_file \
        $num_least_in_group_1 \
        $num_least_in_group_2 \
        $coverage_cutoff \
        $fdr_cutoff \
        $pdui_cutoff \
        $fold_change_cutoff
        """
 }



PDUI_cutoff=0.5

Fold_change_cutoff=0.59