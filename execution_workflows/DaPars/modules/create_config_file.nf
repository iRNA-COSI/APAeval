// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['create_config_file']

process CREATE_CONFIG_FILE {
        tag "$sample"
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
        container "docker.io/apaeval/dapars:latest"

        input:
        path annotated_3utr

        output:
        path "*", emit: ch_dapars_input

        script:
        bedgraphs_dir = "$PWD/${params.outdir}/dapars/sample_bedgraph_files_dir"
        output_dir = "$PWD/${params.outdir}/dapars/"
        config_output = "$PWD/${params.outdir}/dapars/config"
        num_least_in_group1 = options.num_least_in_group1
        num_least_in_group2 = options.num_least_in_group2
        coverage_cutoff = options.coverage_cutoff
        fdr_cutoff = options.fdr_cutoff
        pdui_cutoff = options.pdui_cutoff
        fold_change_cutoff = options.fold_change_cutoff

        """
        create_config_file.py \
        $annotated_3utr \
        $bedgraphs_dir \
        $output_dir \
        $num_least_in_group1 \
        $num_least_in_group2 \
        $coverage_cutoff \
        $fdr_cutoff \
        $pdui_cutoff \
        $fold_change_cutoff
        $config_output
        """
 }