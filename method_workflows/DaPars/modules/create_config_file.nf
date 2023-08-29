// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['create_config_file']

process CREATE_CONFIG_FILE {
        container "docker.io/apaeval/dapars:latest"

        input:
        tuple val(sample), val(bedgraph_file), path(annotated_3utr)
        val run_mode

        output:
        path config_output, emit: ch_dapars_input
        val sample, emit: ch_sample

        script:
        annotated_3utr = "$PWD/${params.outdir}/dapars/final_extracted_3utr.bed"
        bedgraphs_dir = "$PWD/${params.outdir}/dapars/sample_bedgraph_files"
        output_dir = "."
        num_least_in_group1 = options.num_least_in_group1
        num_least_in_group2 = options.num_least_in_group2
        coverage_cutoff = options.coverage_cutoff
        fdr_cutoff = options.fdr_cutoff
        pdui_cutoff = options.pdui_cutoff
        fold_change_cutoff = options.fold_change_cutoff
        if ( run_mode == "identification" || run_mode == "relative_usage_quantification" ) {
            config_output = sample + "_config"
        }
        else {
            config_output = "config"
        }

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
        $fold_change_cutoff \
        $config_output \
        $bedgraph_file \
        $run_mode
        """
 }
