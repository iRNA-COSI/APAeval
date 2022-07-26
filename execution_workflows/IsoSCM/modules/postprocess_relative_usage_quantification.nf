// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def output_files = modules['output_files'].clone()
def relative_usage_quantification_out_suffix = output_files.relative_usage_quantification_out_suffix

process POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION {
        publishDir "${params.outdir}/isoscm/${output_dir}", mode: params.publish_dir_mode
        container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), path(compare_output)

        output:
        path(relative_usage_quantification_output_file), emit: ch_relative_usage_quantification_postprocess_out

        script:
        output_dir = "${output_files.output_dir}"
        relative_usage_quantification_output_file = sample + "_" + relative_usage_quantification_out_suffix
        run_mode = "relative_usage_quantification"
 
        """
        postprocess.py $compare_output $run_mode $relative_usage_quantification_output_file
        """
 }
