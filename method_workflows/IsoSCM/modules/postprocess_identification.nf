// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def output_files = modules['output_files'].clone()
def identification_out_suffix = output_files.identification_out_suffix

process POSTPROCESS_IDENTIFICATION {
        publishDir "${params.outdir}/isoscm/${output_dir}", mode: params.publish_dir_mode
        container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), path(compare_output)

        output:
        path(identification_output_file), emit: ch_identification_postprocess_out

        script:
        output_dir = "${output_files.output_dir}"
        identification_output_file = sample + "_" + identification_out_suffix
        run_mode = 'identification'
 
        """
        postprocess.py $compare_output $run_mode $identification_output_file
        """
 }
