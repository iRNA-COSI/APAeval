// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['output_files']

/*
    Convert IsoSCM compare step output file to differential challenge file
*/
process POSTPROCESSING_DIFFERENTIAL {
    publishDir "${params.outdir}/isoscm/${options.output_dir}", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"

    input:
    path isoscm_compare_output

    output:
    path "*"

    script:
    output_file_name = options.differential_file_out
    file = "$PWD/${params.outdir}/isoscm/isoscm_compare_output.txt"
    """
    postprocess_differential.py $file $output_file_name
    """
}