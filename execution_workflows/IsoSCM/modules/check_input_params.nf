// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def run_mode = modules['run_mode']
def output_files = modules['output_files']

/*
    Check IsoSCM run mode and file name input parameters
*/
process CHECK_INPUT_PARAMS {
    container "docker.io/apaeval/isoscm:latest"

    script:
    run_identification = run_mode.run_identification
    run_relative_usage_quantification = run_mode.run_relative_usage_quantification
    identification_out_suffix = output_files.identification_out_suffix
    relative_usage_quantification_out_suffix = output_files.relative_usage_quantification_out_suffix
    """
    check_input_params.py $identification_out_suffix $relative_usage_quantification_out_suffix $run_identification $run_relative_usage_quantification
    """
}
