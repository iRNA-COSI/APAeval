// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Create files for identification and quantification challenges
*/
process CONVERT_TO_BED {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/challenges_outputs", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    tuple val(sample), path(predict_apa_output_file)

    output:
    path "*"

    script:
    identification_out = sample + "_" + inputs.identification_out_suffix
    """
    convert_to_bed.py $predict_apa_output_file $identification_out $quantification_out
    """

}
