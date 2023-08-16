// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Move bam and bai files into a directory to be given as an input for CSI-UTR at a later step
*/
process CREATE_INPUT_FILE {
    publishDir "${params.outdir}/csi_utr/input_files", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(condition), val(replicate), path(bam_file), path(bai_file)

    output:
    tuple path(bam_file), path(bai_file) , emit: ch_input_files

    script:
    """
    """
}
