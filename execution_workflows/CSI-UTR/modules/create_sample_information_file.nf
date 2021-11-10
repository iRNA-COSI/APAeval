// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Move bam and bai files into a directory to be given as an input for CSI-UTR at a later step
*/
process CREATE_SAMPLE_INFORMATION_FILE {
    publishDir "${params.outdir}/csi_utr", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path sample_information_file , emit: ch_sample_information_file

    script:
    sample_information_file = "sampleInformation.txt"
    """
    create_sample_information_file.py $samplesheet $sample_information_file
    """
}