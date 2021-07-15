// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for APAtrap
*/
process PREPROCESSING {

    publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
    container "faricazjj/apatrap"

    input:
    tuple path(bam_file), path(bai_file)

    output:
    path "$input_bedgraph", emit: ch_3utr_input

    script:
    input_bedgraph = "3utr_input.bedgraph"
    """
    bamCoverage -b $bam_file -o $input_bedgraph
    """

}
