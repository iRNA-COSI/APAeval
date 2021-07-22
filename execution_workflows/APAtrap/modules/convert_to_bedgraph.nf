// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for APAtrap
*/
process CONVERT_TO_BEDGRAPH {

    publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
    container "faricazjj/apatrap"

    input:
    tuple path(bam_file), path(bai_file)

    output:
    path "$bedgraph_file", emit: ch_bedgraph

    script:
    bedgraph_file = "3utr_input.bedgraph"
    """
    bedtools genomecov -ibam $bam_file -bg > $bedgraph_file
    """

}
