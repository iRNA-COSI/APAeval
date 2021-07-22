// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Ensure that the final bedgraph file has leading "chr"
    in all sequence regions
*/
process PREPROCESSING {

    publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path bedgraph_file

    output:
    path "$final_input_bedgraph", emit: ch_3utr_input

    script:
    final_input_bedgraph = "3utr_input_final.bedgraph"
    """
    check_bedgraph.py $bedgraph_file $final_input_bedgraph
    """

}
