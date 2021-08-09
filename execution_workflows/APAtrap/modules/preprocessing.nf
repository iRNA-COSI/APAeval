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
    tuple path(bedgraph_file1), path(bedgraph_file2)

    output:
    path "*.bedgraph", emit: ch_3utr_input

    script:
    final_input_bedgraph1 = "3utr_input_final1.bedgraph"
    final_input_bedgraph2 = "3utr_input_final2.bedgraph"
    """
    check_bedgraph.py $bedgraph_file1 $final_input_bedgraph1
    check_bedgraph.py $bedgraph_file2 $final_input_bedgraph2
    """

}
