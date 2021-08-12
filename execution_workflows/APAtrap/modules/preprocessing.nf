// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Ensure that the final bedgraph file has leading "chr"
    in all sequence regions
*/
process PREPROCESSING {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/sample_bedgraph_files", mode: params.publish_dir_mode
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(sample), path(bedgraph_file)

    output:
    val "apatrap/sample_bedgraph_files", emit: ch_sample_bedgraph_files_dir
    tuple val(sample), path(sample_bedgraph), emit: ch_3utr_input

    script:
    sample_bedgraph = "3utr_input_" + sample + ".bedgraph"
    """
    check_bedgraph.py $bedgraph_file $sample_bedgraph
    """
}
