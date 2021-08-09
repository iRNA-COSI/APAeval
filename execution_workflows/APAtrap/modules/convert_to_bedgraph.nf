// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for APAtrap
*/
process CONVERT_TO_BEDGRAPH {

    publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    tuple path(bam_file1), path(bai_file1), path(bam_file2), path(bai_file2)

    output:
    path "*.bedgraph", emit: ch_bedgraph

    script:
    bedgraph_file1 = "3utr_input1.bedgraph"
    bedgraph_file2 = "3utr_input2.bedgraph"
    """
    bedtools genomecov -ibam $bam_file1 -bg > $bedgraph_file1
    bedtools genomecov -ibam $bam_file2 -bg > $bedgraph_file2
    """

}
