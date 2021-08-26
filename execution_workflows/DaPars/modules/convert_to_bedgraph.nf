// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for DaPars step 2
*/
process CONVERT_TO_BEDGRAPH {
    tag "$sample"
    publishDir "${params.outdir}/dapars/sample_bedgraph_files/${condition}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    tuple val(sample), val(condition), path(bam_file), path(bai_file)

    output:
    tuple val(sample), path(bedgraph_file), emit: ch_convert_to_bedgraph_out

    script:
    bedgraph_file = sample + ".bedgraph"
    """
    bedtools genomecov -ibam $bam_file -bg > "sample.bedgraph"
    check_bedgraph.py "sample.bedgraph" $bedgraph_file
    """
}