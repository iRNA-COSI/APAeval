// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for APAtrap
*/
process CONVERT_TO_BEDGRAPH {
    tag "$sample"
    publishDir "${params.outdir}/apatrap/sample_bedgraph_files/${condition}", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    tuple val(condition), val(sample), path(bam_file), path(bai_file)

    output:
    val "apatrap/sample_bedgraph_files", emit: ch_sample_bedgraph_files_dir
    tuple val(sample), path(bedgraph_file), emit: ch_3utr_input

    script:
    bedgraph_file = "preprocess_input_" + sample + ".bedgraph"
    """
    bedtools genomecov -ibam $bam_file -bg > $bedgraph_file
    """

}
