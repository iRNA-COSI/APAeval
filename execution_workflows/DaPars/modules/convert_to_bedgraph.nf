// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert the input .bam file to .bedgraph for DaPars step 2
    Check that the bedgraph file is valid (i.e. has leading 'chr' in the chromosome column)
*/
process CONVERT_TO_BEDGRAPH {
    tag "$sample"
    publishDir "${params.outdir}/dapars/sample_bedgraph_files/${condition}", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    tuple val(sample), val(condition), path(bam_file), path(bai_file)

    output:
    tuple val(sample), val(full_bedgraph_file_path), emit: ch_convert_to_bedgraph_out
    path bedgraph_file, emit: ch_bedgraph_file

    script:
    bedgraph_file = sample + ".bedgraph"
    full_bedgraph_file_path = "$PWD/${params.outdir}/dapars/sample_bedgraph_files/${condition}/" + bedgraph_file
    """
    bedtools genomecov -ibam $bam_file -bg > $bedgraph_file
    check_bedgraph.py $bedgraph_file
    """
}