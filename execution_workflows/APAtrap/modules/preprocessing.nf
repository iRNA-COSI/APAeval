// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Ensure that the final bedgraph file has leading "chr"
    in all sequence regions
*/
process PREPROCESSING {
    publishDir "${params.outdir}/apatrap/genome_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/apatrap:latest"

    input:
    file genome_file

    output:
    path(converted_genome_file), emit: ch_genome_file

    script:
    converted_genome_file = "genemodel.bed"

    """
    gtfToGenePred $genome_file test.genePhred
    genePredToBed test.genePhred $converted_genome_file
    """
}
