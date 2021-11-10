// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Convert provided gtf genome file to bed
*/
process CONVERT_GTF_TO_BED {
    publishDir "${params.outdir}/csi_utr/genome_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/csi-utr:latest"

    input:
    file genome_file

    output:
    path converted_genome_file, emit: ch_genome_file

    script:
    converted_genome_file = "./data/annotations/Mm10.CSIs.annot.bed"
    """
    gtfToGenePred $genome_file test.genePred
    genePredToBed test.genePred $converted_genome_file
    """
}