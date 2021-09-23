// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Convert provided gtf genome file to bed
    Create gene symbol file from the gtf genome file
*/
process PREPROCESSING {
    publishDir "${params.outdir}/dapars/genome_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    file genome_file

    output:
    tuple path(gene_symbol_file), path(converted_genome_file), emit: ch_genome_file

    script:
    converted_genome_file = "genemodel.bed"
    gene_symbol_file = "gene_symbol_file.txt"
    """
    gtfToGenePred $genome_file test.genePhred
    genePredToBed test.genePhred $converted_genome_file
    create_gene_symbol_file.py $genome_file $gene_symbol_file
    """
}