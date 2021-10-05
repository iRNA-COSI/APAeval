// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Create gene symbol file from the gtf genome file
*/
process CREATE_GENE_SYMBOL_FILE {
    publishDir "${params.outdir}/dapars/gene_symbol_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/dapars:latest"

    input:
    file genome_file

    output:
    tuple path(gene_symbol_file), emit: ch_gene_symbol_file

    script:
    gene_symbol_file = "gene_symbol_file.txt"
    """
    create_gene_symbol_file.py $genome_file $gene_symbol_file
    """
}