// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


/*
    Perform STAR alignment
*/
process PREPROCESS_GENOME {
    publishDir "${params.outdir}/isoscm/genome_index_file", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"

    input:
    path gtf_genome_file
    path fasta_genome_file

    output:
    path genome_dir, emit: ch_genome_dir

    script:
    genome_dir = "star_index/"
    """
    STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $genome_dir \
     --genomeFastaFiles $fasta_genome_file \
     --sjdbGTFfile $gtf_genome_file \
     --sjdbOverhang 74
    """
}
