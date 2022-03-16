// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
def options    = initOptions(params.options)


/*
    Perform STAR genome generate to get genome index for STAR alignment in the next step
*/
process STAR_GENOME_GENERATE {
    publishDir "${params.outdir}/isoscm", mode: params.publish_dir_mode
    container "docker.io/apaeval/isoscm:latest"
    label 'process_high'

    input:
    path gtf_genome_file
    path fasta_genome_file

    output:
    path(star_index_dir), emit: ch_star_index_dir 
    val "done", emit: ch_star_genome_indicator

    script:
    star_index_dir = "star_index"

    """
    [ -d $star_index_dir ] || mkdir $star_index_dir

    STAR \
      --runMode genomeGenerate \
      --runThreadN $task.cpus \
      --genomeDir $star_index_dir \
      --genomeFastaFiles $fasta_genome_file \
      --sjdbGTFfile $gtf_genome_file \
      --sjdbOverhang 74
    """
}
