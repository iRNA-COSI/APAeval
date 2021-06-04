// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QAPA_INDEX {

        publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode
        container  "docker.io/apaeval/qapa:1.3.1"

        input:
        tuple path(fasta), path(bed)

        output:
        path "$indexed_fasta", emit: ch_indexed_fasta

        script:
        indexed_fasta  = "output_sequences.fa"
        """
        qapa fasta -f $fasta $bed $indexed_fasta
        """
}
