// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QAPA_PREPARE_UTRLIB {
        publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode

        input:
        tuple path(fasta), path(bed)

        output:
        path "$utr_library", emit: ch_utr_library

        script:
        indexed_fasta  = "output_sequences.fa"
        utr_library    = "utr_library"
        """
        qapa fasta -f $fasta $bed $indexed_fasta
        salmon index -t $indexed_fasta -i $utr_library
        """
}
