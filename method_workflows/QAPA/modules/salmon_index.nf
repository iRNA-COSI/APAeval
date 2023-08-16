// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SALMON_INDEX {
        publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode
        container "quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

        input:
        path indexed_fasta

        output:
        path "$utr_library", emit: ch_utr_library

        script:
        utr_library    = "utr_library"
        """
        salmon index -t $indexed_fasta -i $utr_library
        """
}
