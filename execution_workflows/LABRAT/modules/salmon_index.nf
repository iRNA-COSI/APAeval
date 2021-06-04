// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SALMON_INDEX {
        publishDir "${params.outdir}/labrat/salmon", mode: params.publish_dir_mode
        container "quay.io/biocontainers/salmon:0.14.0--h86b0361_1"

        input:
        path txfasta

        output:
        path "*.idx", emit: ch_txfasta_idx

        script:
        """
        salmon index -t $txfasta -i txfasta.idx --type quasi -k 31 --keepDuplicates
        """
}
