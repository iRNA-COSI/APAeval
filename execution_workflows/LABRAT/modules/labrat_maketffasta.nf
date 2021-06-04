// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LABRAT_MAKETFFASTA {
        publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode
        container "quay.io/biocontainers/labrat:0.2.2--pyhdfd78af_0"

        input:
        tuple path(gff), path(fasta)

        output:
        path "*.fasta", emit: ch_tffasta

        script:
        """
        LABRAT.py --mode makeTFfasta --gff $gff --genomefasta $fasta --lasttwoexons --librarytype RNAseq
        """
}
