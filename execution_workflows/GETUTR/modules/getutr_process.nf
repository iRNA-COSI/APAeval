// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GETUTR_PROCESS {
        publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode
        container "quay.io/biocontainers/labrat:0.2.2--pyhdfd78af_0"

        label "process_long"

        input:
        tuple path(gff), path(fasta)

        output:
        path "*.fasta", emit: ch_tffasta
        path "*.db", emit: ch_gffutilsdb

        script:
        """
        python2.7 GETUTR/getutr.py --mode makeTFfasta --gff $gff --genomefasta $fasta --lasttwoexons --librarytype RNAseq
        """
}
