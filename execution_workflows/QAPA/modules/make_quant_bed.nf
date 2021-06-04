// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MAKE_QUANT_BED {
        tag "$sample"
        publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode
        container "quay.io/biocontainers/python:3.8.3"

        input:
        tuple val(sample), path(qapa_results)

        output:
        path "*", emit: ch_qapa_quant_bed

        script:
        qapa_quant_bed = "qapa_quant.bed"
        """
        make_quant_bed.py $qapa_results $qapa_quant_bed
        """
 }
