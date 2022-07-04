// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MAKE_QUANT_BED {
    tag "$sample"
    publishDir "${params.outdir}/tapas/tapas_quant/", mode: params.publish_dir_mode

    container "quay.io/biocontainers/python:3.8.3"

    input:
    path tapas_quant_txt

    output:
    path "*.bed"

    script:
    """
    make_quant_bed.py $tapas_quant_txt tapas_quant.bed
    """
}
