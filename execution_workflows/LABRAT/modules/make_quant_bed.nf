// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MAKE_QUANT_BED {
    tag "$sample"
    publishDir "${params.outdir}/labrat/salmon/$sample", mode: params.publish_dir_mode

    input:
    tuple path(gff), val(sample), path(labrat_salmon_quant)

    output:
    path "*.bed"

    script:
    """
    make_quant_bed.py $labrat_salmon_quant/quant.sf labrat_quant.bed $gff
    """
}
