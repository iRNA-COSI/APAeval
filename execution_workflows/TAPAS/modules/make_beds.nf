// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process MAKE_BEDS {
    publishDir "${params.outdir}/${params.modules.make_beds.publish_dir}", mode: params.publish_dir_mode
    tag "$sample"

    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(sample), path(tapas_quant_txt)

    output:
    path "*.bed"

    script:
    identification_bed= "$sample" + "${params.identification_bed_suffix}"
    quantification_bed= "$sample" + "${params.quantification_bed_suffix}"
    """
    make_tapas_beds.py $tapas_quant_txt $quantification_bed $identification_bed
    """
}
