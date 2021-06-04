// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QAPA_QUANT {
        tag "$sample"
        publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode
        container  "docker.io/apaeval/qapa:1.3.1"

        input:
        tuple val(sample), path(salmon_quantsf), path(mart_export)

        output:
        tuple val(sample), path("$qapa_results"), emit: ch_qapa_outputs

        script:
        qapa_results   = "qapa_results.txt"
        """
        qapa quant --db $mart_export $salmon_quantsf > $qapa_results
        """
 }
