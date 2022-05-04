// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QAPA_BUILD_REF {

        publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode
        container  "docker.io/apaeval/qapa:1.3.1"

        input:
        path mart_export
        path genepred
        path polyabed

        output:
        path "*.bed", emit: bed

        script:
        output_utrs_bed  = "output_utrs.bed"
        """
        qapa build --db $mart_export -o $polyabed $genepred > $output_utrs_bed
        """
}
