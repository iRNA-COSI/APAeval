// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LABRAT_CALCULATEPSI {
        echo true
        publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode
        container "quay.io/biocontainers/labrat:0.2.2--pyhdfd78af_0"

        input:
        val salmon_dir
        path sampcond
        val conditionA
        val conditionB
        path gff

        output:
        path "*", emit: ch_labrat_diff

        script:
        """
        LABRAT.py --mode calculatepsi --librarytype RNAseq --salmondir $PWD/${params.outdir}/$salmon_dir --sampconds $sampcond --conditionA $conditionA --conditionB $conditionB --gff $gff
        """
}
