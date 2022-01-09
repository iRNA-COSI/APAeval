// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LABRAT_CALCULATEPSI {
        publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode
        container "quay.io/biocontainers/labrat:0.2.2--pyhdfd78af_0"

        cpus 3
        memory '16 GB'

        input:
        val salmon_dir
        path sampcond
        val conditionA
        val conditionB
        path gff
        path db
        val bed

        output:
        path "*.psis" , emit: ch_psis
        path "*.txt"
        path "*.pval" , emit: ch_pval_results

        script:
        def add_to_path = (params.run_quantification) ? "$PWD/${params.outdir}" : ""
        """
        LABRAT.py --mode calculatepsi --librarytype RNAseq --salmondir $add_to_path/$salmon_dir --sampconds $sampcond --conditionA $conditionA --conditionB $conditionB --gff $gff
        """
}
