// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MAKE_DIFFERENTIAL_TSV {
    tag "$samplesheet"
    publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode

    container "quay.io/biocontainers/python:3.8.3"

    input:
    path labrat_diff_results
    
    output:
    path "*"

    script:
    """
    make_differential_tsv.py $labrat_diff_results labrat_differential.tsv
    """
}
