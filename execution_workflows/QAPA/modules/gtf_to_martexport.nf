// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GTFTOMARTEXPORT {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path gtf
    
    output:
    path "*.txt", emit: mart_export
    
    script:
    """
    GTFtoMartExport.py $gtf
    """
}
