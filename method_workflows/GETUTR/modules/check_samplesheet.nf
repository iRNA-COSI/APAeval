// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode
    
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path "*reformat.csv", emit: ch_samplesheet
    
    script:
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}

def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.bam ]
}
