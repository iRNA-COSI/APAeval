/*
 * Check input samplesheet and get read channels
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../modules/functions'

def options    = initOptions(params.options)
params.options = [:]

process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/sample_info", mode: params.publish_dir_mode
    
    container "docker.io/apaeval/isoscm:1.0"
    
    input:
    path samplesheet
    
    output:
    path "*reformat.csv", emit: ch_samplesheet
    
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}


def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.strandinfo, sample.condition ]
}


workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    CHECK_SAMPLESHEET ( samplesheet )
    
    // Process per CSV record
    ch_samplesheet = CHECK_SAMPLESHEET.out.ch_samplesheet
    
    emit:
    ch_samplesheet
}
