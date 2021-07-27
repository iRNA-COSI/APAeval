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
    
    input:
    path samplesheet
    
    output:
    path "*reformat.csv", emit: ch_samplesheet
    
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}


def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.bamdir, sample.strandinfo, sample.gtf, sample.condition ]
}


workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    CHECK_SAMPLESHEET ( samplesheet )
    
    // Process per CSV record
    ch_samplesheet = CHECK_SAMPLESHEET.out.ch_samplesheet
    ch_samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { it -> it.bamdir }
        .unique()
        .set { ch_bamdir }
    ch_bamdir.view()
    
    emit:
    ch_samplesheet
    ch_bamdir
}
