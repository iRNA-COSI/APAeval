/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { CHECK_SAMPLESHEET;
          get_sample_info } from '../modules/check_samplesheet' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    CHECK_SAMPLESHEET ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
}
