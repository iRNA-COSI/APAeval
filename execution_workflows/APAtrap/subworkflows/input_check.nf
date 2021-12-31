/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { CHECK_SAMPLESHEET;
          get_sample_info } from '../modules/check_samplesheet' addParams( options: params.options )
include { CHECK_INPUT_PARAMS      } from '../modules/check_input_params' addParams( options: [:] )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    /*
        Check sample sheet
    */
    CHECK_SAMPLESHEET ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { ch_sample }

    /*
        Check input params
    */
    CHECK_INPUT_PARAMS( samplesheet )

    emit:
    ch_sample // [ sample,  ]
}
