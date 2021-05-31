/*
 * Convert BAM to BigWig
 */

params.qapa_options   = [:]

include { QAPA_PREPARE_UTRLIB   } from '../modules/qapa_prepare_utrlib' addParams( options: params.qapa_options )
include { QAPA                  } from '../modules/qapa'                addParams( options: params.qapa_options )

workflow RUN_QAPA {
    take:
    ch_sample
    
    main:
    ch_sample
       .map { it -> [ it[6], it[7] ] }
       .unique()
       .set { ch_pre_utr_lib }

    /*
     * Run QAPA & SALMON PREPARE UTRLIB
     */
     QAPA_PREPARE_UTRLIB ( ch_pre_utr_lib )
     ch_utr_library = QAPA_PREPARE_UTRLIB.out.ch_utr_library

     ch_utr_library
       .combine ( ch_sample )
       .map { it -> [ it[0], it[1], it[2], it[3], it[9] ] }
       .set { ch_input_qapa }

    /*
     * Run QAPA & SALMON QUANT
     */
    QAPA ( ch_input_qapa )
    ch_qapa_output= QAPA.out.ch_qapa_outputs

    emit:
    ch_qapa_output

}

