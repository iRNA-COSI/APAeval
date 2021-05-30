/*
 * Convert BAM to BigWig
 */

params.qapa_options   = [:]

include { QAPA_INDEX     } from '../modules/qapa_index'     addParams( options: params.qapa_options )
include { SALMON_INDEX   } from '../modules/salmon_index'   addParams( options: [:] )
include { SALMON_QUANT   } from '../modules/salmon_quant'   addParams( options: [:] )
include { QAPA_QUANT     } from '../modules/qapa_quant'     addParams( options: params.qapa_options )
include { MAKE_QUANT_BED } from '../modules/make_quant_bed' addParams( options: [:] )

workflow QAPA_SALMON {
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
     QAPA_INDEX ( ch_pre_utr_lib )
     SALMON_INDEX( QAPA_INDEX.out.ch_indexed_fasta )
     ch_utr_library = SALMON_INDEX.out.ch_utr_library

     ch_utr_library
       .combine ( ch_sample )
       .map { it -> [ it[0], it[1], it[2], it[3], it[9] ] }
       .set { ch_input_salmon_quant }

    /*
     * Run QAPA & SALMON QUANT
     */
     SALMON_QUANT( ch_input_salmon_quant )
     QAPA_QUANT( SALMON_QUANT.out.ch_salmon_quant_outputs )
     ch_qapa_outputs= QAPA_QUANT.out.ch_qapa_outputs
     MAKE_QUANT_BED( ch_qapa_outputs )

//    emit:
//    ch_qapa_output

}

