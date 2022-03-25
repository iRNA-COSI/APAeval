/*
 * Convert BAM to BigWig
 */

params.qapa_options   = [:]

include { GTFTOGENEPRED  } from '../modules/gtf_to_genepred'     addParams( options: [:] )
include { GTFTOMARTEXPORT} from '../modules/gtf_to_martexport'   addParams( options: [:] )
include { QAPA_BUILD_REF } from '../modules/qapa_build_ref'   addParams( options: [:] )
include { QAPA_INDEX     } from '../modules/qapa_index'     addParams( options: params.qapa_options )
include { SALMON_INDEX   } from '../modules/salmon_index'   addParams( options: [:] )
include { SALMON_QUANT   } from '../modules/salmon_quant'   addParams( options: [:] )
include { QAPA_QUANT     } from '../modules/qapa_quant'     addParams( options: params.qapa_options )
include { MAKE_QUANT_BED } from '../modules/make_quant_bed' addParams( options: [:] )

workflow QAPA_SALMON {
    take:
    ch_sample
    ch_gtf
    ch_polyabed
    
    main:
    GTFTOGENEPRED ( ch_gtf )
    GTFTOMARTEXPORT ( ch_gtf )
    genepred = GTFTOGENEPRED.out.genepred
    martexport = GTFTOMARTEXPORT.out.mart_export
    QAPA_BUILD_REF ( martexport, genepred, ch_polyabed )

    ch_sample
       .map { it -> [ it[3], it[4] ] }
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
       .map { it -> [ it[0], it[1], it[2], it[3], it[6] ] }
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

