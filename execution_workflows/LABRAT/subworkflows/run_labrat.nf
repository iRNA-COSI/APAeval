/*
 * RUN LABRAT
 */

include { LABRAT_MAKETFFASTA   } from '../modules/labrat_maketffasta' addParams( options: [:] )
include { SALMON_INDEX     } from '../modules/salmon_index'   addParams( options: [:] )
include { SALMON_QUANT     } from '../modules/salmon_quant'   addParams( options: [:] )

workflow RUN_LABRAT {
    take:
    ch_sample
    
    main:
    ch_sample
       .map { it -> [ it[5], it[6] ] }
       .unique()
       .set { ch_gff_fasta }

    /*
     * Run LABRAT PREPARE TF FASTA
     */
     ch_tffasta = ''
     LABRAT_MAKETFFASTA ( ch_gff_fasta )

    /*
     * Replace LABRAT RUNSALMON with SALMON
     */
     SALMON_INDEX ( LABRAT_MAKETFFASTA.out.ch_tffasta )
     ch_txfasta_idx = SALMON_INDEX.out.ch_txfasta_idx
     ch_txfasta_idx
       .combine ( ch_sample )
       .map { it -> [ it[0], it[1], it[2], it[3], it[9] ] }
       .set { ch_input_salmon_quant }
     SALMON_QUANT ( ch_input_salmon_quant )

     ch_salmon_quant_outputs = SALMON_QUANT.out.ch_salmon_quant_outputs
     ch_salmon_quant_outputs.view()
}

